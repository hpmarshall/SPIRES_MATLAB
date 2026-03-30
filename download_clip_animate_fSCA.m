%% download_clip_animate_fSCA.m
%
% Downloads SPIRES HIST V01 fractional snow covered area (fSCA) data from:
%     ftp://dtn.rc.colorado.edu/shares/snow-today/gridded_data/SPIRES_HIST_V01
%
% Clips the data to the Boise River Basin (BRB) using a user-supplied
% shapefile, and creates an MP4 animation of daily fSCA for a user-specified
% water year.  SNOTEL station locations within the basin are overlaid.
%
% The basin mask is cached to a .mat file so subsequent runs skip the
% expensive point-in-polygon computation.
%
% Requirements:
%     - MATLAB R2018b+ (for ftp, ncread, VideoWriter, shaperead)
%     - Mapping Toolbox (for shaperead; or use readgeotable in R2022+)
%
% Usage:
%     >> download_clip_animate_fSCA(2020, 'BRB_outline.shp')
%     >> download_clip_animate_fSCA(2015, 'BRB_outline.shp', 'fps', 15)
%     >> download_clip_animate_fSCA(2020, 'BRB_outline.shp', 'output_format', 'kmz')
%     >> download_clip_animate_fSCA(2020, 'BRB_outline.shp', 'output_format', 'both')
%
%   The year argument is the WATER YEAR:
%     WY 2020 = Oct 1 2019 – Sep 30 2020
%     WY 2015 = Oct 1 2014 – Sep 30 2015
%
%   output_format: 'mp4' (default video), 'kmz' (Google Earth overlay),
%                  or 'both'
%
% Data Citation:
%     Rittger, K., Lenard, S. J., Palomaki, R. T., Bair, E. H., Dozier, J. &
%     Mankoff, K. (2025). Historical MODIS/Terra L3 Global Daily 500m SIN Grid
%     Snow Cover, Snow Albedo, and Snow Surface Properties. (SPIRES_HIST, V1).
%     [Data Set]. Boulder, CO. NSIDC. https://doi.org/10.7265/a3vr-c014

function download_clip_animate_fSCA(water_year, shapefile_path, varargin)

    %% Parse optional arguments
    p = inputParser;
    addRequired(p, 'water_year', @(x) isnumeric(x) && x >= 2001 && x <= 2025);
    addRequired(p, 'shapefile_path', @ischar);
    addParameter(p, 'output_dir', './spires_output', @ischar);
    addParameter(p, 'fps', 10, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'tile', 'h09v04', @ischar);
    addParameter(p, 'output_format', 'both', ...
                 @(x) any(strcmpi(x, {'mp4','kmz','both'})));
    parse(p, water_year, shapefile_path, varargin{:});

    output_dir    = p.Results.output_dir;
    fps           = p.Results.fps;
    tile          = p.Results.tile;
    output_format = lower(p.Results.output_format);

    %% Water year date range
    % WY N = Oct 1 of year N-1  through  Sep 30 of year N
    wy_start = datenum(water_year - 1, 10, 1);
    wy_end   = datenum(water_year,      9, 30);
    cal_year1 = water_year - 1;   % calendar year containing Oct-Dec
    cal_year2 = water_year;       % calendar year containing Jan-Sep

    %% Constants
    FTP_HOST   = 'dtn.rc.colorado.edu';
    FTP_DIR    = '/shares/snow-today/gridded_data/SPIRES_HIST_V01';
    R_EARTH    = 6371007.181;      % MODIS sphere radius (m)
    TILE_SIZE  = 1111950.5197;     % tile edge length (m)
    NPIX       = 2400;             % pixels per tile edge
    PIXEL_SIZE = TILE_SIZE / NPIX;

    % Parse tile h and v
    h_tile = str2double(tile(2:3));
    v_tile = str2double(tile(5:6));

    %% ----------------------------------------------------------------
    %  SNOTEL stations within / near the Boise River Basin
    %  Coordinates from NRCS (deg-min converted to decimal degrees)
    %  ----------------------------------------------------------------
    snotel = struct( ...
        'name', {{'Bogus Basin', 'Mores Creek Summit', 'Graham Guard Sta.', ...
                  'Jackson Peak', 'Atlanta Summit', 'Trinity Mtn.', ...
                  'Banner Summit', 'Prairie'}}, ...
        'lat', [43+46/60, 43+56/60, 43+57/60, 44+3/60, ...
                43+45/60, 43+38/60, 44+18/60, 43+49/60], ...
        'lon', [-(116+6/60), -(115+40/60), -(115+16/60), -(115+27/60), ...
                -(115+14/60), -(115+26/60), -(115+14/60), -(115+54/60)], ...
        'elev_ft', [6370, 6090, 5680, 7060, 7570, 7790, 7040, 5580], ...
        'id',      [978,  637,  496,  550,  306,  830,  312,  700] ...
    );

    fprintf('==================================================================\n');
    fprintf('  SPIRES fSCA - Download, Clip, and Animate (MATLAB)\n');
    fprintf('==================================================================\n');
    fprintf('  Water Year: %d  (Oct 1 %d - Sep 30 %d)\n', ...
            water_year, cal_year1, cal_year2);
    fprintf('  Tile:       %s\n', tile);
    fprintf('  Shapefile:  %s\n', shapefile_path);
    fprintf('  Output:     %s\n', output_dir);
    fprintf('==================================================================\n\n');

    if ~exist(output_dir, 'dir'); mkdir(output_dir); end

    %% ====================================================================
    %  STEP 1: Read shapefile / load or create basin mask
    %  ====================================================================
    fprintf('[1/4] Basin mask...\n');

    [~, shp_stem] = fileparts(shapefile_path);
    mask_file = fullfile(output_dir, ...
                         sprintf('basin_mask_%s_%s.mat', shp_stem, tile));

    if exist(mask_file, 'file')
        fprintf('  Loading cached mask: %s\n', mask_file);
        tmp = load(mask_file);
        mask       = tmp.mask;
        row_min    = tmp.row_min;
        row_max    = tmp.row_max;
        col_min    = tmp.col_min;
        col_max    = tmp.col_max;
        nrows      = tmp.nrows;
        ncols      = tmp.ncols;
        all_poly_x_sin = tmp.all_poly_x_sin;
        all_poly_y_sin = tmp.all_poly_y_sin;
        utm_zone   = tmp.utm_zone;
        x_origin   = tmp.x_origin;
        y_origin   = tmp.y_origin;
        n_pixels = sum(mask(:));
        fprintf('  Basin contains %d pixels (%.1f km^2)\n', n_pixels, n_pixels*0.25);
    else
        fprintf('  Building mask from shapefile (will be cached for future runs)...\n');

        S = shaperead(shapefile_path);

        bbox = [S(1).BoundingBox];
        utm_zone = 11;
        if bbox(1,1) > 1000
            fprintf('  Detected projected coordinates (UTM Zone %dN assumed)\n', utm_zone);
            is_utm = true;
        else
            fprintf('  Detected geographic coordinates, will reproject to UTM Zone %dN\n', utm_zone);
            is_utm = false;
        end

        all_poly_x_sin = {};
        all_poly_y_sin = {};
        all_x_sin = [];
        all_y_sin = [];

        for k = 1:length(S)
            px = S(k).X;
            py = S(k).Y;
            nan_idx = [0, find(isnan(px)), length(px)+1];
            for p_idx = 1:(length(nan_idx)-1)
                i1 = nan_idx(p_idx) + 1;
                i2 = nan_idx(p_idx+1) - 1;
                if i2 < i1; continue; end
                part_x = px(i1:i2);
                part_y = py(i1:i2);
                if is_utm
                    [lat, lon] = utm2latlon(part_x, part_y, utm_zone);
                else
                    lat = part_y; lon = part_x;
                end
                [sx, sy] = latlon2sin(lat, lon, R_EARTH);
                all_poly_x_sin{end+1} = sx; %#ok<AGROW>
                all_poly_y_sin{end+1} = sy; %#ok<AGROW>
                all_x_sin = [all_x_sin; sx(:)]; %#ok<AGROW>
                all_y_sin = [all_y_sin; sy(:)]; %#ok<AGROW>
            end
        end

        bbox_sin = [min(all_x_sin), min(all_y_sin), ...
                    max(all_x_sin), max(all_y_sin)];

        x_origin = -R_EARTH * pi + h_tile * TILE_SIZE;
        y_origin =  R_EARTH * pi / 2 - v_tile * TILE_SIZE;

        col_min = max(1, floor((bbox_sin(1) - x_origin) / PIXEL_SIZE) - 1);
        col_max = min(NPIX, ceil((bbox_sin(3) - x_origin) / PIXEL_SIZE) + 2);
        row_min = max(1, floor((y_origin - bbox_sin(4)) / PIXEL_SIZE) - 1);
        row_max = min(NPIX, ceil((y_origin - bbox_sin(2)) / PIXEL_SIZE) + 2);
        nrows = row_max - row_min + 1;
        ncols = col_max - col_min + 1;

        fprintf('  Subset: rows [%d:%d], cols [%d:%d]\n', row_min, row_max, col_min, col_max);
        fprintf('  Mask dimensions: %d x %d\n', nrows, ncols);

        mask = false(nrows, ncols);
        for r = 1:nrows
            for c = 1:ncols
                mpx = x_origin + (col_min + c - 1 - 0.5) * PIXEL_SIZE;
                mpy = y_origin - (row_min + r - 1 - 0.5) * PIXEL_SIZE;
                for poly_idx = 1:length(all_poly_x_sin)
                    if inpolygon(mpx, mpy, all_poly_x_sin{poly_idx}, ...
                                          all_poly_y_sin{poly_idx})
                        mask(r, c) = true;
                        break;
                    end
                end
            end
        end

        n_pixels = sum(mask(:));
        fprintf('  Basin contains %d pixels (%.1f km^2)\n', n_pixels, n_pixels*0.25);

        fprintf('  Saving mask cache: %s\n', mask_file);
        save(mask_file, 'mask', 'row_min', 'row_max', 'col_min', 'col_max', ...
             'nrows', 'ncols', 'all_poly_x_sin', 'all_poly_y_sin', ...
             'utm_zone', 'x_origin', 'y_origin');
    end

    %% ====================================================================
    %  STEP 2: Download data (both calendar years spanning the water year)
    %  ====================================================================
    fprintf('\n[2/4] Downloading SPIRES data for WY%d...\n', water_year);

    download_dir = fullfile(output_dir, 'downloads');
    if ~exist(download_dir, 'dir'); mkdir(download_dir); end

    nc_files_y1 = discover_and_download(FTP_HOST, FTP_DIR, cal_year1, tile, download_dir);
    nc_files_y2 = discover_and_download(FTP_HOST, FTP_DIR, cal_year2, tile, download_dir);
    nc_files = unique([nc_files_y1, nc_files_y2]);

    if isempty(nc_files)
        error('No data files were downloaded. Check FTP connection and directory structure.');
    end

    %% ====================================================================
    %  STEP 3: Extract fSCA for the water year window
    %  ====================================================================
    fprintf('\n[3/4] Extracting fSCA data for WY%d...\n', water_year);

    [dates, nc_files_for_dates, nc_time_indices, fsca_varname] = ...
        extract_fsca_wy(nc_files, wy_start, wy_end, ...
                         row_min, row_max, col_min, col_max, ...
                         nrows, ncols);

    if isempty(dates)
        error('No fSCA data extracted. Check netCDF variable names.');
    end

    fprintf('  Found %d time steps\n', length(dates));
    fprintf('  Date range: %s to %s\n', ...
            datestr(dates(1), 'yyyy-mm-dd'), datestr(dates(end), 'yyyy-mm-dd'));

    %% ====================================================================
    %  STEP 4: Reproject sinusoidal → UTM and create outputs
    %  ====================================================================

    % --- Read _FillValue and valid_range from first netCDF file ---
    fill_val = NaN; valid_min_attr = 0; valid_max_attr = 1;
    try
        first_info = ncinfo(nc_files_for_dates{1});
        for vi = 1:length(first_info.Variables)
            if strcmp(first_info.Variables(vi).Name, fsca_varname)
                for aa = 1:length(first_info.Variables(vi).Attributes)
                    aname = first_info.Variables(vi).Attributes(aa).Name;
                    aval  = first_info.Variables(vi).Attributes(aa).Value;
                    if strcmp(aname,'_FillValue');    fill_val = double(aval); end
                    if strcmp(aname,'missing_value'); fill_val = double(aval); end
                    if strcmp(aname,'valid_min');     valid_min_attr = double(aval); end
                    if strcmp(aname,'valid_max');     valid_max_attr = double(aval); end
                    if strcmp(aname,'valid_range')
                        valid_min_attr = double(aval(1));
                        valid_max_attr = double(aval(2));
                    end
                end
                break;
            end
        end
    catch; end
    is_pct = valid_max_attr > 1.5;
    fprintf('  Data attributes: _FillValue=%g, valid_range=[%g, %g]\n', ...
            fill_val, valid_min_attr, valid_max_attr);

    % --- Sinusoidal pixel center coordinates (for the mask subset) ---
    sin_x_vec = x_origin + ((col_min:col_max) - 0.5) * PIXEL_SIZE;
    sin_y_vec = y_origin - ((row_min:row_max) - 0.5) * PIXEL_SIZE;

    % --- Compute lat/lon for EVERY basin pixel to get true UTM extent ---
    fprintf('\n[4/5] Reprojecting sinusoidal data to UTM...\n');
    [sin_x_grid, sin_y_grid] = meshgrid(sin_x_vec, sin_y_vec);
    lat_grid = rad2deg(sin_y_grid ./ R_EARTH);
    cos_lat  = cos(sin_y_grid ./ R_EARTH);
    cos_lat(abs(cos_lat) < 1e-10) = 1e-10;
    lon_grid = rad2deg(sin_x_grid ./ (R_EARTH .* cos_lat));
    [utm_e_grid, utm_n_grid] = latlon2utm(lat_grid, lon_grid, utm_zone);

    % Determine UTM bounding box from actual basin pixel locations
    utm_e_valid = utm_e_grid(mask);
    utm_n_valid = utm_n_grid(mask);
    utm_e_min = min(utm_e_valid) - PIXEL_SIZE;
    utm_e_max = max(utm_e_valid) + PIXEL_SIZE;
    utm_n_min = min(utm_n_valid) - PIXEL_SIZE;
    utm_n_max = max(utm_n_valid) + PIXEL_SIZE;

    % Build regular UTM output grid (~500 m resolution)
    utm_res  = 500;
    utm_e_out = utm_e_min : utm_res : utm_e_max;
    utm_n_out = utm_n_max : -utm_res : utm_n_min;
    [UTM_E_OUT, UTM_N_OUT] = meshgrid(utm_e_out, utm_n_out);
    nr_out = length(utm_n_out);
    nc_out = length(utm_e_out);

    % --- Build reprojection lookup using FULL TILE indexing ---
    %     The sinusoidal projection is nonlinear: UTM cells in the north
    %     of the basin may need source pixels from columns beyond our
    %     original col_min:col_max subset. We index into the FULL TILE
    %     and expand the data read range accordingly.
    fprintf('  Building reprojection lookup table (%d x %d output grid)...\n', nr_out, nc_out);
    [lat_out, lon_out] = utm2latlon(UTM_E_OUT(:), UTM_N_OUT(:), utm_zone);
    [sin_x_out, sin_y_out] = latlon2sin(lat_out, lon_out, R_EARTH);

    % Map to FULL TILE row/col (1-based, 1 to NPIX)
    tile_col_idx = round((sin_x_out - x_origin) / PIXEL_SIZE + 0.5);
    tile_row_idx = round((y_origin - sin_y_out) / PIXEL_SIZE + 0.5);

    % Determine expanded read range needed
    valid_tile = tile_col_idx >= 1 & tile_col_idx <= NPIX & ...
                 tile_row_idx >= 1 & tile_row_idx <= NPIX;
    read_col_min = max(1,    min(tile_col_idx(valid_tile)));
    read_col_max = min(NPIX, max(tile_col_idx(valid_tile)));
    read_row_min = max(1,    min(tile_row_idx(valid_tile)));
    read_row_max = min(NPIX, max(tile_row_idx(valid_tile)));
    read_nrows = read_row_max - read_row_min + 1;
    read_ncols = read_col_max - read_col_min + 1;

    fprintf('  Original subset: rows [%d:%d] cols [%d:%d] (%dx%d)\n', ...
            row_min, row_max, col_min, col_max, nrows, ncols);
    fprintf('  Expanded read:   rows [%d:%d] cols [%d:%d] (%dx%d)\n', ...
            read_row_min, read_row_max, read_col_min, read_col_max, ...
            read_nrows, read_ncols);

    % Convert tile indices to local (expanded-read) indices
    local_col_idx = tile_col_idx - read_col_min + 1;
    local_row_idx = tile_row_idx - read_row_min + 1;
    valid_lookup = valid_tile & ...
                   local_col_idx >= 1 & local_col_idx <= read_ncols & ...
                   local_row_idx >= 1 & local_row_idx <= read_nrows;

    local_col_idx = reshape(local_col_idx, nr_out, nc_out);
    local_row_idx = reshape(local_row_idx, nr_out, nc_out);
    valid_lookup  = reshape(valid_lookup,  nr_out, nc_out);

    % Also need a mask for the expanded read range
    % (the original mask covers col_min:col_max, row_min:row_max;
    %  we embed it into the expanded range)
    mask_expanded = false(read_nrows, read_ncols);
    r_off = row_min - read_row_min;
    c_off = col_min - read_col_min;
    mask_expanded(r_off+1 : r_off+nrows, c_off+1 : c_off+ncols) = mask;

    % Build linear index for fast lookup
    src_lin_idx = zeros(nr_out, nc_out);
    for r = 1:nr_out
        for c = 1:nc_out
            if valid_lookup(r,c)
                lr = local_row_idx(r,c);
                lc = local_col_idx(r,c);
                if mask_expanded(lr, lc)
                    src_lin_idx(r,c) = sub2ind([read_nrows, read_ncols], lr, lc);
                end
            end
        end
    end
    mask_out = src_lin_idx > 0;

    % --- Re-read data with expanded range and reproject ---
    fprintf('  Re-reading data with expanded range and reprojecting %d frames...\n', length(dates));
    data_utm = NaN(nr_out, nc_out, length(dates));

    for t = 1:length(dates)
        % Re-read this frame with the expanded spatial range
        nc_path_t = nc_files_for_dates{t};
        nc_tidx_t = nc_time_indices(t);

        % Read the expanded subset
        fr_exp = read_one_frame(nc_path_t, fsca_varname, nc_tidx_t, ...
                                read_row_min, read_col_min, read_nrows, read_ncols);

        % Apply mask, fill, valid range
        fr_exp(~mask_expanded) = NaN;
        if isfinite(fill_val); fr_exp(fr_exp == fill_val) = NaN; end
        fr_exp(fr_exp < valid_min_attr | fr_exp > valid_max_attr) = NaN;

        % Detect unobserved pixels
        try
            gs = read_one_frame(nc_path_t, 'grain_size', nc_tidx_t, ...
                                read_row_min, read_col_min, read_nrows, read_ncols);
            fr_exp((fr_exp == 0) & (gs == 0)) = NaN;
        catch; end

        if is_pct; fr_exp = fr_exp / 100.0; end

        % Reproject
        fr_out = NaN(nr_out, nc_out);
        fr_out(mask_out) = fr_exp(src_lin_idx(mask_out));
        data_utm(:,:,t) = fr_out;
    end

    % --- Lat/lon bounding box (for KMZ) ---
    lat_north = max(lat_out);  lat_south = min(lat_out);
    lon_west  = min(lon_out);  lon_east  = max(lon_out);

    % --- Basin polygons in UTM ---
    poly_utm_e = {};  poly_utm_n = {};
    for poly_idx = 1:length(all_poly_x_sin)
        [plat, plon] = sin2latlon(all_poly_x_sin{poly_idx}, ...
                                   all_poly_y_sin{poly_idx}, R_EARTH);
        [pe, pn] = latlon2utm(plat, plon, utm_zone);
        poly_utm_e{end+1} = pe; %#ok<AGROW>
        poly_utm_n{end+1} = pn; %#ok<AGROW>
    end

    % --- SNOTEL in UTM ---
    [snotel_e, snotel_n] = latlon2utm(snotel.lat, snotel.lon, utm_zone);

    % ---- Generate MP4 ----
    if any(strcmp(output_format, {'mp4','both'}))
        fprintf('\n[5a] Creating MP4 animation (UTM Zone %dN)...\n', utm_zone);
        output_mp4 = fullfile(output_dir, sprintf('BRB_fSCA_WY%d.mp4', water_year));
        create_movie(dates, data_utm, mask_out, output_mp4, water_year, fps, ...
                     utm_e_min, utm_e_max, utm_n_min, utm_n_max, ...
                     poly_utm_e, poly_utm_n, utm_zone, ...
                     snotel, snotel_e, snotel_n);
        fprintf('  MP4 saved: %s\n', output_mp4);
    end

    % ---- Generate KMZ for Google Earth ----
    if any(strcmp(output_format, {'kmz','both'}))
        fprintf('\n[5b] Creating KMZ for Google Earth...\n');
        output_kmz = fullfile(output_dir, sprintf('BRB_fSCA_WY%d.kmz', water_year));
        create_kmz(dates, data_utm, mask_out, output_kmz, water_year, ...
                   lat_north, lat_south, lon_west, lon_east, ...
                   snotel);
        fprintf('  KMZ saved: %s\n', output_kmz);
    end

    fprintf('\n==================================================================\n');
    fprintf('  COMPLETE — outputs in: %s\n', output_dir);
    fprintf('==================================================================\n');
end


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [lat, lon] = utm2latlon(easting, northing, zone)
    a  = 6378137.0;  f = 1/298.257223563;
    e  = sqrt(2*f - f^2);  ep = e/sqrt(1-e^2);  k0 = 0.9996;
    x = easting - 500000.0;  y = northing;
    M = y/k0;
    mu = M / (a*(1 - e^2/4 - 3*e^4/64 - 5*e^6/256));
    e1 = (1-sqrt(1-e^2))/(1+sqrt(1-e^2));
    phi1 = mu + (3*e1/2-27*e1^3/32).*sin(2*mu) ...
              + (21*e1^2/16-55*e1^4/32).*sin(4*mu) ...
              + (151*e1^3/96).*sin(6*mu);
    N1 = a./sqrt(1-e^2*sin(phi1).^2);
    T1 = tan(phi1).^2;  C1 = ep^2*cos(phi1).^2;
    R1 = a*(1-e^2)./(1-e^2*sin(phi1).^2).^1.5;
    D  = x./(N1*k0);
    lat_rad = phi1 - (N1.*tan(phi1)./R1).*( ...
        D.^2/2 - (5+3*T1+10*C1-4*C1.^2-9*ep^2).*D.^4/24 ...
        + (61+90*T1+298*C1+45*T1.^2-252*ep^2-3*C1.^2).*D.^6/720);
    lon_rad = (D-(1+2*T1+C1).*D.^3/6 ...
        + (5-2*C1+28*T1-3*C1.^2+8*ep^2+24*T1.^2).*D.^5/120)./cos(phi1);
    lon0 = deg2rad((zone-1)*6-180+3);
    lat = rad2deg(lat_rad);  lon = rad2deg(lon_rad+lon0);
end

function [easting, northing] = latlon2utm(lat, lon, zone)
    a = 6378137.0; f = 1/298.257223563;
    e = sqrt(2*f-f^2); ep = e/sqrt(1-e^2); k0 = 0.9996;
    lon0 = (zone-1)*6-180+3;
    lat_rad = deg2rad(lat); lon_rad = deg2rad(lon); lon0_rad = deg2rad(lon0);
    N = a./sqrt(1-e^2*sin(lat_rad).^2);
    T = tan(lat_rad).^2; C = ep^2*cos(lat_rad).^2;
    A = (lon_rad-lon0_rad).*cos(lat_rad);
    M = a*((1-e^2/4-3*e^4/64-5*e^6/256).*lat_rad ...
          -(3*e^2/8+3*e^4/32+45*e^6/1024).*sin(2*lat_rad) ...
          +(15*e^4/256+45*e^6/1024).*sin(4*lat_rad) ...
          -(35*e^6/3072).*sin(6*lat_rad));
    easting  = k0*N.*(A+(1-T+C).*A.^3/6 ...
               +(5-18*T+T.^2+72*C-58*ep^2).*A.^5/120)+500000.0;
    northing = k0*(M+N.*tan(lat_rad).*( ...
               A.^2/2+(5-T+9*C+4*C.^2).*A.^4/24 ...
               +(61-58*T+T.^2+600*C-330*ep^2).*A.^6/720));
end

function [x, y] = latlon2sin(lat, lon, R)
    lat_rad = deg2rad(lat); lon_rad = deg2rad(lon);
    x = R.*lon_rad.*cos(lat_rad); y = R.*lat_rad;
end

function [lat, lon] = sin2latlon(x, y, R)
    lat_rad = y./R; lat = rad2deg(lat_rad);
    cos_lat = cos(lat_rad); cos_lat(abs(cos_lat)<1e-10) = 1e-10;
    lon = rad2deg(x./(R.*cos_lat));
end


function nc_files = discover_and_download(ftp_host, ftp_dir, year, tile, download_dir)
    nc_files = {};
    try
        f = ftp(ftp_host);
    catch ME
        warning('FTP connection failed: %s', ME.message);
        local = dir(fullfile(download_dir, '*.nc'));
        for i = 1:length(local)
            nc_files{end+1} = fullfile(download_dir, local(i).name); %#ok<AGROW>
        end
        return;
    end
    try
        cd(f, ftp_dir);
        top_items = dir(f);
        top_names = {top_items.name}; top_isdir = [top_items.isdir];
        year_str = num2str(year); downloaded = false;

        % Strategy 1: top-level files
        for i = 1:length(top_names)
            name = top_names{i};
            if ~top_isdir(i) && contains(name,tile) && contains(name,year_str) ...
                    && endsWith(name,'.nc','IgnoreCase',true)
                lp = fullfile(download_dir,name);
                if ~exist(lp,'file'); fprintf('  Downloading: %s\n',name); mget(f,name,download_dir);
                else; fprintf('  Already downloaded: %s\n',name); end
                nc_files{end+1}=lp; downloaded=true; %#ok<AGROW>
            end
        end
        if downloaded; close(f); return; end

        % Strategy 2: tile subdir
        if any(strcmp(top_names,tile)&top_isdir)
            cd(f,tile); sub=dir(f); sn={sub.name}; sd=[sub.isdir];
            for i=1:length(sn)
                if ~sd(i)&&contains(sn{i},year_str)&&endsWith(sn{i},'.nc','IgnoreCase',true)
                    lp=fullfile(download_dir,sn{i});
                    if ~exist(lp,'file'); fprintf('  Downloading: %s/%s\n',tile,sn{i}); mget(f,sn{i},download_dir); end
                    nc_files{end+1}=lp; downloaded=true; %#ok<AGROW>
                end
            end
            if downloaded; close(f); return; end
            if any(strcmp(sn,year_str))
                cd(f,year_str); yi=dir(f); yn={yi.name}; yd=[yi.isdir];
                for i=1:length(yn)
                    if ~yd(i)&&endsWith(yn{i},'.nc','IgnoreCase',true)
                        lp=fullfile(download_dir,yn{i});
                        if ~exist(lp,'file'); mget(f,yn{i},download_dir); end
                        nc_files{end+1}=lp; %#ok<AGROW>
                    end
                end
            end
            cd(f,ftp_dir);
        end
        if ~isempty(nc_files); close(f); return; end

        % Strategy 3: year subdir
        if any(strcmp(top_names,year_str)&top_isdir)
            cd(f,year_str); sub=dir(f); sn={sub.name}; sd=[sub.isdir];
            for i=1:length(sn)
                if ~sd(i)&&contains(sn{i},tile)&&endsWith(sn{i},'.nc','IgnoreCase',true)
                    lp=fullfile(download_dir,sn{i});
                    if ~exist(lp,'file'); mget(f,sn{i},download_dir); end
                    nc_files{end+1}=lp; %#ok<AGROW>
                end
            end
            if ~isempty(nc_files); close(f); return; end
            if any(strcmp(sn,tile))
                cd(f,tile); ti=dir(f); tn={ti.name}; td=[ti.isdir];
                for i=1:length(tn)
                    if ~td(i)&&endsWith(tn{i},'.nc','IgnoreCase',true)
                        lp=fullfile(download_dir,tn{i});
                        if ~exist(lp,'file'); mget(f,tn{i},download_dir); end
                        nc_files{end+1}=lp; %#ok<AGROW>
                    end
                end
            end
        end

        % Strategy 4: any tile .nc at top
        if isempty(nc_files)
            cd(f,ftp_dir);
            for i=1:length(top_names)
                name=top_names{i};
                if ~top_isdir(i)&&contains(name,tile)&&endsWith(name,'.nc','IgnoreCase',true)
                    lp=fullfile(download_dir,name);
                    if ~exist(lp,'file'); mget(f,name,download_dir); end
                    nc_files{end+1}=lp; %#ok<AGROW>
                end
            end
        end
        close(f);
    catch ME
        try close(f); catch; end
        warning('FTP error: %s', ME.message);
        local=dir(fullfile(download_dir,'*.nc'));
        for i=1:length(local); nc_files{end+1}=fullfile(download_dir,local(i).name); end %#ok<AGROW>
    end
    if isempty(nc_files)
        fprintf('  WARNING: No files found for %s/%d.\n', tile, year);
    end
end


function [dates, nc_paths_out, time_indices_out, fsca_varname] = ...
        extract_fsca_wy(nc_files, wy_start, wy_end, ...
                         row_min, row_max, col_min, col_max, nrows, ncols)
    % Scan netCDF files and return dates, file paths, and time indices
    % for all frames within the water year window.
    dates = []; nc_paths_out = {}; time_indices_out = [];
    fsca_varname = '';

    for f_idx = 1:length(nc_files)
        nc_path = nc_files{f_idx}; [~,fname] = fileparts(nc_path);
        fprintf('  Scanning: %s\n', fname);
        try info = ncinfo(nc_path); catch ME
            fprintf('    Skipping: %s\n', ME.message); continue; end
        var_names = {info.Variables.Name};

        % Find fSCA variable
        fvar = '';
        if any(strcmp(var_names,'snow_fraction')); fvar = 'snow_fraction'; end
        if isempty(fvar)
            other_exact = {'fSCA','fsca','FSCA','snow_fraction_on_ground', ...
                           'Snow_Fraction','fractional_snow_cover','snow_cover_fraction'};
            for i=1:length(other_exact)
                if any(strcmp(var_names,other_exact{i})); fvar=other_exact{i}; break; end
            end
        end
        if isempty(fvar) && any(strcmp(var_names,'viewable_snow_fraction'))
            fvar = 'viewable_snow_fraction';
        end
        if isempty(fvar)
            for i=1:length(var_names)
                if contains(lower(var_names{i}),'snow')&&contains(lower(var_names{i}),'frac')
                    fvar=var_names{i}; break; end
            end
        end
        if isempty(fvar)
            fprintf('    WARNING: No fSCA variable found.\n'); continue;
        end
        if isempty(fsca_varname); fsca_varname = fvar; end

        var_info = [];
        for i=1:length(info.Variables)
            if strcmp(info.Variables(i).Name,fvar); var_info=info.Variables(i); break; end
        end
        ndims_var = length(var_info.Size);

        if ndims_var == 3
            dim_names = {var_info.Dimensions.Name};
            dim_sizes = [var_info.Dimensions.Length];
            time_dim = 0;
            for d=1:3
                if contains(lower(dim_names{d}),'time')||contains(lower(dim_names{d}),'day')
                    time_dim=d; break; end
            end
            if time_dim==0; time_dim=3; end

            time_var_name = '';
            for tname={'time','Time','date','day'}
                if any(strcmp(var_names,tname{1})); time_var_name=tname{1}; break; end
            end
            nt = dim_sizes(time_dim);

            for t=1:nt
                if ~isempty(time_var_name)
                    try
                        tv=ncread(nc_path,time_var_name,t,1);
                        ti=ncinfo(nc_path,time_var_name); tu='';
                        for a=1:length(ti.Attributes)
                            if strcmp(ti.Attributes(a).Name,'units'); tu=ti.Attributes(a).Value; break; end
                        end
                        if contains(tu,'days since')
                            pts=strsplit(tu); dt=datenum(pts{3})+double(tv);
                        elseif contains(tu,'seconds since')
                            pts=strsplit(tu); dt=datenum(pts{3})+double(tv)/86400;
                        elseif tv>100000; dt=double(tv);
                        elseif tv>2000000; dt=datenum(floor(tv/1000),1,mod(tv,1000));
                        else; dt=wy_start+t-1;
                        end
                    catch; dt=wy_start+t-1;
                    end
                else; dt=wy_start+t-1;
                end

                if dt>=wy_start && dt<=wy_end
                    dates = [dates; dt]; %#ok<AGROW>
                    nc_paths_out{end+1} = nc_path; %#ok<AGROW>
                    time_indices_out = [time_indices_out; t]; %#ok<AGROW>
                end
            end

        elseif ndims_var == 2 || (ndims_var == 3 && min([var_info.Dimensions.Length]) == 1)
            % Single-day file or 3D with singleton time
            [~,fn] = fileparts(nc_path);
            tok = regexp(fn,'(20\d{6})','tokens');
            if ~isempty(tok)
                dt = datenum(tok{1}{1},'yyyymmdd');
            else
                tok = regexp(fn,'(\d{7})','tokens');
                if ~isempty(tok)
                    dt = datenum(str2double(tok{1}{1}(1:4)),1,str2double(tok{1}{1}(5:7)));
                else; dt = wy_start + length(dates);
                end
            end
            if dt>=wy_start && dt<=wy_end
                dates = [dates; dt]; %#ok<AGROW>
                nc_paths_out{end+1} = nc_path; %#ok<AGROW>
                time_indices_out = [time_indices_out; 1]; %#ok<AGROW>
            end
        end
    end

    % Sort by date
    if ~isempty(dates)
        [dates, si] = sort(dates);
        nc_paths_out = nc_paths_out(si);
        time_indices_out = time_indices_out(si);
    end
end


function create_movie(dates, data, mask, output_path, water_year, fps, ...
                      utm_e_min, utm_e_max, utm_n_min, utm_n_max, ...
                      poly_utm_e, poly_utm_n, utm_zone, ...
                      snotel, snotel_e, snotel_n)
    if ~exist(fileparts(output_path),'dir'); mkdir(fileparts(output_path)); end
    nt=length(dates);
    fprintf('  Creating animation: %d frames at %d fps\n', nt, fps);

    % Build plot_context struct for plot_fsca_frame (Mode 2)
    buf = 2000;
    ctx.ek         = [utm_e_min utm_e_max] / 1000;
    ctx.nk         = [utm_n_min utm_n_max] / 1000;
    ctx.poly_utm_e = poly_utm_e;
    ctx.poly_utm_n = poly_utm_n;
    ctx.snotel_e   = snotel_e;
    ctx.snotel_n   = snotel_n;
    ctx.in_view    = snotel_e >= (utm_e_min-buf) & snotel_e <= (utm_e_max+buf) & ...
                     snotel_n >= (utm_n_min-buf) & snotel_n <= (utm_n_max+buf);
    ctx.utm_zone   = utm_zone;

    v = VideoWriter(output_path, 'MPEG-4');
    v.FrameRate = fps; v.Quality = 90; open(v);
    fig = figure('Position', [100 100 1000 800], 'Visible', 'off');

    for t = 1:nt
        if mod(t,30)==1 || t==nt
            fprintf('  Frame %d/%d: %s\n', t, nt, datestr(dates(t),'yyyy-mm-dd'));
        end

        plot_fsca_frame('frame_data',   data(:,:,t), ...
                        'frame_date',   dates(t), ...
                        'plot_context', ctx, ...
                        'water_year',   water_year, ...
                        'fig_handle',   fig, ...
                        'visible',      'off');

        writeVideo(v, getframe(fig));
    end
    close(v); close(fig);
    fprintf('  Animation saved: %s\n', output_path);
end


function create_kmz(dates, data, mask, output_path, water_year, ...
                    lat_north, lat_south, lon_west, lon_east, snotel)
    % Create a KMZ file with time-stamped ground overlays for Google Earth.
    %
    % Each day becomes a transparent PNG georeferenced via a GroundOverlay
    % with a TimeSpan.  SNOTEL sites are included as Placemarks.  The KML
    % and PNGs are bundled into a single .kmz (zipped KML).

    [out_dir, out_stem] = fileparts(output_path);
    tmp_dir = fullfile(out_dir, [out_stem '_tmp']);
    img_dir = fullfile(tmp_dir, 'images');
    if exist(tmp_dir, 'dir'); rmdir(tmp_dir, 's'); end
    mkdir(img_dir);

    nt = length(dates);
    fprintf('  Rendering %d overlay images...\n', nt);

    % Build an RGBA colormap: 0 = transparent, then blue ramp for fSCA
    ncolors = 256;
    cmap_rgb = parula(ncolors);

    % Render each frame as a transparent PNG
    img_names = cell(nt, 1);
    for t = 1:nt
        if mod(t, 30) == 1 || t == nt
            fprintf('    Image %d/%d: %s\n', t, nt, datestr(dates(t), 'yyyy-mm-dd'));
        end

        frame = data(:, :, t);

        % Map fSCA [0,1] to color index [1, ncolors]
        cidx = round(frame * (ncolors - 1)) + 1;
        cidx(cidx < 1) = 1;
        cidx(cidx > ncolors) = ncolors;

        % Build RGBA image
        [nr, nc] = size(frame);
        rgba = zeros(nr, nc, 4, 'uint8');

        for r = 1:nr
            for c = 1:nc
                if mask(r,c) && isfinite(frame(r,c))
                    idx = cidx(r,c);
                    rgba(r,c,1) = uint8(cmap_rgb(idx,1) * 255);
                    rgba(r,c,2) = uint8(cmap_rgb(idx,2) * 255);
                    rgba(r,c,3) = uint8(cmap_rgb(idx,3) * 255);
                    if frame(r,c) < 0.01
                        rgba(r,c,4) = 0;    % essentially snow-free = transparent
                    else
                        rgba(r,c,4) = 200;  % semi-transparent overlay
                    end
                end
                % else alpha stays 0 (fully transparent outside basin)
            end
        end

        fname = sprintf('fsca_%s.png', datestr(dates(t), 'yyyymmdd'));
        imwrite(rgba(:,:,1:3), fullfile(img_dir, fname), ...
                'Alpha', rgba(:,:,4));
        img_names{t} = fname;
    end

    % --- Write KML ---
    fprintf('  Writing KML...\n');
    kml_path = fullfile(tmp_dir, 'doc.kml');
    fid = fopen(kml_path, 'w');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n');
    fprintf(fid, '<Document>\n');
    fprintf(fid, '  <name>Boise River Basin fSCA — WY %d</name>\n', water_year);
    fprintf(fid, '  <description>Daily fractional snow covered area from SPIReS HIST V01.\nData: Rittger et al. (2025), NSIDC, doi:10.7265/a3vr-c014</description>\n');

    % Color scale legend as a ScreenOverlay
    fprintf(fid, '  <ScreenOverlay>\n');
    fprintf(fid, '    <name>fSCA Color Scale</name>\n');
    fprintf(fid, '    <description>0 (transparent) to 1 (yellow/bright)</description>\n');
    fprintf(fid, '    <visibility>1</visibility>\n');
    fprintf(fid, '  </ScreenOverlay>\n');

    % GroundOverlay for each day
    for t = 1:nt
        dt_str  = datestr(dates(t), 'yyyy-mm-dd');
        if t < nt
            dt_end = datestr(dates(t+1), 'yyyy-mm-dd');
        else
            dt_end = datestr(dates(t) + 1, 'yyyy-mm-dd');
        end

        fprintf(fid, '  <GroundOverlay>\n');
        fprintf(fid, '    <name>fSCA %s</name>\n', dt_str);
        fprintf(fid, '    <TimeSpan>\n');
        fprintf(fid, '      <begin>%sT00:00:00Z</begin>\n', dt_str);
        fprintf(fid, '      <end>%sT00:00:00Z</end>\n', dt_end);
        fprintf(fid, '    </TimeSpan>\n');
        fprintf(fid, '    <Icon><href>images/%s</href></Icon>\n', img_names{t});
        fprintf(fid, '    <LatLonBox>\n');
        fprintf(fid, '      <north>%.6f</north>\n', lat_north);
        fprintf(fid, '      <south>%.6f</south>\n', lat_south);
        fprintf(fid, '      <east>%.6f</east>\n',   lon_east);
        fprintf(fid, '      <west>%.6f</west>\n',   lon_west);
        fprintf(fid, '    </LatLonBox>\n');
        fprintf(fid, '  </GroundOverlay>\n');
    end

    % SNOTEL Placemarks
    fprintf(fid, '  <Folder>\n');
    fprintf(fid, '    <name>SNOTEL Stations</name>\n');
    fprintf(fid, '    <Style id="snotelStyle">\n');
    fprintf(fid, '      <IconStyle>\n');
    fprintf(fid, '        <color>ff0000ff</color>\n');
    fprintf(fid, '        <scale>0.8</scale>\n');
    fprintf(fid, '        <Icon><href>http://maps.google.com/mapfiles/kml/shapes/triangle.png</href></Icon>\n');
    fprintf(fid, '      </IconStyle>\n');
    fprintf(fid, '      <LabelStyle><scale>0.7</scale></LabelStyle>\n');
    fprintf(fid, '    </Style>\n');

    for si = 1:length(snotel.name)
        fprintf(fid, '    <Placemark>\n');
        fprintf(fid, '      <name>%s (#%d)</name>\n', snotel.name{si}, snotel.id(si));
        fprintf(fid, '      <description>Elevation: %d ft</description>\n', snotel.elev_ft(si));
        fprintf(fid, '      <styleUrl>#snotelStyle</styleUrl>\n');
        fprintf(fid, '      <Point><coordinates>%.6f,%.6f,0</coordinates></Point>\n', ...
                snotel.lon(si), snotel.lat(si));
        fprintf(fid, '    </Placemark>\n');
    end
    fprintf(fid, '  </Folder>\n');

    fprintf(fid, '</Document>\n');
    fprintf(fid, '</kml>\n');
    fclose(fid);

    % --- Zip into KMZ ---
    fprintf('  Packaging KMZ...\n');
    zip_path = [output_path '.zip'];
    zip(zip_path, {'doc.kml', 'images'}, tmp_dir);

    % Rename .zip to .kmz
    if exist(output_path, 'file'); delete(output_path); end
    movefile(zip_path, output_path);

    % Clean up temp directory
    rmdir(tmp_dir, 's');
    fprintf('  KMZ created: %s\n', output_path);
end


function frame = read_one_frame(nc_path, varname, time_idx, ...
                                 read_row_min, read_col_min, read_nrows, read_ncols)
    % Read one 2D frame from a netCDF file, handling both 2D and 3D variables.
    info_v = ncinfo(nc_path, varname);
    ndv = length(info_v.Size);

    if ndv == 3
        % Determine time dimension
        dim_names = {info_v.Dimensions.Name};
        time_dim = 0;
        for d = 1:3
            if contains(lower(dim_names{d}), 'time') || contains(lower(dim_names{d}), 'day')
                time_dim = d; break;
            end
        end
        if time_dim == 0; time_dim = 3; end
        spatial_dims = setdiff(1:3, time_dim);

        s = ones(1,3); c = [info_v.Dimensions.Length];
        s(time_dim) = time_idx; c(time_dim) = 1;
        s(spatial_dims(1)) = read_row_min; c(spatial_dims(1)) = read_nrows;
        s(spatial_dims(2)) = read_col_min; c(spatial_dims(2)) = read_ncols;
        frame = double(squeeze(ncread(nc_path, varname, s, c)));
    else
        frame = double(ncread(nc_path, varname, ...
                       [read_row_min, read_col_min], [read_nrows, read_ncols]));
    end

    if size(frame, 1) ~= read_nrows
        frame = frame';
    end
end
