function [fig, ax] = plot_fsca_frame(varargin)
%PLOT_FSCA_FRAME  Plot a single day of fractional snow covered area.
%
%  This function can be used in two ways:
%
%  MODE 1 — Quick look from a single daily netCDF file:
%     plot_fsca_frame('nc_file', 'SPIRES_HIST_h09v04_MOD09GA061_20200315_V1.0.nc', ...
%                     'shapefile', 'BRB_outline.shp')
%
%     % With a cached mask (much faster after first run):
%     plot_fsca_frame('nc_file', 'SPIRES_HIST_h09v04_MOD09GA061_20200315_V1.0.nc', ...
%                     'mask_file', './spires_output/basin_mask_BRB_outline_h09v04.mat')
%
%     % For a multi-day file, specify which time index to plot:
%     plot_fsca_frame('nc_file', 'SPIRES_HIST_h09v04_2020.nc', ...
%                     'day_index', 45, ...
%                     'shapefile', 'BRB_outline.shp')
%
%     % Save to PNG:
%     plot_fsca_frame('nc_file', 'SPIRES_HIST_h09v04_MOD09GA061_20200315_V1.0.nc', ...
%                     'shapefile', 'BRB_outline.shp', ...
%                     'save_png', 'test_frame.png')
%
%  MODE 2 — Called with pre-processed data (used by the animation script):
%     plot_fsca_frame('frame_data', data_utm(:,:,t), ...
%                     'frame_date', dates(t), ...
%                     'plot_context', ctx)
%
%  Optional parameters (both modes):
%     'fig_handle'    — existing figure handle to draw into (default: new)
%     'visible'       — 'on' or 'off' (default: 'on' in Mode 1, 'off' in 2)
%     'save_png'      — file path to save the frame as a PNG (default: '')
%     'water_year'    — water year label for the title (default: auto)
%
%  The function returns [fig, ax] handles so you can further customize.
%
%  WORKFLOW:  Edit the APPEARANCE section below until a single frame looks
%  perfect, then run download_clip_animate_fSCA.m to produce the movie.

    %% ================================================================
    %  APPEARANCE — edit these to taste
    %  ================================================================
    COLORMAP       = [0.96 0.96 0.86; parula(256)];  % first row = NaN color
    FSCA_RANGE     = [0 1];            % colorbar limits
    BASIN_COLOR    = 'k';              % basin boundary line color
    BASIN_WIDTH    = 1.2;              % basin boundary line width
    SNOTEL_MARKER  = '^';              % marker shape
    SNOTEL_SIZE    = 8;                % marker size
    SNOTEL_FACE    = [1 0.2 0.2];      % marker fill color (red)
    SNOTEL_EDGE    = 'k';              % marker edge color
    SNOTEL_FONT    = 7;                % label font size
    TITLE_FONT     = 14;
    LABEL_FONT     = 11;
    FIG_SIZE       = [100 100 1000 800];
    % ================================================================

    %% Parse inputs
    ip = inputParser;
    % Mode 1 inputs
    addParameter(ip, 'nc_file',      '',  @ischar);
    addParameter(ip, 'day_index',    1,   @isnumeric);
    addParameter(ip, 'shapefile',    '',  @ischar);
    addParameter(ip, 'mask_file',    '',  @ischar);
    % Mode 2 inputs
    addParameter(ip, 'frame_data',   [],  @isnumeric);
    addParameter(ip, 'frame_date',   NaN, @isnumeric);
    addParameter(ip, 'plot_context', [],  @isstruct);
    % Common inputs
    addParameter(ip, 'fig_handle',   [],  @(x) isempty(x) || ishgfigure(x));
    addParameter(ip, 'visible',      '',  @ischar);
    addParameter(ip, 'save_png',     '',  @ischar);
    addParameter(ip, 'water_year',   0,   @isnumeric);
    parse(ip, varargin{:});
    args = ip.Results;

    %% Constants
    R_EARTH    = 6371007.181;
    TILE_SIZE  = 1111950.5197;
    NPIX       = 2400;
    PIXEL_SIZE = TILE_SIZE / NPIX;
    H_TILE     = 9;
    V_TILE     = 4;
    UTM_ZONE   = 11;
    UTM_RES    = 500;   % output grid resolution (m)

    %% ================================================================
    %  SNOTEL stations (edit / add / remove as needed)
    %  ================================================================
    snotel_names   = {'Bogus Basin','Mores Creek Summit','Graham Guard Sta.', ...
                      'Jackson Peak','Atlanta Summit','Trinity Mtn.', ...
                      'Banner Summit','Prairie'};
    snotel_lat     = [43+46/60, 43+56/60, 43+57/60, 44+3/60, ...
                      43+45/60, 43+38/60, 44+18/60, 43+49/60];
    snotel_lon     = [-(116+6/60), -(115+40/60), -(115+16/60), -(115+27/60), ...
                      -(115+14/60), -(115+26/60), -(115+14/60), -(115+54/60)];
    snotel_id      = [978, 637, 496, 550, 306, 830, 312, 700]; %#ok<NASGU>
    snotel_elev_ft = [6370, 6090, 5680, 7060, 7570, 7790, 7040, 5580]; %#ok<NASGU>

    %% ================================================================
    %  Determine mode and prepare data
    %  ================================================================
    if ~isempty(args.plot_context)
        % ----- MODE 2: pre-processed data from the animation script -----
        ctx       = args.plot_context;
        frame_utm = args.frame_data;
        dt        = args.frame_date;
        wy        = args.water_year;
        ek        = ctx.ek;
        nk        = ctx.nk;
        poly_utm_e = ctx.poly_utm_e;
        poly_utm_n = ctx.poly_utm_n;
        snotel_e   = ctx.snotel_e;
        snotel_n   = ctx.snotel_n;
        in_view    = ctx.in_view;
        utm_zone_  = ctx.utm_zone;
        if isempty(args.visible); vis = 'off'; else; vis = args.visible; end

    else
        % ----- MODE 1: read directly from a netCDF file -----
        if isempty(args.nc_file)
            error('plot_fsca_frame:noInput', ...
                  'Provide either ''nc_file'' (Mode 1) or ''frame_data''+''plot_context'' (Mode 2).');
        end
        if isempty(args.visible); vis = 'on'; else; vis = args.visible; end
        utm_zone_ = UTM_ZONE;

        % --- Load or build basin mask ---
        if ~isempty(args.mask_file) && exist(args.mask_file, 'file')
            fprintf('Loading cached mask: %s\n', args.mask_file);
            m = load(args.mask_file);
        elseif ~isempty(args.shapefile)
            fprintf('Building basin mask from shapefile...\n');
            m = build_mask(args.shapefile, H_TILE, V_TILE, ...
                           R_EARTH, TILE_SIZE, PIXEL_SIZE, NPIX, UTM_ZONE);
            % Cache it next to the shapefile
            [sp, sn] = fileparts(args.shapefile);
            if isempty(sp); sp = '.'; end
            cache = fullfile(sp, sprintf('basin_mask_%s_h%02dv%02d.mat', ...
                    sn, H_TILE, V_TILE));
            save(cache, '-struct', 'm');
            fprintf('Mask cached: %s\n', cache);
        else
            error('plot_fsca_frame:noMask', ...
                  'Provide ''shapefile'' or ''mask_file'' for Mode 1.');
        end

        mask     = m.mask;
        row_min  = m.row_min;  row_max  = m.row_max;
        col_min  = m.col_min;  col_max  = m.col_max;
        nrows    = m.nrows;    ncols    = m.ncols;
        x_origin = m.x_origin; y_origin = m.y_origin;
        all_poly_x_sin = m.all_poly_x_sin;
        all_poly_y_sin = m.all_poly_y_sin;

        % --- Read fSCA variable from netCDF ---
        fprintf('Reading: %s\n', args.nc_file);
        info = ncinfo(args.nc_file);
        var_names = {info.Variables.Name};

        % Search for fSCA variable — prefer exact 'snow_fraction' first
        fsca_var = '';

        % Priority 1: exact match for 'snow_fraction' (on-the-ground, gap-filled)
        if any(strcmp(var_names, 'snow_fraction'))
            fsca_var = 'snow_fraction';
        end

        % Priority 2: other exact name matches
        if isempty(fsca_var)
            other_exact = {'fSCA','fsca','FSCA','snow_fraction_on_ground', ...
                           'Snow_Fraction','fractional_snow_cover', ...
                           'snow_cover_fraction'};
            for ii = 1:length(other_exact)
                if any(strcmp(var_names, other_exact{ii}))
                    fsca_var = other_exact{ii}; break;
                end
            end
        end

        % Priority 3: fall back to viewable_snow_fraction
        if isempty(fsca_var) && any(strcmp(var_names, 'viewable_snow_fraction'))
            fsca_var = 'viewable_snow_fraction';
            fprintf('  NOTE: Using viewable_snow_fraction (on-ground version not found)\n');
        end

        % Priority 4: any variable with snow+frac in the name
        if isempty(fsca_var)
            for ii = 1:length(var_names)
                if contains(lower(var_names{ii}),'snow') && contains(lower(var_names{ii}),'frac')
                    fsca_var = var_names{ii}; break;
                end
            end
        end

        if isempty(fsca_var)
            error('No fSCA variable found. Available: %s', strjoin(var_names,', '));
        end
        fprintf('  Using variable: %s\n', fsca_var);

        vi = [];
        for ii = 1:length(info.Variables)
            if strcmp(info.Variables(ii).Name, fsca_var); vi = info.Variables(ii); break; end
        end
        fprintf('  Variable: %s, size: [%s]\n', fsca_var, num2str([vi.Dimensions.Length]));

        tidx = args.day_index;
        dt = NaN;

        if length(vi.Size) == 3
            % --- Multi-day file: extract one time slice ---
            dim_names = {vi.Dimensions.Name};
            time_dim = 0;
            for d = 1:3
                if contains(lower(dim_names{d}),'time') || contains(lower(dim_names{d}),'day')
                    time_dim = d; break;
                end
            end
            if time_dim == 0; time_dim = 3; end
            spatial_dims = setdiff(1:3, time_dim);

            s = ones(1,3); c = [vi.Dimensions.Length];
            s(time_dim) = tidx; c(time_dim) = 1;
            s(spatial_dims(1)) = row_min; c(spatial_dims(1)) = nrows;
            s(spatial_dims(2)) = col_min; c(spatial_dims(2)) = ncols;
            slice = squeeze(ncread(args.nc_file, fsca_var, s, c));
            if size(slice,1) ~= nrows; slice = slice'; end

            % Try to read date from time variable
            for tname = {'time','Time','date','day'}
                if any(strcmp(var_names, tname{1}))
                    try
                        tv = ncread(args.nc_file, tname{1}, tidx, 1);
                        ti = ncinfo(args.nc_file, tname{1}); tu = '';
                        for a = 1:length(ti.Attributes)
                            if strcmp(ti.Attributes(a).Name,'units'); tu = ti.Attributes(a).Value; break; end
                        end
                        if contains(tu,'days since')
                            pts = strsplit(tu); dt = datenum(pts{3}) + double(tv);
                        elseif contains(tu,'seconds since')
                            pts = strsplit(tu); dt = datenum(pts{3}) + double(tv)/86400;
                        end
                    catch; end
                    break;
                end
            end

        elseif length(vi.Size) == 2
            % --- Single-day file (e.g. SPIRES_HIST_h09v04_MOD09GA061_20200315_V1.0.nc) ---
            slice = ncread(args.nc_file, fsca_var, [row_min col_min], [nrows ncols]);
            if size(slice,1) ~= nrows; slice = slice'; end
        end

        % --- Read a companion variable to detect unobserved pixels ---
        %     In SPIRES daily files, unobserved pixels (between swath lines)
        %     are stored as 0 in all variables, NOT as _FillValue.
        %     We detect these by checking: if snow_fraction==0 AND
        %     grain_size==0, the pixel was not observed (a truly snow-free
        %     observed pixel would have grain_size == _FillValue/255).
        unobs_mask = false(size(slice));
        try
            if any(strcmp(var_names, 'grain_size'))
                gs = ncread(args.nc_file, 'grain_size', [row_min col_min 1], [nrows ncols 1]);
                gs = squeeze(gs);
                if size(gs,1) ~= nrows; gs = gs'; end
                % Unobserved: snow_fraction==0 AND grain_size==0
                unobs_mask = (slice == 0) & (gs == 0);
            elseif any(strcmp(var_names, 'viewable_snow_fraction'))
                vsf = ncread(args.nc_file, 'viewable_snow_fraction', [row_min col_min 1], [nrows ncols 1]);
                vsf = squeeze(vsf);
                if size(vsf,1) ~= nrows; vsf = vsf'; end
                % Unobserved: both snow_fraction and viewable_snow_fraction are 0
                unobs_mask = (slice == 0) & (vsf == 0);
            end
        catch
            % If companion read fails, fall back to no unobs detection
        end

        % --- Extract date from filename if not found in file ---
        if isnan(dt)
            [~, fn] = fileparts(args.nc_file);
            % Try YYYYMMDD pattern (matches filenames like ..._20200315_...)
            tok = regexp(fn, '(20\d{6})', 'tokens');
            if ~isempty(tok)
                dt = datenum(tok{1}{1}, 'yyyymmdd');
            else
                % Try YYYYDDD pattern
                tok = regexp(fn, '(20\d{4})', 'tokens');
                if ~isempty(tok)
                    yr_t = str2double(tok{1}{1}(1:4));
                    doy_t = str2double(tok{1}{1}(5:7));
                    if doy_t >= 1 && doy_t <= 366
                        dt = datenum(yr_t, 1, doy_t);
                    end
                end
            end
        end

        % --- Read variable attributes for proper no-data handling ---
        fill_val = NaN;
        valid_min = 0;
        valid_max = 1;  % default: fSCA should be 0–1
        for a = 1:length(vi.Attributes)
            aname = vi.Attributes(a).Name;
            aval  = vi.Attributes(a).Value;
            if strcmp(aname, '_FillValue');   fill_val  = double(aval); end
            if strcmp(aname, 'missing_value'); fill_val = double(aval); end
            if strcmp(aname, 'valid_min');    valid_min = double(aval); end
            if strcmp(aname, 'valid_max');    valid_max = double(aval); end
            if strcmp(aname, 'valid_range')
                valid_min = double(aval(1));
                valid_max = double(aval(2));
            end
        end

        % If valid_max > 1, data is likely stored as percentage or scaled int
        if valid_max > 1.5
            is_pct = true;
        else
            is_pct = false;
        end

        % --- Diagnostic: show value distribution BEFORE masking ---
        raw_in_basin = slice(mask);
        raw_finite   = raw_in_basin(isfinite(raw_in_basin));
        n_unobs_basin = sum(unobs_mask(mask));
        fprintf('  --- Data diagnostics (within basin, before filtering) ---\n');
        fprintf('  _FillValue = %g,  valid_range = [%g, %g]\n', fill_val, valid_min, valid_max);
        fprintf('  Total pixels in basin: %d\n', numel(raw_in_basin));
        fprintf('  NaN pixels:  %d\n', sum(isnan(raw_in_basin)));
        fprintf('  Unobserved pixels (swath gaps): %d (%.1f%%)\n', ...
                n_unobs_basin, 100*n_unobs_basin/numel(raw_in_basin));
        if ~isempty(raw_finite)
            fprintf('  Finite min:  %g\n', min(raw_finite));
            fprintf('  Finite max:  %g\n', max(raw_finite));
            fprintf('  Finite mean (all): %g\n', mean(raw_finite));
            % Mean excluding unobserved
            obs_vals = raw_in_basin(~unobs_mask(mask) & isfinite(raw_in_basin));
            if ~isempty(obs_vals)
                fprintf('  Finite mean (observed only): %g\n', mean(obs_vals));
            end
        end
        fprintf('  ---------------------------------------------------\n');

        % --- Apply mask and filter to valid range ---
        slice(~mask) = NaN;

        % Mark unobserved pixels as NaN (the key fix!)
        slice(unobs_mask) = NaN;

        % Replace fill values
        if isfinite(fill_val)
            slice(slice == fill_val) = NaN;
        end

        % Replace anything outside the valid range with NaN
        slice(slice < valid_min | slice > valid_max) = NaN;

        % Convert percentage to fraction if needed
        if is_pct
            slice = slice / 100.0;
        end

        % --- Reproject sinusoidal → UTM ---
        %     The sinusoidal projection is nonlinear: at higher latitudes,
        %     the same UTM easting maps to a different sinusoidal column.
        %     We must index into the FULL TILE (cols 1–2400) and read an
        %     expanded data range to avoid clipping the NW/NE corners.
        fprintf('  Reprojecting to UTM Zone %dN...\n', UTM_ZONE);
        sin_x_vec = x_origin + ((col_min:col_max) - 0.5) * PIXEL_SIZE;
        sin_y_vec = y_origin - ((row_min:row_max) - 0.5) * PIXEL_SIZE;

        % Get UTM extent from basin mask pixels
        [sin_x_g, sin_y_g] = meshgrid(sin_x_vec, sin_y_vec);
        lat_g = rad2deg(sin_y_g ./ R_EARTH);
        cos_l = cos(sin_y_g ./ R_EARTH); cos_l(abs(cos_l)<1e-10) = 1e-10;
        lon_g = rad2deg(sin_x_g ./ (R_EARTH .* cos_l));
        [utm_eg, utm_ng] = latlon2utm(lat_g, lon_g, UTM_ZONE);

        utm_e_valid = utm_eg(mask); utm_n_valid = utm_ng(mask);
        ue_min = min(utm_e_valid) - PIXEL_SIZE;
        ue_max = max(utm_e_valid) + PIXEL_SIZE;
        un_min = min(utm_n_valid) - PIXEL_SIZE;
        un_max = max(utm_n_valid) + PIXEL_SIZE;

        % Build regular UTM output grid
        utm_e_out = ue_min : UTM_RES : ue_max;
        utm_n_out = un_max : -UTM_RES : un_min;
        [UE, UN] = meshgrid(utm_e_out, utm_n_out);
        nr2 = length(utm_n_out); nc2 = length(utm_e_out);

        % Inverse-project UTM grid to FULL TILE sinusoidal indices
        [la2, lo2] = utm2latlon(UE(:), UN(:), UTM_ZONE);
        [sx2, sy2] = latlon2sin(la2, lo2, R_EARTH);

        tile_col = round((sx2 - x_origin) / PIXEL_SIZE + 0.5);
        tile_row = round((y_origin - sy2) / PIXEL_SIZE + 0.5);

        % Determine expanded read range
        valid_tile = tile_col >= 1 & tile_col <= NPIX & ...
                     tile_row >= 1 & tile_row <= NPIX;
        rd_col_min = max(1,    min(tile_col(valid_tile)));
        rd_col_max = min(NPIX, max(tile_col(valid_tile)));
        rd_row_min = max(1,    min(tile_row(valid_tile)));
        rd_row_max = min(NPIX, max(tile_row(valid_tile)));
        rd_nrows = rd_row_max - rd_row_min + 1;
        rd_ncols = rd_col_max - rd_col_min + 1;

        fprintf('  Original mask subset: cols [%d:%d], rows [%d:%d]\n', ...
                col_min, col_max, row_min, row_max);
        fprintf('  Expanded read range:  cols [%d:%d], rows [%d:%d]\n', ...
                rd_col_min, rd_col_max, rd_row_min, rd_row_max);

        % Convert to local (expanded-read) indices
        local_col = tile_col - rd_col_min + 1;
        local_row = tile_row - rd_row_min + 1;
        vl = valid_tile & local_col >= 1 & local_col <= rd_ncols & ...
                          local_row >= 1 & local_row <= rd_nrows;
        local_col = reshape(local_col, nr2, nc2);
        local_row = reshape(local_row, nr2, nc2);
        vl = reshape(vl, nr2, nc2);

        % Build lookup — test basin membership using the sinusoidal
        % coordinates of each UTM output cell (not the pre-built mask,
        % which only covers the tight subset and misses expanded pixels)
        fprintf('  Building lookup with basin polygon test...\n');
        sli = zeros(nr2, nc2);
        sx2_grid = reshape(sx2, nr2, nc2);
        sy2_grid = reshape(sy2, nr2, nc2);
        for rr = 1:nr2
            for cc = 1:nc2
                if vl(rr,cc)
                    % Test if this sinusoidal point is inside the basin
                    for pp = 1:length(all_poly_x_sin)
                        if inpolygon(sx2_grid(rr,cc), sy2_grid(rr,cc), ...
                                     all_poly_x_sin{pp}, all_poly_y_sin{pp})
                            sli(rr,cc) = sub2ind([rd_nrows rd_ncols], ...
                                                  local_row(rr,cc), local_col(rr,cc));
                            break;
                        end
                    end
                end
            end
        end

        % Re-read the data with expanded range
        fprintf('  Re-reading data with expanded range [%dx%d]...\n', rd_nrows, rd_ncols);
        info_v = ncinfo(args.nc_file, fsca_var);
        if length(info_v.Size) == 3
            dim_names_v = {info_v.Dimensions.Name};
            td = 0;
            for d = 1:3
                if contains(lower(dim_names_v{d}),'time') || contains(lower(dim_names_v{d}),'day')
                    td = d; break;
                end
            end
            if td == 0; td = 3; end
            sd = setdiff(1:3, td);
            ss = ones(1,3); cc_rd = [info_v.Dimensions.Length];
            ss(td) = tidx; cc_rd(td) = 1;
            ss(sd(1)) = rd_row_min; cc_rd(sd(1)) = rd_nrows;
            ss(sd(2)) = rd_col_min; cc_rd(sd(2)) = rd_ncols;
            slice_exp = double(squeeze(ncread(args.nc_file, fsca_var, ss, cc_rd)));
        else
            slice_exp = double(ncread(args.nc_file, fsca_var, ...
                              [rd_row_min rd_col_min], [rd_nrows rd_ncols]));
        end
        if size(slice_exp,1) ~= rd_nrows; slice_exp = slice_exp'; end

        % Apply filtering on expanded data (basin clipping done by lookup)
        if isfinite(fill_val); slice_exp(slice_exp == fill_val) = NaN; end
        slice_exp(slice_exp < valid_min | slice_exp > valid_max) = NaN;

        % Detect unobserved pixels on expanded data
        try
            if any(strcmp(var_names, 'grain_size'))
                if length(info_v.Size) == 3
                    gs_exp = double(squeeze(ncread(args.nc_file, 'grain_size', ss, cc_rd)));
                else
                    gs_exp = double(ncread(args.nc_file, 'grain_size', ...
                                   [rd_row_min rd_col_min], [rd_nrows rd_ncols]));
                end
                if size(gs_exp,1) ~= rd_nrows; gs_exp = gs_exp'; end
                slice_exp((slice_exp == 0) & (gs_exp == 0)) = NaN;
            end
        catch; end

        if is_pct; slice_exp = slice_exp / 100.0; end

        % Reproject
        frame_utm = NaN(nr2, nc2);
        frame_utm(sli > 0) = slice_exp(sli(sli > 0));

        ek = [ue_min ue_max] / 1000;
        nk = [un_min un_max] / 1000;

        % Basin polygons in UTM
        poly_utm_e = {}; poly_utm_n = {};
        for pp = 1:length(all_poly_x_sin)
            [plat, plon] = sin2latlon(all_poly_x_sin{pp}, all_poly_y_sin{pp}, R_EARTH);
            [pe, pn] = latlon2utm(plat, plon, UTM_ZONE);
            poly_utm_e{end+1} = pe; poly_utm_n{end+1} = pn; %#ok<AGROW>
        end

        % SNOTEL in UTM
        [snotel_e, snotel_n] = latlon2utm(snotel_lat, snotel_lon, UTM_ZONE);
        buf = 2000;
        in_view = snotel_e >= (ue_min - buf) & snotel_e <= (ue_max + buf) & ...
                  snotel_n >= (un_min - buf) & snotel_n <= (un_max + buf);

        % Determine water year from date
        if ~isnan(dt)
            [yy, mm] = datevec(dt);
            if mm >= 10; wy = yy + 1; else; wy = yy; end
        else
            wy = 0;
        end
    end

    if args.water_year > 0; wy = args.water_year; end

    %% ================================================================
    %  PLOT THE FRAME — edit this section to change the look
    %  ================================================================
    if isempty(args.fig_handle)
        fig = figure('Position', FIG_SIZE, 'Visible', vis);
    else
        fig = args.fig_handle;
        clf(fig);
    end
    ax = axes('Parent', fig);

    % ---- fSCA image ----
    imagesc(ax, ek, [nk(2) nk(1)], frame_utm);
    set(ax, 'YDir', 'normal');
    caxis(ax, FSCA_RANGE);
    colormap(ax, COLORMAP);
    hold(ax, 'on');

    % ---- Basin boundary ----
    for pp = 1:length(poly_utm_e)
        plot(ax, poly_utm_e{pp}/1000, poly_utm_n{pp}/1000, ...
             BASIN_COLOR, 'LineWidth', BASIN_WIDTH);
    end

    % ---- SNOTEL stations ----
    for si = 1:length(snotel_names)
        if in_view(si)
            plot(ax, snotel_e(si)/1000, snotel_n(si)/1000, ...
                 SNOTEL_MARKER, 'MarkerSize', SNOTEL_SIZE, ...
                 'MarkerEdgeColor', SNOTEL_EDGE, ...
                 'MarkerFaceColor', SNOTEL_FACE, 'LineWidth', 1.0);
            text(ax, snotel_e(si)/1000 + 0.8, snotel_n(si)/1000, ...
                 snotel_names{si}, ...
                 'FontSize', SNOTEL_FONT, 'FontWeight', 'bold', ...
                 'Color', [0.1 0.1 0.1], ...
                 'BackgroundColor', [1 1 1 0.7], ...
                 'EdgeColor', 'none', 'Margin', 0.5);
        end
    end

    % ---- Axes ----
    axis(ax, 'equal');
    xlim(ax, ek);
    ylim(ax, nk);
    xlabel(ax, sprintf('Easting (km, UTM Zone %dN)', utm_zone_), ...
           'FontSize', LABEL_FONT);
    ylabel(ax, 'Northing (km)', 'FontSize', LABEL_FONT);

    % ---- Title ----
    if ~isnan(dt)
        date_str = datestr(dt, 'mmmm dd, yyyy');
    else
        date_str = sprintf('Day index %d', args.day_index);
    end
    if wy > 0
        title(ax, sprintf('Boise River Basin - Fractional Snow Covered Area  (WY %d)\n%s', ...
              wy, date_str), 'FontSize', TITLE_FONT, 'FontWeight', 'bold');
    else
        title(ax, sprintf('Boise River Basin - Fractional Snow Covered Area\n%s', ...
              date_str), 'FontSize', TITLE_FONT, 'FontWeight', 'bold');
    end

    % ---- Colorbar ----
    cb = colorbar(ax);
    cb.Label.String = 'Fractional Snow Covered Area';
    cb.Label.FontSize = LABEL_FONT;

    % ---- Basin-mean annotation ----
    vd = frame_utm(isfinite(frame_utm));
    if ~isempty(vd)
        text(ax, ek(1) + 0.02*diff(ek), nk(1) + 0.05*diff(nk), ...
             sprintf('Basin mean fSCA: %.3f', mean(vd)), ...
             'FontSize', LABEL_FONT, ...
             'BackgroundColor', [1 1 1 0.8], 'EdgeColor', [0.5 0.5 0.5]);
    end

    hold(ax, 'off');
    drawnow;

    % ---- Optionally save to PNG ----
    if ~isempty(args.save_png)
        exportgraphics(fig, args.save_png, 'Resolution', 150);
        fprintf('Saved: %s\n', args.save_png);
    end
end


%% ====================================================================
%  LOCAL HELPERS (only used in Mode 1, self-contained)
%  ====================================================================

function m = build_mask(shapefile_path, h_tile, v_tile, R_EARTH, TS, PS, NP, utm_zone)
    S = shaperead(shapefile_path);
    bbox = [S(1).BoundingBox];
    is_utm = bbox(1,1) > 1000;

    all_poly_x_sin = {}; all_poly_y_sin = {};
    all_x = []; all_y = [];
    for k = 1:length(S)
        px = S(k).X; py = S(k).Y;
        ni = [0, find(isnan(px)), length(px)+1];
        for pp = 1:(length(ni)-1)
            i1 = ni(pp)+1; i2 = ni(pp+1)-1;
            if i2 < i1; continue; end
            ppx = px(i1:i2); ppy = py(i1:i2);
            if is_utm
                [la, lo] = utm2latlon(ppx, ppy, utm_zone);
            else
                la = ppy; lo = ppx;
            end
            [sx, sy] = latlon2sin(la, lo, R_EARTH);
            all_poly_x_sin{end+1} = sx; %#ok<AGROW>
            all_poly_y_sin{end+1} = sy; %#ok<AGROW>
            all_x = [all_x; sx(:)]; %#ok<AGROW>
            all_y = [all_y; sy(:)]; %#ok<AGROW>
        end
    end

    x_origin = -R_EARTH * pi + h_tile * TS;
    y_origin =  R_EARTH * pi / 2 - v_tile * TS;

    col_min = max(1, floor((min(all_x) - x_origin) / PS) - 1);
    col_max = min(NP, ceil((max(all_x) - x_origin) / PS) + 2);
    row_min = max(1, floor((y_origin - max(all_y)) / PS) - 1);
    row_max = min(NP, ceil((y_origin - min(all_y)) / PS) + 2);
    nrows = row_max - row_min + 1;
    ncols = col_max - col_min + 1;

    fprintf('  Mask dimensions: %d x %d\n', nrows, ncols);
    mask = false(nrows, ncols);
    for r = 1:nrows
        for c = 1:ncols
            mpx = x_origin + (col_min + c - 1 - 0.5) * PS;
            mpy = y_origin - (row_min + r - 1 - 0.5) * PS;
            for pp = 1:length(all_poly_x_sin)
                if inpolygon(mpx, mpy, all_poly_x_sin{pp}, all_poly_y_sin{pp})
                    mask(r, c) = true;
                    break;
                end
            end
        end
    end
    fprintf('  Basin contains %d pixels\n', sum(mask(:)));

    m = struct('mask', mask, 'row_min', row_min, 'row_max', row_max, ...
               'col_min', col_min, 'col_max', col_max, ...
               'nrows', nrows, 'ncols', ncols, ...
               'all_poly_x_sin', {all_poly_x_sin}, ...
               'all_poly_y_sin', {all_poly_y_sin}, ...
               'utm_zone', utm_zone, ...
               'x_origin', x_origin, 'y_origin', y_origin);
end


function [lat, lon] = utm2latlon(easting, northing, zone)
    a = 6378137; f = 1/298.257223563;
    e = sqrt(2*f - f^2); ep = e/sqrt(1 - e^2); k0 = 0.9996;
    x = easting - 500000; y = northing;
    M = y / k0;
    mu = M / (a * (1 - e^2/4 - 3*e^4/64 - 5*e^6/256));
    e1 = (1 - sqrt(1 - e^2)) / (1 + sqrt(1 - e^2));
    phi1 = mu + (3*e1/2 - 27*e1^3/32).*sin(2*mu) ...
              + (21*e1^2/16 - 55*e1^4/32).*sin(4*mu) ...
              + (151*e1^3/96).*sin(6*mu);
    N1 = a./sqrt(1 - e^2*sin(phi1).^2);
    T1 = tan(phi1).^2;
    C1 = ep^2*cos(phi1).^2;
    R1 = a*(1 - e^2)./(1 - e^2*sin(phi1).^2).^1.5;
    D  = x./(N1*k0);
    lat_rad = phi1 - (N1.*tan(phi1)./R1) .* ( ...
        D.^2/2 - (5 + 3*T1 + 10*C1 - 4*C1.^2 - 9*ep^2).*D.^4/24 ...
        + (61 + 90*T1 + 298*C1 + 45*T1.^2 - 252*ep^2 - 3*C1.^2).*D.^6/720);
    lon_rad = (D - (1 + 2*T1 + C1).*D.^3/6 ...
        + (5 - 2*C1 + 28*T1 - 3*C1.^2 + 8*ep^2 + 24*T1.^2).*D.^5/120) ...
        ./ cos(phi1);
    lon0 = deg2rad((zone - 1)*6 - 180 + 3);
    lat = rad2deg(lat_rad);
    lon = rad2deg(lon_rad + lon0);
end


function [easting, northing] = latlon2utm(lat, lon, zone)
    a = 6378137; f = 1/298.257223563;
    e = sqrt(2*f - f^2); ep = e/sqrt(1 - e^2); k0 = 0.9996;
    lon0 = (zone - 1)*6 - 180 + 3;
    lr = deg2rad(lat); lnr = deg2rad(lon); l0 = deg2rad(lon0);
    N = a./sqrt(1 - e^2*sin(lr).^2);
    T = tan(lr).^2; C = ep^2*cos(lr).^2;
    A = (lnr - l0).*cos(lr);
    M = a*((1 - e^2/4 - 3*e^4/64 - 5*e^6/256).*lr ...
          - (3*e^2/8 + 3*e^4/32 + 45*e^6/1024).*sin(2*lr) ...
          + (15*e^4/256 + 45*e^6/1024).*sin(4*lr) ...
          - (35*e^6/3072).*sin(6*lr));
    easting  = k0*N.*(A + (1 - T + C).*A.^3/6 ...
               + (5 - 18*T + T.^2 + 72*C - 58*ep^2).*A.^5/120) + 500000;
    northing = k0*(M + N.*tan(lr).*(A.^2/2 + (5 - T + 9*C + 4*C.^2).*A.^4/24 ...
               + (61 - 58*T + T.^2 + 600*C - 330*ep^2).*A.^6/720));
end


function [x, y] = latlon2sin(lat, lon, R)
    lr = deg2rad(lat); lnr = deg2rad(lon);
    x = R .* lnr .* cos(lr);
    y = R .* lr;
end


function [lat, lon] = sin2latlon(x, y, R)
    lr = y ./ R;
    lat = rad2deg(lr);
    cl = cos(lr); cl(abs(cl) < 1e-10) = 1e-10;
    lon = rad2deg(x ./ (R .* cl));
end


function tf = ishgfigure(h)
    tf = isgraphics(h, 'figure');
end
