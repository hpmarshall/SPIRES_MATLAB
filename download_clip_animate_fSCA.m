%% download_clip_animate_fSCA.m
%
% Downloads SPIRES HIST V01 fractional snow covered area (fSCA) data from:
%     ftp://dtn.rc.colorado.edu/shares/snow-today/gridded_data/SPIRES_HIST_V01
%
% Clips the data to the Boise River Basin (BRB) using a user-supplied
% shapefile, and creates an MP4 animation of daily fSCA for a user-specified
% year.
%
% Requirements:
%     - MATLAB R2018b+ (for ftp, ncread, VideoWriter, shaperead)
%     - Mapping Toolbox (for shaperead; or use readgeotable in R2022+)
%
% Usage:
%     >> download_clip_animate_fSCA(2020, 'BRB_outline.shp')
%     >> download_clip_animate_fSCA(2015, 'BRB_outline.shp', 'fps', 15)
%
% Data Citation:
%     Rittger, K., Lenard, S. J., Palomaki, R. T., Bair, E. H., Dozier, J. &
%     Mankoff, K. (2025). Historical MODIS/Terra L3 Global Daily 500m SIN Grid
%     Snow Cover, Snow Albedo, and Snow Surface Properties. (SPIRES_HIST, V1).
%     [Data Set]. Boulder, CO. NSIDC. https://doi.org/10.7265/a3vr-c014

function download_clip_animate_fSCA(year, shapefile_path, varargin)

    p = inputParser;
    addRequired(p, 'year', @(x) isnumeric(x) && x >= 2000 && x <= 2025);
    addRequired(p, 'shapefile_path', @ischar);
    addParameter(p, 'output_dir', './spires_output', @ischar);
    addParameter(p, 'fps', 10, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'tile', 'h09v04', @ischar);
    parse(p, year, shapefile_path, varargin{:});

    output_dir = p.Results.output_dir;
    fps        = p.Results.fps;
    tile       = p.Results.tile;

    FTP_HOST   = 'dtn.rc.colorado.edu';
    FTP_DIR    = '/shares/snow-today/gridded_data/SPIRES_HIST_V01';
    R_EARTH    = 6371007.181;
    TILE_SIZE  = 1111950.5197;
    NPIX       = 2400;
    PIXEL_SIZE = TILE_SIZE / NPIX;
    h_tile = str2double(tile(2:3));
    v_tile = str2double(tile(5:6));

    fprintf('Year: %d, Tile: %s, Shapefile: %s\n', year, tile, shapefile_path);

    %% Step 1: Read shapefile and create basin mask
    S = shaperead(shapefile_path);
    bbox = [S(1).BoundingBox];
    utm_zone = 11;
    is_utm = bbox(1,1) > 1000;

    all_poly_x_sin = {}; all_poly_y_sin = {};
    all_x_sin = []; all_y_sin = [];
    for k = 1:length(S)
        px = S(k).X; py = S(k).Y;
        nan_idx = [0, find(isnan(px)), length(px)+1];
        for pp = 1:(length(nan_idx)-1)
            i1 = nan_idx(pp)+1; i2 = nan_idx(pp+1)-1;
            if i2 < i1; continue; end
            if is_utm
                [lat, lon] = utm2latlon(px(i1:i2), py(i1:i2), utm_zone);
            else
                lat = py(i1:i2); lon = px(i1:i2);
            end
            [sx, sy] = latlon2sin(lat, lon, R_EARTH);
            all_poly_x_sin{end+1} = sx; all_poly_y_sin{end+1} = sy;
            all_x_sin = [all_x_sin; sx(:)]; all_y_sin = [all_y_sin; sy(:)];
        end
    end

    x_origin = -R_EARTH * pi + h_tile * TILE_SIZE;
    y_origin =  R_EARTH * pi / 2 - v_tile * TILE_SIZE;
    col_min = max(1, floor((min(all_x_sin)-x_origin)/PIXEL_SIZE)-1);
    col_max = min(NPIX, ceil((max(all_x_sin)-x_origin)/PIXEL_SIZE)+2);
    row_min = max(1, floor((y_origin-max(all_y_sin))/PIXEL_SIZE)-1);
    row_max = min(NPIX, ceil((y_origin-min(all_y_sin))/PIXEL_SIZE)+2);
    nrows = row_max - row_min + 1; ncols = col_max - col_min + 1;

    mask = false(nrows, ncols);
    for r = 1:nrows
        for c = 1:ncols
            mpx = x_origin + (col_min+c-1-0.5)*PIXEL_SIZE;
            mpy = y_origin - (row_min+r-1-0.5)*PIXEL_SIZE;
            for pp = 1:length(all_poly_x_sin)
                if inpolygon(mpx, mpy, all_poly_x_sin{pp}, all_poly_y_sin{pp})
                    mask(r,c) = true; break;
                end
            end
        end
    end
    fprintf('Basin: %d pixels\n', sum(mask(:)));

    %% Step 2-3: Download and extract (abbreviated — see full version)
    % ... FTP download and ncread logic ...

    %% Step 4: Plot in lat/lon
    x_left  = x_origin + (col_min-1)*PIXEL_SIZE;
    x_right = x_origin + col_max*PIXEL_SIZE;
    y_top   = y_origin - (row_min-1)*PIXEL_SIZE;
    y_bot   = y_origin - row_max*PIXEL_SIZE;
    [lat_top, lon_left] = sin2latlon(x_left, y_top, R_EARTH);
    [lat_bot, lon_right] = sin2latlon(x_right, y_bot, R_EARTH);

    % ... create_movie with lat/lon axes ...
end

function [lat, lon] = utm2latlon(easting, northing, zone)
    a=6378137; f=1/298.257223563; e=sqrt(2*f-f^2); ep=e/sqrt(1-e^2); k0=0.9996;
    x=easting-500000; y=northing; M=y/k0;
    mu=M/(a*(1-e^2/4-3*e^4/64-5*e^6/256));
    e1=(1-sqrt(1-e^2))/(1+sqrt(1-e^2));
    p=mu+(3*e1/2-27*e1^3/32).*sin(2*mu)+(21*e1^2/16-55*e1^4/32).*sin(4*mu)+(151*e1^3/96).*sin(6*mu);
    N1=a./sqrt(1-e^2*sin(p).^2); T1=tan(p).^2; C1=ep^2*cos(p).^2;
    R1=a*(1-e^2)./(1-e^2*sin(p).^2).^1.5; D=x./(N1*k0);
    lat=rad2deg(p-(N1.*tan(p)./R1).*(D.^2/2-(5+3*T1+10*C1-4*C1.^2-9*ep^2).*D.^4/24 ...
        +(61+90*T1+298*C1+45*T1.^2-252*ep^2-3*C1.^2).*D.^6/720));
    lon=rad2deg((D-(1+2*T1+C1).*D.^3/6+(5-2*C1+28*T1-3*C1.^2+8*ep^2+24*T1.^2).*D.^5/120) ...
        ./cos(p)+deg2rad((zone-1)*6-180+3));
end
function [x,y] = latlon2sin(lat,lon,R)
    lr=deg2rad(lat); lnr=deg2rad(lon); x=R.*lnr.*cos(lr); y=R.*lr;
end
function [lat,lon] = sin2latlon(x,y,R)
    lr=y./R; lat=rad2deg(lr);
    cl=cos(lr); cl(abs(cl)<1e-10)=1e-10; lon=rad2deg(x./(R.*cl));
end
