#!/usr/bin/env python3
"""
download_clip_animate_fSCA.py

Downloads SPIRES HIST V01 fractional snow covered area (fSCA) data from:
    ftp://dtn.rc.colorado.edu/shares/snow-today/gridded_data/SPIRES_HIST_V01

Clips the data to the Boise River Basin (BRB) using a user-supplied shapefile,
and creates an MP4 animation of daily fSCA for a user-specified year.

Requirements:
    pip install netCDF4 numpy matplotlib shapefile pyproj imageio[ffmpeg]

Usage:
    python download_clip_animate_fSCA.py --year 2020 --shapefile BRB_outline.shp

Data Citation:
    Rittger, K., Lenard, S. J., Palomaki, R. T., Bair, E. H., Dozier, J. &
    Mankoff, K. (2025). Historical MODIS/Terra L3 Global Daily 500m SIN Grid
    Snow Cover, Snow Albedo, and Snow Surface Properties. (SPIRES_HIST, Version 1).
    [Data Set]. Boulder, Colorado USA. National Snow and Ice Data Center.
    https://doi.org/10.7265/a3vr-c014
"""

import os
import sys
import glob
import argparse
import ftplib
import math
from datetime import datetime, timedelta

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

try:
    import netCDF4 as nc
except ImportError:
    sys.exit("ERROR: netCDF4 not installed. Run: pip install netCDF4")

try:
    import shapefile
except ImportError:
    sys.exit("ERROR: pyshp not installed. Run: pip install pyshp")

try:
    import imageio.v2 as imageio
except ImportError:
    sys.exit("ERROR: imageio not installed. Run: pip install imageio[ffmpeg]")


# ============================================================================
# CONSTANTS
# ============================================================================
FTP_HOST = "dtn.rc.colorado.edu"
FTP_DIR = "/shares/snow-today/gridded_data/SPIRES_HIST_V01"

# MODIS sinusoidal projection parameters
R_EARTH = 6371007.181        # MODIS sphere radius (m)
TILE_SIZE = 1111950.5197     # tile edge length (m)
NPIX = 2400                  # pixels per tile edge
PIXEL_SIZE = TILE_SIZE / NPIX

# Tile h09v04 (contains the Boise River Basin)
H_TILE = 9
V_TILE = 4


# ============================================================================
# COORDINATE CONVERSION FUNCTIONS
# ============================================================================
def latlon_to_sinusoidal(lat, lon):
    """Convert geographic lat/lon (degrees) to MODIS sinusoidal x, y (meters)."""
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)
    x = R_EARTH * lon_rad * math.cos(lat_rad)
    y = R_EARTH * lat_rad
    return x, y


def sinusoidal_to_latlon(x, y):
    """Convert MODIS sinusoidal x, y (meters) to geographic lat/lon (degrees)."""
    lat_rad = y / R_EARTH
    lat = math.degrees(lat_rad)
    if abs(math.cos(lat_rad)) < 1e-10:
        lon = 0.0
    else:
        lon = math.degrees(x / (R_EARTH * math.cos(lat_rad)))
    return lat, lon


def tile_pixel_to_sinusoidal(h, v, row, col):
    """Convert tile h, v and pixel row, col to sinusoidal x, y (meters)."""
    x_origin = -R_EARTH * math.pi + h * TILE_SIZE
    y_origin = R_EARTH * math.pi / 2 - v * TILE_SIZE
    x = x_origin + (col + 0.5) * PIXEL_SIZE
    y = y_origin - (row + 0.5) * PIXEL_SIZE
    return x, y


# ============================================================================
# SHAPEFILE HANDLING
# ============================================================================
def read_shapefile_boundary(shp_path):
    """
    Read the shapefile and return the bounding box and polygon vertices
    in MODIS sinusoidal coordinates.

    The shapefile is assumed to be in UTM Zone 11N (EPSG:32611).
    If a .prj file is present, the script will attempt to detect the CRS.
    """
    sf = shapefile.Reader(shp_path)
    bbox = sf.bbox  # [xmin, ymin, xmax, ymax]

    # Detect CRS from .prj file if available
    prj_path = shp_path.replace('.shp', '.prj')
    is_utm = True  # default assumption for Idaho shapefiles
    utm_zone = 11

    if os.path.exists(prj_path):
        with open(prj_path, 'r') as f:
            prj_text = f.read()
        if 'GEOGCS' in prj_text and 'PROJCS' not in prj_text:
            is_utm = False
        elif 'UTM' in prj_text:
            # Try to extract zone number
            import re
            zone_match = re.search(r'Zone_(\d+)', prj_text)
            if zone_match:
                utm_zone = int(zone_match.group(1))

    # If coordinates are large (> 1000), assume projected (UTM)
    if bbox[0] > 1000:
        is_utm = True

    # Collect all polygon vertices
    all_polygons_sin = []
    for shape in sf.shapes():
        parts = list(shape.parts) + [len(shape.points)]
        for i in range(len(parts) - 1):
            ring = shape.points[parts[i]:parts[i+1]]
            ring_sin = []
            for pt in ring:
                if is_utm:
                    lat, lon = utm_to_latlon(pt[0], pt[1], zone=utm_zone)
                else:
                    lat, lon = pt[1], pt[0]
                x, y = latlon_to_sinusoidal(lat, lon)
                ring_sin.append((x, y))
            all_polygons_sin.append(ring_sin)

    # Compute bounding box in sinusoidal
    all_x = [p[0] for ring in all_polygons_sin for p in ring]
    all_y = [p[1] for ring in all_polygons_sin for p in ring]
    bbox_sin = [min(all_x), min(all_y), max(all_x), max(all_y)]

    return bbox_sin, all_polygons_sin


def utm_to_latlon(easting, northing, zone=11, northern=True):
    """Convert UTM coordinates to lat/lon (WGS84)."""
    a = 6378137.0
    f = 1 / 298.257223563
    e = math.sqrt(2 * f - f * f)
    e_prime = e / math.sqrt(1 - e * e)
    k0 = 0.9996
    x = easting - 500000.0
    y = northing
    if not northern:
        y -= 10000000.0

    M = y / k0
    mu = M / (a * (1 - e**2 / 4 - 3 * e**4 / 64 - 5 * e**6 / 256))
    e1 = (1 - math.sqrt(1 - e**2)) / (1 + math.sqrt(1 - e**2))

    phi1 = (mu + (3 * e1 / 2 - 27 * e1**3 / 32) * math.sin(2 * mu)
            + (21 * e1**2 / 16 - 55 * e1**4 / 32) * math.sin(4 * mu)
            + (151 * e1**3 / 96) * math.sin(6 * mu))

    N1 = a / math.sqrt(1 - e**2 * math.sin(phi1)**2)
    T1 = math.tan(phi1)**2
    C1 = e_prime**2 * math.cos(phi1)**2
    R1 = a * (1 - e**2) / (1 - e**2 * math.sin(phi1)**2)**1.5
    D = x / (N1 * k0)

    lat = phi1 - (N1 * math.tan(phi1) / R1) * (
        D**2 / 2
        - (5 + 3 * T1 + 10 * C1 - 4 * C1**2 - 9 * e_prime**2) * D**4 / 24
        + (61 + 90 * T1 + 298 * C1 + 45 * T1**2 - 252 * e_prime**2
           - 3 * C1**2) * D**6 / 720
    )

    lon = (D - (1 + 2 * T1 + C1) * D**3 / 6
           + (5 - 2 * C1 + 28 * T1 - 3 * C1**2 + 8 * e_prime**2
              + 24 * T1**2) * D**5 / 120) / math.cos(phi1)

    lon0 = math.radians((zone - 1) * 6 - 180 + 3)

    return math.degrees(lat), math.degrees(lon + lon0)


def point_in_polygon(x, y, polygon):
    """Ray-casting algorithm for point-in-polygon test."""
    n = len(polygon)
    inside = False
    j = n - 1
    for i in range(n):
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi):
            inside = not inside
        j = i
    return inside


def create_basin_mask(bbox_sin, polygons_sin, h, v):
    """
    Create a boolean mask for the basin within the MODIS tile.
    Returns the mask and the pixel row/col bounds for the subset.
    """
    x_origin = -R_EARTH * math.pi + h * TILE_SIZE
    y_origin = R_EARTH * math.pi / 2 - v * TILE_SIZE

    # Convert bbox to pixel indices (with buffer)
    col_min = max(0, int((bbox_sin[0] - x_origin) / PIXEL_SIZE) - 2)
    col_max = min(NPIX - 1, int((bbox_sin[2] - x_origin) / PIXEL_SIZE) + 2)
    row_min = max(0, int((y_origin - bbox_sin[3]) / PIXEL_SIZE) - 2)
    row_max = min(NPIX - 1, int((y_origin - bbox_sin[1]) / PIXEL_SIZE) + 2)

    nrows = row_max - row_min + 1
    ncols = col_max - col_min + 1
    mask = np.zeros((nrows, ncols), dtype=bool)

    print(f"  Creating basin mask: rows [{row_min}:{row_max}], cols [{col_min}:{col_max}]")
    print(f"  Mask dimensions: {nrows} x {ncols}")

    for r in range(nrows):
        for c in range(ncols):
            x, y = tile_pixel_to_sinusoidal(h, v, row_min + r, col_min + c)
            for poly in polygons_sin:
                if point_in_polygon(x, y, poly):
                    mask[r, c] = True
                    break

    n_pixels = np.sum(mask)
    print(f"  Basin contains {n_pixels} pixels ({n_pixels * 0.25:.1f} km²)")

    return mask, row_min, row_max, col_min, col_max


# ============================================================================
# FTP DOWNLOAD
# ============================================================================
def list_ftp_files(ftp_host, ftp_dir, pattern=None):
    """List files on the FTP server matching an optional pattern."""
    ftp = ftplib.FTP(ftp_host)
    ftp.login()  # anonymous login
    ftp.cwd(ftp_dir)
    files = ftp.nlst()
    ftp.quit()

    if pattern:
        files = [f for f in files if pattern in f]
    return sorted(files)


def download_ftp_file(ftp_host, ftp_dir, filename, local_dir):
    """Download a single file from FTP."""
    local_path = os.path.join(local_dir, filename)
    if os.path.exists(local_path):
        print(f"  Already downloaded: {filename}")
        return local_path

    print(f"  Downloading: {filename} ...", end=" ", flush=True)
    ftp = ftplib.FTP(ftp_host)
    ftp.login()
    ftp.cwd(ftp_dir)

    with open(local_path, 'wb') as f:
        ftp.retrbinary(f'RETR {filename}', f.write)

    ftp.quit()
    print("done.")
    return local_path


def discover_and_download_files(year, local_dir, tile="h09v04"):
    """
    Discover and download SPIRES HIST netCDF files for the given year and tile.

    The FTP directory structure is explored to find files matching the year
    and tile. File naming may follow patterns like:
        SPIReS_HIST_h09v04_YYYY.nc     (annual file with time dimension)
        SPIReS_h09v04_YYYYDDD.nc       (daily files by DOY)
        or organized in subdirectories by year/tile

    This function adapts to the actual directory structure found.
    """
    os.makedirs(local_dir, exist_ok=True)

    print(f"\nConnecting to FTP: {FTP_HOST}")
    ftp = ftplib.FTP(FTP_HOST)
    ftp.login()

    # Explore the directory structure
    print(f"Exploring: {FTP_DIR}")
    ftp.cwd(FTP_DIR)
    top_items = []
    ftp.retrlines('LIST', top_items.append)

    # Parse directory listing
    dirs = []
    files = []
    for item in top_items:
        parts = item.split()
        name = parts[-1]
        if item.startswith('d'):
            dirs.append(name)
        else:
            files.append(name)

    print(f"  Found {len(dirs)} directories, {len(files)} files at top level")

    downloaded_files = []

    # Strategy 1: Look for files directly matching tile and year
    year_str = str(year)
    matching_files = [f for f in files if tile in f and year_str in f]

    if matching_files:
        print(f"  Found {len(matching_files)} files matching {tile} and {year}")
        for fname in matching_files:
            path = download_ftp_file(FTP_HOST, FTP_DIR, fname, local_dir)
            downloaded_files.append(path)
        ftp.quit()
        return downloaded_files

    # Strategy 2: Check subdirectories (by tile, then year, or vice versa)
    # Try tile subdirectory first
    if tile in dirs:
        tile_dir = f"{FTP_DIR}/{tile}"
        ftp.cwd(tile_dir)
        sub_items = ftp.nlst()

        # Look for year-specific files or subdirectories
        year_files = [f for f in sub_items if year_str in f and '.nc' in f]
        if year_files:
            print(f"  Found {len(year_files)} files in {tile}/ for {year}")
            for fname in year_files:
                path = download_ftp_file(FTP_HOST, tile_dir, fname, local_dir)
                downloaded_files.append(path)
        elif year_str in sub_items:
            # Year subdirectory within tile
            year_dir = f"{tile_dir}/{year_str}"
            ftp.cwd(year_dir)
            nc_files = [f for f in ftp.nlst() if '.nc' in f]
            print(f"  Found {len(nc_files)} files in {tile}/{year}/")
            for fname in nc_files:
                path = download_ftp_file(FTP_HOST, year_dir, fname, local_dir)
                downloaded_files.append(path)
        ftp.quit()
        if downloaded_files:
            return downloaded_files

    # Strategy 3: Check year subdirectory, then tile
    if year_str in dirs:
        year_dir = f"{FTP_DIR}/{year_str}"
        ftp.cwd(year_dir)
        sub_items = ftp.nlst()

        tile_files = [f for f in sub_items if tile in f and '.nc' in f]
        if tile_files:
            print(f"  Found {len(tile_files)} files in {year}/ for {tile}")
            for fname in tile_files:
                path = download_ftp_file(FTP_HOST, year_dir, fname, local_dir)
                downloaded_files.append(path)
        elif tile in sub_items:
            tile_dir = f"{year_dir}/{tile}"
            ftp.cwd(tile_dir)
            nc_files = [f for f in ftp.nlst() if '.nc' in f]
            print(f"  Found {len(nc_files)} files in {year}/{tile}/")
            for fname in nc_files:
                path = download_ftp_file(FTP_HOST, tile_dir, fname, local_dir)
                downloaded_files.append(path)
        ftp.quit()
        if downloaded_files:
            return downloaded_files

    # Strategy 4: Download any netCDF files with tile name
    matching = [f for f in files if tile in f and '.nc' in f.lower()]
    if matching:
        print(f"  Found {len(matching)} tile files at top level")
        for fname in matching:
            path = download_ftp_file(FTP_HOST, FTP_DIR, fname, local_dir)
            downloaded_files.append(path)

    ftp.quit()

    if not downloaded_files:
        print(f"\n  WARNING: Could not find files automatically.")
        print(f"  Directory contents: {dirs[:10]} ... {files[:10]}")
        print(f"  You may need to manually explore the FTP directory and adjust paths.")
        print(f"  Try browsing: ftp://{FTP_HOST}{FTP_DIR}")

    return downloaded_files


# ============================================================================
# NETCDF DATA EXTRACTION
# ============================================================================
def extract_fsca_from_netcdf(nc_files, year, row_min, row_max, col_min, col_max):
    """
    Extract fSCA data from downloaded netCDF files for the specified year.

    Handles both:
      - Annual files with a time dimension
      - Individual daily files

    Returns:
        dates: list of datetime objects
        data:  3D numpy array [time, rows, cols] of fSCA values
    """
    dates = []
    data_list = []

    for nc_path in sorted(nc_files):
        print(f"  Reading: {os.path.basename(nc_path)}")

        ds = nc.Dataset(nc_path, 'r')

        # Discover the fSCA variable name
        fsca_var = None
        possible_names = [
            'snow_fraction', 'fSCA', 'fsca', 'FSCA',
            'snow_fraction_on_ground', 'Snow_Fraction',
            'fractional_snow_cover', 'snow_cover_fraction',
            'viewable_snow_fraction'
        ]

        for name in possible_names:
            if name in ds.variables:
                fsca_var = name
                break

        # If not found, search for any variable with 'snow' and 'frac' in name
        if fsca_var is None:
            for var_name in ds.variables:
                if 'snow' in var_name.lower() and 'frac' in var_name.lower():
                    fsca_var = var_name
                    break

        # Last resort: list all variables for debugging
        if fsca_var is None:
            print(f"    Available variables: {list(ds.variables.keys())}")
            # Try the first 2D or 3D variable that isn't a coordinate
            coord_names = {'x', 'y', 'lat', 'lon', 'latitude', 'longitude',
                          'time', 'crs', 'sinusoidal'}
            for var_name, var in ds.variables.items():
                if var_name.lower() not in coord_names and len(var.dimensions) >= 2:
                    fsca_var = var_name
                    print(f"    Using variable: {fsca_var}")
                    break

        if fsca_var is None:
            print(f"    WARNING: No suitable variable found, skipping file.")
            ds.close()
            continue

        var = ds.variables[fsca_var]
        print(f"    Variable: {fsca_var}, shape: {var.shape}, dims: {var.dimensions}")

        # Check for time dimension
        if len(var.shape) == 3:
            # 3D variable: [time, y, x]
            time_var = None
            for tname in ['time', 'Time', 'date', 'day']:
                if tname in ds.variables:
                    time_var = ds.variables[tname]
                    break

            if time_var is not None:
                # Parse time values
                time_units = getattr(time_var, 'units', None)
                time_calendar = getattr(time_var, 'calendar', 'standard')

                if time_units:
                    time_dates = nc.num2date(time_var[:], time_units,
                                            calendar=time_calendar)
                else:
                    # Assume day-of-year or ordinal dates
                    time_vals = time_var[:]
                    time_dates = []
                    for tv in time_vals:
                        if tv > 100000:  # Likely MATLAB datenum or ordinal
                            dt = datetime(year, 1, 1) + timedelta(days=int(tv) - 1)
                        elif tv > 2000:  # Likely YYYYDDD format
                            dt = datetime.strptime(str(int(tv)), '%Y%j')
                        else:  # Day of year
                            dt = datetime(year, 1, 1) + timedelta(days=int(tv) - 1)
                        time_dates.append(dt)

                # Filter to requested year
                for t_idx, dt in enumerate(time_dates):
                    if hasattr(dt, 'year'):
                        dt_year = dt.year
                    else:
                        dt_year = year

                    if dt_year == year:
                        slice_data = var[t_idx,
                                        row_min:row_max + 1,
                                        col_min:col_max + 1]
                        data_list.append(np.array(slice_data))
                        if hasattr(dt, 'strftime'):
                            dates.append(dt)
                        else:
                            dates.append(datetime(year, 1, 1)
                                        + timedelta(days=t_idx))
            else:
                # No time variable found; assume sequential days
                nt = var.shape[0]
                for t_idx in range(nt):
                    slice_data = var[t_idx,
                                    row_min:row_max + 1,
                                    col_min:col_max + 1]
                    data_list.append(np.array(slice_data))
                    dates.append(datetime(year, 1, 1) + timedelta(days=t_idx))

        elif len(var.shape) == 2:
            # Single 2D image (daily file)
            slice_data = var[row_min:row_max + 1, col_min:col_max + 1]
            data_list.append(np.array(slice_data))

            # Try to extract date from filename
            fname = os.path.basename(nc_path)
            date_found = False
            # Try YYYYDDD pattern
            for i in range(len(fname) - 6):
                chunk = fname[i:i+7]
                if chunk.isdigit() and chunk[:4] == str(year):
                    try:
                        dt = datetime.strptime(chunk, '%Y%j')
                        dates.append(dt)
                        date_found = True
                        break
                    except ValueError:
                        pass
            if not date_found:
                # Try YYYYMMDD
                for i in range(len(fname) - 7):
                    chunk = fname[i:i+8]
                    if chunk.isdigit() and chunk[:4] == str(year):
                        try:
                            dt = datetime.strptime(chunk, '%Y%m%d')
                            dates.append(dt)
                            date_found = True
                            break
                        except ValueError:
                            pass
            if not date_found:
                dates.append(datetime(year, 1, 1)
                            + timedelta(days=len(dates)))

        ds.close()

    if data_list:
        data = np.stack(data_list, axis=0)
    else:
        data = np.array([])

    return dates, data


# ============================================================================
# ANIMATION
# ============================================================================
def create_animation(dates, data, mask, output_path, year,
                     bbox_sin, polygons_sin,
                     row_min, row_max, col_min, col_max,
                     fps=10):
    """
    Create an MP4 animation of daily fSCA clipped to the basin boundary.
    """
    print(f"\nCreating animation: {output_path}")
    print(f"  {len(dates)} frames at {fps} fps")

    # Apply mask: set pixels outside basin to NaN
    mask_3d = np.broadcast_to(mask[np.newaxis, :, :], data.shape)
    data_masked = np.where(mask_3d, data, np.nan)

    # Handle fill/missing values
    data_masked = np.where((data_masked < 0) | (data_masked > 1.0),
                           np.nan, data_masked)

    # If data is in percentage (0-100), convert to fraction
    valid = data_masked[np.isfinite(data_masked)]
    if len(valid) > 0 and np.nanmax(valid) > 1.5:
        data_masked = data_masked / 100.0

    # Setup colormap
    cmap = plt.cm.Blues.copy()
    cmap.set_bad(color='#F5F5DC')  # beige for no-data/outside basin

    # Compute geographic extent for axis labels
    x_left = -R_EARTH * math.pi + H_TILE * TILE_SIZE + col_min * PIXEL_SIZE
    x_right = -R_EARTH * math.pi + H_TILE * TILE_SIZE + (col_max + 1) * PIXEL_SIZE
    y_top = R_EARTH * math.pi / 2 - V_TILE * TILE_SIZE - row_min * PIXEL_SIZE
    y_bottom = R_EARTH * math.pi / 2 - V_TILE * TILE_SIZE - (row_max + 1) * PIXEL_SIZE

    # Convert corners to lat/lon for axis labels
    lat_top, lon_left = sinusoidal_to_latlon(x_left, y_top)
    lat_bottom, lon_right = sinusoidal_to_latlon(x_right, y_bottom)

    # Prepare frames
    frame_dir = os.path.join(os.path.dirname(output_path), 'frames_temp')
    os.makedirs(frame_dir, exist_ok=True)

    writer = imageio.get_writer(output_path, fps=fps, codec='libx264',
                                 quality=8)

    for i, (date, frame_data) in enumerate(zip(dates, data_masked)):
        if i % 30 == 0:
            print(f"  Frame {i+1}/{len(dates)}: {date.strftime('%Y-%m-%d')}")

        fig, ax = plt.subplots(1, 1, figsize=(10, 8))

        im = ax.imshow(frame_data,
                       extent=[lon_left, lon_right, lat_bottom, lat_top],
                       cmap=cmap, vmin=0, vmax=1,
                       interpolation='nearest', aspect='auto')

        # Overlay basin boundary
        for poly in polygons_sin:
            # Convert sinusoidal polygon to lat/lon
            lons_poly = []
            lats_poly = []
            for px, py in poly:
                plat, plon = sinusoidal_to_latlon(px, py)
                lons_poly.append(plon)
                lats_poly.append(plat)
            ax.plot(lons_poly, lats_poly, 'k-', linewidth=1.0)

        ax.set_xlabel('Longitude (°)')
        ax.set_ylabel('Latitude (°)')
        ax.set_title(f'Boise River Basin — Fractional Snow Covered Area\n'
                     f'{date.strftime("%B %d, %Y")}',
                     fontsize=14, fontweight='bold')

        cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label('Fractional Snow Covered Area', fontsize=11)

        # Add mean fSCA text
        valid_data = frame_data[np.isfinite(frame_data)]
        if len(valid_data) > 0:
            mean_fsca = np.nanmean(valid_data)
            ax.text(0.02, 0.02,
                    f'Basin mean fSCA: {mean_fsca:.3f}',
                    transform=ax.transAxes, fontsize=11,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        fig.tight_layout()

        # Save frame to buffer
        fig.canvas.draw()
        frame_image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
        frame_image = frame_image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        writer.append_data(frame_image)

        plt.close(fig)

    writer.close()
    print(f"  Animation saved: {output_path}")

    # Cleanup temp frames
    import shutil
    if os.path.exists(frame_dir):
        shutil.rmtree(frame_dir)


# ============================================================================
# MAIN
# ============================================================================
def main():
    parser = argparse.ArgumentParser(
        description='Download SPIRES fSCA, clip to Boise River Basin, and animate.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python download_clip_animate_fSCA.py --year 2020 --shapefile BRB_outline.shp
    python download_clip_animate_fSCA.py --year 2015 --shapefile BRB_outline.shp --fps 15
    python download_clip_animate_fSCA.py --year 2020 --shapefile BRB_outline.shp --output_dir ./results
        """)

    parser.add_argument('--year', type=int, required=True,
                        help='Year to process (e.g., 2020). Valid range: 2000-2025.')
    parser.add_argument('--shapefile', type=str, required=True,
                        help='Path to the Boise River Basin shapefile (.shp).')
    parser.add_argument('--output_dir', type=str, default='./spires_output',
                        help='Directory for downloads and outputs (default: ./spires_output).')
    parser.add_argument('--fps', type=int, default=10,
                        help='Frames per second for animation (default: 10).')
    parser.add_argument('--tile', type=str, default='h09v04',
                        help='MODIS tile ID (default: h09v04 for Boise River Basin).')

    args = parser.parse_args()

    # Validate year
    if args.year < 2000 or args.year > 2025:
        sys.exit(f"ERROR: Year must be between 2000 and 2025 (got {args.year})")

    # Validate shapefile
    if not os.path.exists(args.shapefile):
        sys.exit(f"ERROR: Shapefile not found: {args.shapefile}")

    print("=" * 70)
    print("  SPIRES fSCA — Download, Clip, and Animate")
    print("=" * 70)
    print(f"  Year:      {args.year}")
    print(f"  Tile:      {args.tile}")
    print(f"  Shapefile: {args.shapefile}")
    print(f"  Output:    {args.output_dir}")
    print("=" * 70)

    # Step 1: Read shapefile and create basin mask
    print("\n[1/4] Reading shapefile and creating basin mask...")
    bbox_sin, polygons_sin = read_shapefile_boundary(args.shapefile)
    mask, row_min, row_max, col_min, col_max = create_basin_mask(
        bbox_sin, polygons_sin, H_TILE, V_TILE)

    # Step 2: Download data
    print(f"\n[2/4] Downloading SPIRES data for {args.year}...")
    download_dir = os.path.join(args.output_dir, 'downloads')
    nc_files = discover_and_download_files(args.year, download_dir, tile=args.tile)

    if not nc_files:
        sys.exit("\nERROR: No data files were downloaded. Check FTP connection and "
                 "directory structure.")

    # Step 3: Extract fSCA
    print(f"\n[3/4] Extracting fSCA data...")
    dates, data = extract_fsca_from_netcdf(
        nc_files, args.year, row_min, row_max, col_min, col_max)

    if len(dates) == 0:
        sys.exit("\nERROR: No fSCA data extracted. Check variable names in netCDF files.")

    print(f"  Extracted {len(dates)} time steps")
    print(f"  Date range: {dates[0].strftime('%Y-%m-%d')} to "
          f"{dates[-1].strftime('%Y-%m-%d')}")

    # Step 4: Create animation
    print(f"\n[4/4] Creating animation...")
    output_mp4 = os.path.join(args.output_dir,
                              f'BRB_fSCA_{args.year}.mp4')
    create_animation(dates, data, mask, output_mp4, args.year,
                     bbox_sin, polygons_sin,
                     row_min, row_max, col_min, col_max,
                     fps=args.fps)

    print("\n" + "=" * 70)
    print("  COMPLETE")
    print(f"  Animation: {output_mp4}")
    print("=" * 70)


if __name__ == '__main__':
    main()
