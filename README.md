# Boise River Basin Fractional Snow Covered Area Visualization

Tools for downloading, clipping, reprojecting, and animating SPIReS HIST V01
fractional snow covered area (fSCA) data over the Boise River Basin, Idaho.

## Data Source

SPIRES HIST V01 from NSIDC:
- FTP: `ftp://dtn.rc.colorado.edu/shares/snow-today/gridded_data/SPIRES_HIST_V01`
- DOI: https://doi.org/10.7265/a3vr-c014
- Citation: Rittger, K., et al. (2025). Historical MODIS/Terra L3 Global Daily
  500m SIN Grid Snow Cover, Snow Albedo, and Snow Surface Properties.
  (SPIRES_HIST, Version 1). NSIDC.

## Files

- `download_clip_animate_fSCA.py` — Python script (standalone)
- `download_clip_animate_fSCA.m` — MATLAB main script (orchestrator)
- `plot_fsca_frame.m` — MATLAB helper for visualizing a single day's fSCA

## Requirements

### Python
```
NOTE: PYTHON IS CURRENTLY NOT WORKING!!!! PLEASE USE MATLAB.
pip install netCDF4 numpy matplotlib pyshp imageio[ffmpeg]
```

### MATLAB
- R2018b+ with Mapping Toolbox (for `shaperead`)

## Usage

### Python
```bash
python download_clip_animate_fSCA.py --year 2020 --shapefile BRB_outline.shp
```

### MATLAB — Quick look at one day
```matlab
plot_fsca_frame('nc_file', 'SPIRES_HIST_h09v04_MOD09GA061_20150201_V1.0.nc', ...
                'shapefile', 'BRB_outline.shp')
```

### MATLAB — Full water year animation
```matlab
download_clip_animate_fSCA(2015, 'BRB_outline.shp')
download_clip_animate_fSCA(2020, 'BRB_outline.shp', 'output_format', 'kmz')
```

## MODIS Tile

The Boise River Basin falls entirely within MODIS sinusoidal tile **h09v04**.

## SNOTEL Stations

Eight SNOTEL stations within/near the basin are plotted:
Bogus Basin (#978), Mores Creek Summit (#637), Graham Guard Sta. (#496),
Jackson Peak (#550), Atlanta Summit (#306), Trinity Mtn. (#830),
Banner Summit (#312), Prairie (#700).
