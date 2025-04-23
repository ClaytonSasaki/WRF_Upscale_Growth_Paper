#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Oct 24 2024

@author: crs326

Description:
------------
This script plots a base map with terrain and geographic features over South America. 
It also highlights a specific MCS initiation region and optionally includes additional 
elements such as radar locations.
"""

# Use non-interactive backend for environments without display (e.g., remote servers)
import matplotlib
matplotlib.use('Agg') 

import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature, LAND, OCEAN, COASTLINE, BORDERS, LAKES, RIVERS
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from wrf import (to_np, getvar, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, CoordPair, GeoBounds)

# Set plot-wide font size
matplotlib.rcParams.update({'font.size': 25})

# ----------------- USER-DEFINED SETTINGS ------------------

filter_label = '_filtered_init_loc_and_start_type'
MCS_type = 'all_MCS'               # Options: 'all_MCS', 'robustMCS'
MCS_init_area = 'large_area1'      # Choose region for MCS initiation
env_offset_MCS = False             # If True, offset environment data from MCS start time
hours_offset = -1                  # Hours offset for environmental sampling
seperation_plot = False            # Whether to plot MCS growth rate separation
growth_type = 'slow'               # MCS growth type: 'slow', 'rapid'.
plot_dir_shear_label = False       # Label shear direction (optional)
shear_type = '0-3km'               # Shear depth (e.g., '0-3km')
events_removed = True              # Whether to exclude certain MCS events

# ----------------------------------------------------------

# Set file label for type of MCS
if MCS_type == 'all_MCS':
    MCS_file_label = 'MCS'
elif MCS_type == 'robustMCS':
    MCS_file_label = 'robustMCS'
else:
    raise ValueError("Invalid MCS_type. Choose 'all_MCS' or 'robustMCS'.")
    
if env_offset_MCS == True:
    offset_label = '_%dhrPrior' %(abs(hours_offset))
else:
    offset_label = ''

# Define MCS initiation area bounding box
if MCS_init_area == 'large_area1':
    lon_min = -66.0
    lon_max = -57.0
    lat_min = -36.0
    lat_max = -28.0

elif MCS_init_area == 'large_area2':
    lon_min = -66.0
    lon_max = -57.0
    lat_min = -36.0
    lat_max = -29.5
    
elif MCS_init_area == 'Zhang_area':
    lon_min = -70.0
    lon_max = -59.0
    lat_min = -36.0
    lat_max = -26.5
    
elif MCS_init_area == 'SDC_area1':
    lon_min = -66.25
    lon_max = -61.1
    lat_min = -34.0
    lat_max = -29.5

else:
    print('Please add the matching lats and lons to search!')
    
# ---------- Load WRF Terrain Data ----------

# get the file path to plot terrain
path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'
wrf_file_path = path_wrf + '20181101/wrfout_d01_2018-11-01_00:00:00' # date doesn't matter as just plotting terrain

# get the netCDF
wrf_ncfile = Dataset(wrf_file_path,'r')

# Get the terrain height
terrain = getvar(wrf_ncfile, 'HGT')

# Get the latitude and longitude points
lats, lons = latlon_coords(terrain)

# Get the cartopy mapping object
cart_proj = get_cartopy(terrain)

# ---------- Set Plot Area and Figure ----------

fig, axs = plt.subplots(1, 1, figsize=(32, 32), subplot_kw={'projection': cart_proj})

axs.set_xlim(cartopy_xlim(terrain))
axs.set_ylim(cartopy_ylim(terrain))

# ---------- Add Map Features and Terrain ----------

axs.add_feature(OCEAN, zorder=2)
axs.add_feature(LAKES, alpha=0.5, zorder=2)
axs.add_feature(RIVERS, zorder=2)
axs.add_feature(BORDERS, edgecolor='gray', zorder=2)

# Make  filled contours for the terrain
terrain_plot = axs.contourf(to_np(lons), to_np(lats), to_np(terrain), levels=[0, 500, 1000], colors=['papayawhip', 'navajowhite','burlywood'], extend='max', transform=crs.PlateCarree(), zorder=1)

# Add a color bar
#                fig.subplots_adjust(right=0.8)
#                cbar_ax1 = fig.add_axes([0.15, 0.05, 0.25, 0.05])
#                fig.colorbar(terrain_plot,cax=cbar_ax1, label='Elevation (m)', orientation='horizontal')

#cbar = plt.colorbar(terrain_plot, ax=axs, shrink=.2, orientation='horizontal')
#cbar.set_label('Elevation (m)')

# Add the gridlines
gl = axs.gridlines(crs=crs.PlateCarree(), linewidth=1, color='gray', alpha=0.7, draw_labels=True, linestyle='--', zorder=6)
gl.top_labels = False
gl.right_labels = False
gl.left_labels = True
gl.bottom_labels = True

# ---------- MCS Initiation Box and Reference Points ----------

# plot MCS initation box
axs.plot([lon_min, lon_max], [lat_max, lat_max], 'k-', lw=6, transform=crs.PlateCarree(), zorder=6)
axs.plot([lon_min, lon_min], [lat_min, lat_max], 'k-', lw=6, transform=crs.PlateCarree(), zorder=6)
axs.plot([lon_min, lon_max], [lat_min, lat_min], 'k-', lw=6, transform=crs.PlateCarree(), zorder=6)
axs.plot([lon_max, lon_max], [lat_min, lat_max], 'k-', lw=6, transform=crs.PlateCarree(), zorder=6)

axs.scatter(-64.19, -31.44, s=400, color='red', transform=crs.PlateCarree(), zorder=5) # Cordoba radar
#axs.scatter(-64.212, -31.298, s=240, color='red', transform=crs.PlateCarree(), zorder=5)
#axs.scatter(-63.726, -29.906, s=240, color='red', transform=crs.PlateCarree(), zorder=5)

# ---------- Finalize and Save Figure ----------

wrf_ncfile.close()

general_outpath = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/MCSinit_InitialGrowthRate_shear_direction_MCS_centroids_plus0.5deg'

specific_outpath = '/data/'

plt.savefig(general_outpath + '/blank_map_whole_domain_with_region_withOnlyCOR.png', dpi=600)

print('saved')

plt.close()
    