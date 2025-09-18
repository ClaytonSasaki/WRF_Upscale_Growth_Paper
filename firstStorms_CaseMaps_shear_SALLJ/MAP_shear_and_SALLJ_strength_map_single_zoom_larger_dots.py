#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on wed Oct 10 15:00:00 2022

@author: crs326

Makes maps that plots wind shear at a chosen height and shows points that meet the SALLJ criteria and the height at which the criteria is met from the PNNL WRF. The same as shear_and_SALLJ_strength_map.py but with one close zoom plotted and updated for readability.

Based upon REFD_MAX_and_shear_map.py and SALLJ_spatial_map_zoom_efficient_v3.py
"""

import matplotlib
matplotlib.use('Agg') 

import math
import os

import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature, LAND, OCEAN, COASTLINE, BORDERS, LAKES, RIVERS
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.colors as mc
import metpy.calc as mpcalc
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from numpy import exp, where, ma, cos, sin, pi, amax, amin
import pandas as pd
from scipy import stats
from scipy.interpolate import interp1d
#import time
import xarray

# ------------------------- cmaps ------------------------------

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, interpz3d, ll_to_xy, get_basemap, xy_to_ll, CoordPair, GeoBounds)

import matplotlib.colors as colors

#start_time = time.time()
matplotlib.rcParams.update({'font.size': 47})

# used for contourf of terrain 
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap = plt.get_cmap('terrain')
new_cmap = truncate_colormap(cmap, 0.21, 1)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap)

# used for SALLJ
def truncate_colormap2(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap2 = plt.get_cmap('Greys')
new_cmap2 = truncate_colormap(cmap2, 0.2, 1.0)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap2)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap2)

# used for contour of wind
def truncate_colormap2(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap2 = plt.get_cmap('Greens')
new_cmap2 = truncate_colormap(cmap2, 0.1, 1.0)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap2)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap2)

cdict1 = {'charteruse', 'lime', 'mediumsegreen', 'forestgreen', 'darkgreen'}
new_green_cmap = mc.LinearSegmentedColormap('newGreens', cdict1)

# ----------------- reading in file/plotting terrain -----------------

# Open the NetCDF file
path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'

#pressure to plot
#level = 850

#change the start and end days to be plotted
start_date = 20181129
end_date = 20181129

# 2: high - allows for jets centered higher in the atmosphere, and a bit weaker compared to mean flow
#crit = [2 ,3600, 6000, 23.326, 11.663]
crit = [2 ,3200, 5700, 19.4384, 7.77538]

crit_num = crit[0]
max_search_hgt = crit[1]
min_search_hgt = crit[2]
max_wind_threshold = crit[3]
decrease_to_min_threshold = crit[4]

######################################


# get dates from times to only run for correct day folders
input_start_day = pd.to_datetime(str(start_date), format='%Y%m%d', errors='ignore')
input_end_day = pd.to_datetime(str(end_date), format='%Y%m%d', errors='ignore')


# loop through files (times)
for folder_date in sorted(os.listdir(path_wrf)):

    # convert the folder name to datetime object
    folder_date_dt = pd.to_datetime(folder_date, format='%Y%m%d', errors='ignore')

    # only check folder if it is named with a date
    if folder_date.startswith('20') and (folder_date_dt >= input_start_day) and (folder_date_dt <= input_end_day):

        folder_path = os.path.join(path_wrf, folder_date)

        sorted_hour_files = sorted(os.listdir(folder_path))

        # go through hourly files for the chosen date
        for hourly_file in sorted_hour_files:

            if hourly_file.startswith('wrfout'):

                # for each new hour set count of SALLJ points to 0
                count_points_SALLJ_hour = 0
                count_total_hour = 0

                file_path = os.path.join(folder_path, hourly_file)

                file_hour = int(hourly_file[22:24])

                # get the netCDF
                ncfile = Dataset(file_path,'r')
                file_name = hourly_file 

                print(str(file_name[16:24]))

                # Get the terrain height
                terrain = getvar(ncfile, 'HGT')

                # Get the latitude and longitude points
                lats, lons = latlon_coords(terrain)

                # Get the cartopy mapping object
                cart_proj = get_cartopy(terrain)

                ############################# read in some other variables #############################
                
                # use for whole area
                radar_refl_max = getvar(ncfile, 'REFD_MAX')
                pres = getvar(ncfile, 'pressure')
                hght = getvar(ncfile, 'height', msl=False, units='m') # in m NOTE: this is AGL not MSL!!
                u = getvar(ncfile, 'ua', units='kt') # in kts
                v = getvar(ncfile, 'va', units='kt') # in kts
                speed, drct = getvar(ncfile, 'wspd_wdir', units='kt') # in kts

                # ------------------------- smooth relfectivity variable ------------------------------
                
                smooth_radar_refl_max = smooth2d(radar_refl_max, 10)

                # ------------------------- calculate SALLJ presence/strength ------------------------------------

                # interpolate variables to 250 meter vertical spacing up to the search heigth for the wind minimum 
                interp_levels_min = np.arange(0,min_search_hgt+250,250)

                pres_below_min = interplevel(pres, hght, interp_levels_min)
                hght_below_min = interplevel(hght, hght, interp_levels_min)
                speed_below_min = interplevel(speed, hght, interp_levels_min)
                drct_below_min = interplevel(drct, hght, interp_levels_min)

                # interpolate variables to 250 meter vertical spacing up to the search heigth for the wind maximum
                interp_levels_max = np.arange(0,max_search_hgt+250,250)

                pres_below_max = interplevel(pres, hght, interp_levels_max)
                hght_below_max = interplevel(hght, hght, interp_levels_max)
                speed_below_max = interplevel(speed, hght, interp_levels_max)
                drct_below_max = interplevel(drct, hght, interp_levels_max)

                ################ with xarray ###################

                # get max wind
                max_wind = speed_below_max.max(dim='level')

                #print('max_wind', max_wind[point_to_chceck_xy[1],point_to_chceck_xy[0]].values)

                # get height of max wind
                level_max_wind = hght_below_max.isel(level=speed_below_max.argmax('level'))
            #                level_max_wind = speed_below_max.idxmax('level')
            #                level_max_wind = speed_below_max[np.where(np.equal(speed_below_max,max_wind))].get_index('level')

                #print('level_max_wind', level_max_wind.values)

                # get pressure at max wind
                pres_max_wind = pres_below_max.isel(level=speed_below_max.argmax('level'))

                # get direction at max wind
                drct_at_max_wind = drct_below_max.isel(level=speed_below_max.argmax('level'))

                #print('drct_at_max_wind', drct_at_max_wind)

                # set dim of level of max wind array to 'level' (height)
                level_max_wind_subtract = xarray.DataArray(level_max_wind, dims=['south_north', 'west_east'])

                #print('level_max_wind_subtract', level_max_wind_subtract)

                #print('hght_below_max', hght_below_max.values)

                # subtract the heights of max wind from the heights
                hght_minus_level_max_wind = hght_below_min - level_max_wind_subtract

                #print('hght_minus_level_max_wind', hght_minus_level_max_wind)

                speed_below_min_masked_below_max = speed_below_min.where(hght_minus_level_max_wind > 0., np.nan) 

                #print('speed_below_min_masked_below_max', speed_below_min_masked_below_max)

                # get min wind above max
                min_wind = speed_below_min_masked_below_max.min(dim='level')

                #print('min_wind', min_wind[point_to_chceck_xy[1],point_to_chceck_xy[0]].values)

                # checks if max_wind meets threshold and keeps value if it does meet the threshold
                max_wind_meeting_threshold = max_wind.where(max_wind > max_wind_threshold, np.nan)

                #print('max_wind_meeting_threshold', max_wind_meeting_threshold[point_to_chceck_xy[1],point_to_chceck_xy[0]].values)

                #print('max_wind_meeting_threshold', max_wind_meeting_threshold)

                # calculates decrease to min wind
                decrease_to_min = max_wind_meeting_threshold - min_wind

                #print('decrease_to_min', decrease_to_min)

                # checks if decrease_to_min meets threshold and keeps value if it does meet the threshold
                decrease_to_min_meeting_threshold = decrease_to_min.where(decrease_to_min > decrease_to_min_threshold, np.nan)

                #print('decrease_to_min_meeting_threshold', decrease_to_min_meeting_threshold)

                # checks to see if the values met the other criteria (was not replaced with nan) and if it does leave the value
                drct_at_max_wind_meeting_threshold = drct_at_max_wind.where(np.isnan(decrease_to_min_meeting_threshold) == False, np.nan)

                #print('drct_at_max_wind_meeting_threshold', drct_at_max_wind_meeting_threshold[point_to_chceck_xy[1],point_to_chceck_xy[0]].values)

                # checks to see if wind at max_wind is from a northerly direction and keeps value if it is
                drct_at_max_wind_meeting_threshold = drct_at_max_wind_meeting_threshold.where((drct_at_max_wind_meeting_threshold <= 45) | (drct_at_max_wind_meeting_threshold >= 315), np.nan)

                #print('drct_at_max_wind_meeting_threshold', drct_at_max_wind_meeting_threshold)

                # get pressure of max_wind of points that meet SALLJ criteria
                max_wind_SALLJs = max_wind.where(np.isnan(drct_at_max_wind_meeting_threshold) == False, np.nan)

                #print('max_wind_SALLJs', max_wind_SALLJs[point_to_chceck_xy[1],point_to_chceck_xy[0]].values)
                
                # ------------------------- calculate wind shear ------------------------------
                
#                level = 700
#                
#                u_level = np.squeeze(interplevel(u, pres, level))
#                v_level = np.squeeze(interplevel(v, pres, level))
#                
#                u_sfc = np.squeeze(u[3, :, :])
#                v_sfc = np.squeeze(v[3, :, :])
#                
#                u_diff = u_level - u_sfc
#                v_diff = v_level - v_sfc

                level2 = 6000
                
                u_level2 = np.squeeze(interplevel(u, hght, level2))
                v_level2 = np.squeeze(interplevel(v, hght, level2))
                
                level1 = 2000
                
                u_level1 = np.squeeze(interplevel(u, hght, level1))
                v_level1 = np.squeeze(interplevel(v, hght, level1))
                
                u_diff = u_level2 - u_level1
                v_diff = v_level2 - v_level1
                
                mag_diff = np.sqrt(u_diff**2 + v_diff**2)
                
                cmap=plt.cm.jet
                bounds = [0, 10, 20, 30, 40, 50]
                norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
#                bounds = [0, 10, 20, 30, 40, 50]
#                norm = matplotlib.colors.BoundaryNorm(bounds, new_cmap2.N)

#                smooth_u_level = smooth2d(u_level, 10)
#                smooth_v_level = smooth2d(v_level, 10)
#            
#                smooth_u_sfc = smooth2d(u_sfc, 10)
#                smooth_v_sfc = smooth2d(v_sfc, 10)

                smooth_u_level2 = smooth2d(u_level2, 10)
                smooth_v_level2 = smooth2d(v_level2, 10)

                smooth_u_level1 = smooth2d(u_level1, 10)
                smooth_v_level1 = smooth2d(v_level1, 10)
            
                # Set the map bound
                
#                area = ['SDC_extra_zoom', 'SDC_zoom', 'broad_zoom']
#                lat_bottom_left = [-33.0, -33.0, -38.0]
#                lon_bottom_left = [-65.5, -66.0, -70.0]
#                
#                lat_top_right = [-30.5, -29.0, -26.0]
#                lon_top_right = [-64.0, -63.0, -54.0]
#
#                barb_interval = [5, 10, 25]
#
#                fig, axs = plt.subplots(2, 2, figsize=(20, 12), subplot_kw={'projection': cart_proj})
                
                area = 'broad_zoom'
                lat_bottom_left = -38.0
                lon_bottom_left = -70.0
                
                lat_top_right = -26.0
                lon_top_right = -54.0

                barb_interval = 35

                fig, axs = plt.subplots(1, 1, figsize=(32, 32), subplot_kw={'projection': cart_proj})


                bounds = GeoBounds(CoordPair(lat=lat_bottom_left, lon=lon_bottom_left),
                                       CoordPair(lat=lat_top_right, lon=lon_top_right))

                axs.set_xlim(cartopy_xlim(terrain, geobounds=bounds))
                axs.set_ylim(cartopy_ylim(terrain, geobounds=bounds))

                ########### plot features, terrain, and gridlines ###########

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
                gl = axs.gridlines(crs=crs.PlateCarree(), linewidth=1, color='gray', alpha=0.5, draw_labels=True, linestyle='--', zorder=6)
                gl.top_labels = False
                gl.right_labels = False
                gl.left_labels = True
                gl.bottom_labels = True

                ########### plot relfectivity ###########

#                    refl_max = axs.contour(to_np(lons), to_np(lats), to_np(smooth_radar_refl_max), levels=np.arange(-10,80,10), linestyles='dashed', extend='max', transform=crs.PlateCarree(), cmap='gist_ncar', zorder=3)

                # blue
#                refl_max = axs.contour(to_np(lons), to_np(lats), to_np(smooth_radar_refl_max), levels=np.arange(0,60,10), linestyles='solid', transform=crs.PlateCarree(), colors='blue', linewidths=np.arange(7,1,-1), zorder=3)

                # blue_few
#                refl_max = axs.contour(to_np(lons), to_np(lats), to_np(smooth_radar_refl_max), levels=np.arange(0,60,20), linestyles='solid', transform=crs.PlateCarree(), colors='blue', linewidths=np.arange(6,.9,-1.7), zorder=3)
    
                # colors
                refl_max = axs.contour(to_np(lons), to_np(lats), to_np(smooth_radar_refl_max), levels=np.arange(0,60,15), linestyles='solid', transform=crs.PlateCarree(), colors=['blue','green','orange','red','darkviolet'], linewidths=np.arange(6,1,-.5), zorder=3)

                # Add a color bar
#                fig.subplots_adjust(right=0.8)
#                cbar_ax2 = fig.add_axes([0.15, 0.14, 0.2, 0.05])
#                fig.colorbar(refl_max,cax=cbar_ax2, label='Max Reflectivity (dBZ)', orientation='horizontal')

                ########### plot SALLJ presence/strength ###########

                SALLJ_strength = axs.contourf(lons, lats, max_wind_SALLJs, levels=np.arange(21,45,2), extend='max', cmap='Purples', zorder=2, transform=crs.PlateCarree())

                fig.subplots_adjust(right=0.8)
                cbar_ax3 = fig.add_axes([0.15, 0.85, 0.7, 0.05])
                fig.colorbar(SALLJ_strength,cax=cbar_ax3, label='SALLJ peak wind (kts)', orientation='horizontal')

                ########### plot wind and wind shear barbs ###########

                shear_barbs = axs.barbs(to_np(lons[::barb_interval, ::barb_interval]), to_np(lats[::barb_interval, ::barb_interval]),  to_np(u_diff[::barb_interval, ::barb_interval]), to_np(v_diff[::barb_interval, ::barb_interval]), to_np(mag_diff[::barb_interval, ::barb_interval]), transform=crs.PlateCarree(), cmap=cmap, norm=norm, length=10, linewidth=4.5, zorder=7)

#                level_barbs = axs.barbs(to_np(lons[::barb_interval, ::barb_interval]), to_np(lats[::barb_interval, ::barb_interval]),  to_np(smooth_u_level[::barb_interval, ::barb_interval]), to_np(smooth_v_level[::barb_interval, ::barb_interval]), transform=crs.PlateCarree(), length=10, linewidth=4.5, zorder=6, color='lightslategray')
#
#                sfc_barbs = axs.barbs(to_np(lons[::barb_interval, ::barb_interval]), to_np(lats[::barb_interval, ::barb_interval]),  to_np(smooth_u_sfc[::barb_interval, ::barb_interval]), to_np(smooth_v_sfc[::barb_interval, ::barb_interval]), transform=crs.PlateCarree(), length=10, linewidth=4.5, zorder=6, color='black')
                
                level2_barbs = axs.barbs(to_np(lons[::barb_interval, ::barb_interval]), to_np(lats[::barb_interval, ::barb_interval]),  to_np(smooth_u_level2[::barb_interval, ::barb_interval]), to_np(smooth_v_level2[::barb_interval, ::barb_interval]), transform=crs.PlateCarree(), length=10, linewidth=4.5, zorder=6, color='lightslategray')
                
                level1_barbs = axs.barbs(to_np(lons[::barb_interval, ::barb_interval]), to_np(lats[::barb_interval, ::barb_interval]),  to_np(smooth_u_level1[::barb_interval, ::barb_interval]), to_np(smooth_v_level1[::barb_interval, ::barb_interval]), transform=crs.PlateCarree(), length=10, linewidth=4.5, zorder=6, color='black')
                #######################################################

                #plt.title("%s" % (area))

                axs.scatter(-64.212, -31.298, s=240, color='red', transform=crs.PlateCarree(), zorder=5)
                axs.scatter(-63.726, -29.906, s=240, color='red', transform=crs.PlateCarree(), zorder=5)

                #print 'saving'
                    
                #plt.suptitle("Max Reflectivity and sfc-%d wind shear: %sZ" % (level, file_name[11:24]))

                #plt.savefig('/home/disk/meso-home/crs326/Documents/Research/WRF_Paper/PNNL_WRF/WRF_MAPS/SALLJ_strength_and_shear_maps/PNNL_WRF_%sZ_relaxed_SALLJ_strength_%d_shear_barbs_smooth_level_barbs_and_smooth_sfc_barbs_%s.png' %(file_name[11:24], level, 'many_area'), dpi=200)
#                plt.savefig('/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/firstStorms_CaseMaps/6PNNL_WRF_%sZ_relaxed_SALLJ_strength_%d_shear_barbs_smooth_level_barbs_and_smooth_sfc_barbs_blue_few_%s_larger_dots.png' %(file_name[11:24], level, '1_panel'), dpi=600)
                plt.savefig('/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/firstStorms_CaseMaps_shear_SALLJ/6PNNL_WRF_%sZ_relaxed_SALLJ_strength_%d-%d_shear_barbs_smooth_level_barbs_and_smooth_sfc_barbs_blue_few_%s_larger_dots.png' %(file_name[11:24], level1, level2, '1_panel'), dpi=600)

                print('saved')

                plt.close()
                ncfile.close()

