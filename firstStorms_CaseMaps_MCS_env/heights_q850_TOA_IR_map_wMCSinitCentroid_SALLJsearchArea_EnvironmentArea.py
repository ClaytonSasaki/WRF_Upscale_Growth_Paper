#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on wed Oct 1 2024

@author: crs326

Makes maps for times when MCS initated within the chosen region. Plots wind shear barbs at a chosen height, points that meet the SALLJ criteria, and simulated maximum reflectivity. 

Also plots outline where MCS centroid was searched (black box), MCS centroid (black star), and where SALLJ coverage is calculated (dashed black box)

script based on shear_and_SALLJ_strength_map_single_zoom_larger_dots.py

which was based upon shear_and_SALLJ_strength_map.py, REFD_MAX_and_shear_map.py and SALLJ_spatial_map_zoom_efficient_v3.py
"""

import matplotlib
matplotlib.use('Agg') 

import math
import os

import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature, LAND, OCEAN, COASTLINE, BORDERS, LAKES, RIVERS
from datetime import datetime, timedelta
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
                 cartopy_ylim, latlon_coords, interplevel, ll_to_xy, get_basemap, xy_to_ll, CoordPair, GeoBounds)

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

cmap2 = plt.get_cmap('Purples')
new_cmap2 = truncate_colormap(cmap2, 0.25, 1.0)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap2)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap2)

cdict1 = {'charteruse', 'lime', 'mediumsegreen', 'forestgreen', 'darkgreen'}
new_green_cmap = mc.LinearSegmentedColormap('newGreens', cdict1)

def SatVap(tempc,phase="liquid"):
    """Calculate saturation vapour pressure over liquid water and/or ice.

    INPUTS: 
    tempc: (C)
    phase: ['liquid'],'ice'. If 'liquid', do simple dew point. If 'ice',
    return saturation vapour pressure as follows:

    Tc>=0: es = es_liquid
    Tc <0: es = es_ice


    RETURNS: e_sat  (Pa)

    SOURCE: http://cires.colorado.edu/~voemel/vp.html (#2:
    CIMO guide (WMO 2008), modified to return values in Pa)

    This formulation is chosen because of its appealing simplicity, 
    but it performs very well with respect to the reference forms
    at temperatures above -40 C. At some point I'll implement Goff-Gratch
    (from the same resource).
    """

    over_liquid=6.112*exp(17.67*tempc/(tempc+243.12))*100.
    over_ice=6.112*exp(22.46*tempc/(tempc+272.62))*100.
    # return where(tempc<0,over_ice,over_liquid)

    if phase=="liquid":
    # return 6.112*exp(17.67*tempc/(tempc+243.12))*100.
        return over_liquid
    elif phase=="ice":
    # return 6.112*exp(22.46*tempc/(tempc+272.62))*100.
        return where(tempc<0,over_ice,over_liquid)
    else:
        raise NotImplementedError

Epsilon=0.622         # Epsilon=R_s_da/R_s_v; The ratio of the gas constants

def MixRatio(e,p):
    """Mixing ratio of water vapour
    INPUTS
    e (Pa) Water vapor pressure
    p (Pa) Ambient pressure

    RETURNS
    qv (kg kg^-1) Water vapor mixing ratio`
    """

    return Epsilon*e/(p-e)

Rs_da=287.05          # Specific gas const for dry air, J kg^{-1} K^{-1}
Rs_v=461.51           # Specific gas const for water vapour, J kg^{-1} K^{-1}

#################### VARIABLES TO CHANGE ##########################

MCS_type = 'all_MCS' # 'all_MCS' or 'robustMCS'

MCS_init_location_filter = True

MCS_start_type_filter = True

offset_MCS_and_conditions = False

hours_offset = -1

MCS_init_area = 'large_area1' # 'large_area1', 'large_area2', 'Zhang_area', 'SDC_area1'

# NOTE: SALLJ_search_text variables defined below not currently used in file name so will produce will in new folder but with the same name if EITHER SALLJ_search OPTION BELOW IS CHANGED!!
SALLJ_search_area = '2deg4degOffset1degNFromCentroid' # '2deg4degOffset1degNFromCentroid', '1deg3degBottomCentroid', '60-65W28-30SFixed'
plot_SALLJ_search_area = True

env_search_area = '2.00fromMCScentroid' # '0.75fromMCScentroid', '2.00fromMCScentroid'
plot_env_search_area = True

plot_radar_refl = False
plot_ctt = True
plot_shear_barbs = False
plot_SALLJ_strength = False
plot_q_850 = True
plot_MUCAPE = False
plot_heights = True

######### get MCS initation times and centroid locations #########

if MCS_type == 'all_MCS':
    MCS_tracked_file = 'mcs_tracks_pf_20181015_20190430.nc'
    MCS_file_label = 'MCS'
elif MCS_type == 'robustMCS':
    MCS_tracked_file = 'robust_mcs_tracks_20181015_20190430.nc'
    MCS_file_label = 'robustMCS'
else:
    print('MUST choose valid MCS_type')

# read in MCS tracks
MCS_tracks_file_path = '/home/disk/monsoon/relampago/analysis//mcs_tracks/wrf/stats/%s' %(MCS_tracked_file)
MCS_tracks_ncfile = Dataset(MCS_tracks_file_path,'r')

# gets MCS times in byte format
MCS_datetime_bytes = MCS_tracks_ncfile.variables['datetimestring']

# gets MCS iniation times still in byte format
MCS_datetime_initiation_bytes = MCS_datetime_bytes[:,0,:]

# concatenates the bytes that make up one MCS initation time into a single byte string and decodes to a string
MCS_datetime_initiation_str = np.array([bytes(bytes_for_one_time).decode('UTF-8') for bytes_for_one_time in MCS_datetime_initiation_bytes])

# another way to do the same as the line above with the join function instead
#print([b''.join(row).decode('UTF-8') for row in datetime_bytes])

# gets the MCS center lons and lats (nmaxpf chosen to be 0, NOTE: not sure what nmaxpf is?????????)
MCS_center_lons = MCS_tracks_ncfile.variables['meanlon']
MCS_center_lats = MCS_tracks_ncfile.variables['meanlat']

# get MCS initation center lons and lats
MCS_center_lons_initiation = MCS_center_lons[:,0]
MCS_center_lats_initiation = MCS_center_lats[:,0]

#print(MCS_time_initation_str)
#print(MCS_center_lons_initation)
#print(MCS_center_lats_initation)
#print('dates of MCS init \n', MCS_datetime_initiation_str)
#print('lons of MCS init \n', MCS_center_lons_initiation)
#print('lats of MCS init \n', MCS_center_lats_initiation)

######### create filter MCS tracks by centroid location ######### 

# save MCS initation times and centroid locations of those within defined box
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
        
condition_lat_min = lat_min <= MCS_center_lats_initiation
condition_lat_max = MCS_center_lats_initiation <= lat_max
condition_lon_min = lon_min <= MCS_center_lons_initiation
condition_lon_max = MCS_center_lons_initiation <= lon_max

######### creat filter MCS tracks by MCS start type ######### 

MCS_start_type = np.array(MCS_tracks_ncfile.variables['starttrackresult'])

condition_start_type = MCS_start_type == 10

######### fitler MCS tracks by chosen filters ######### 

if(MCS_init_location_filter == True) and (MCS_start_type_filter == False):

    mask = condition_lat_min & condition_lat_max & condition_lon_min & condition_lon_max
    
    filter_label = '_filtered_init_loc'
    
if(MCS_init_location_filter == False) and (MCS_start_type_filter == True):

    mask = condition_start_type
    
    filter_label = '_filtered_start_type'
    
if(MCS_init_location_filter == True) and (MCS_start_type_filter == True):

    mask = condition_lat_min & condition_lat_max & condition_lon_min & condition_lon_max & condition_start_type
    
    filter_label = '_filtered_init_loc_and_start_type'
    
if(MCS_init_location_filter == False) and (MCS_start_type_filter == False):

    mask = np.full(shape=len(MCS_datetime_initiation_str), fill_value=True, dtype=bool)
    
    filter_label = ''
                                   
MCS_datetime_initiation_str_filtered = MCS_datetime_initiation_str[mask]
MCS_center_lons_initiation_filtered = MCS_center_lons_initiation[mask]
MCS_center_lats_initiation_filtered = MCS_center_lats_initiation[mask]

# number of MCSs that meet chosen conditions
num_MCS = len(MCS_datetime_initiation_str_filtered)
print('# of MCSs that meet chosen conditions: ', num_MCS)

MCS_duration = np.array(MCS_tracks_ncfile.variables['length'])

print('MCS_duration', MCS_duration)
print('len(MCS_duration)', len(MCS_duration))

MCS_duration_filtered = MCS_duration[mask]

MCS_majoraxislength_init = np.array(MCS_tracks_ncfile.variables['majoraxislength'])[:,0]

MCS_majoraxislength_init_2 = np.array(MCS_tracks_ncfile.variables['majoraxislength'])[:,1]

MCS_majoraxislength_growth = MCS_majoraxislength_init_2 - MCS_majoraxislength_init

MCS_majoraxislength_growth_filtered = MCS_majoraxislength_growth[mask]

##### compare MCS times to SALLJ times ###
#
#print('dates of MCS init filtered by location and start type: ', MCS_datetime_init_filt_by_loc_and_start_type_str)
#
#SALLJ_time_list = pickle.load(open("/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/Composite_shear_maps_by_SALLJ_presence/3relaxed_SALLJ_times_%s_%s_%s.dat" %('VMRS', '20181101', '20190430'), "rb"))
#
#SALLJ_time_list_arr = np.array([str(SALLJ_time) for SALLJ_time in SALLJ_time_list])
#
#SALLJ_time_list_arr_dt = [datetime.strptime(date, '%Y-%m-%d %H:%M:%S') for date in SALLJ_time_list_arr]
#
#print('dates of SALLJ at VMRS: ', SALLJ_time_list_arr_dt)

# close the MCS tracks file
MCS_tracks_ncfile.close()



# ----------------- reading in file/plotting terrain -----------------

# Open the NetCDF file
path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'

#pressure to plot
#level = 850

#change the start and end days to be plotted
start_date = 20181022
end_date = 20181022

# 2: high - allows for jets centered higher in the atmosphere, and a bit weaker compared to mean flow
#crit = [2 ,3600, 6000, 23.326, 11.663]
crit = [2 ,3200, 5700, 19.4384, 7.77538]

crit_num = crit[0]
max_search_hgt = crit[1]
min_search_hgt = crit[2]
max_wind_threshold = crit[3]
decrease_to_min_threshold = crit[4]

######################################

# go through list of times/centroids for each MCS to get corresponding environmental conditions
for count, (MCS_datetime, MCS_center_lon, MCS_center_lat) in enumerate(zip(MCS_datetime_initiation_str_filtered, MCS_center_lons_initiation_filtered, MCS_center_lats_initiation_filtered)):

    # get the file necessary
    path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'                  
                                   
    MCS_datetime_dt = datetime.strptime(MCS_datetime, '%Y-%m-%d_%H:%M:')
    
    if offset_MCS_and_conditions == True:
        
        conditions_datetime_dt = MCS_datetime_dt + timedelta(hours = hours_offset)
        
    else: #offset_MCS_and_conditions == False
        
        conditions_datetime_dt = MCS_datetime_dt
    
    MCS_datetime_str = conditions_datetime_dt.strftime('%Y%m%d')

    wrf_file_name = '/wrfout_d01_%s-%s-%s_%s:00:00' %(conditions_datetime_dt.strftime('%Y'), conditions_datetime_dt.strftime('%m'), conditions_datetime_dt.strftime('%d'), conditions_datetime_dt.strftime('%H'))

    wrf_file_path = path_wrf + MCS_datetime_str + wrf_file_name
    
    print('path of MCS within region: ', wrf_file_path)

    # get the netCDF
    wrf_ncfile = Dataset(wrf_file_path,'r')

    # Get the terrain height
    terrain = getvar(wrf_ncfile, 'HGT')

    # Get the latitude and longitude points
    lats, lons = latlon_coords(terrain)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(terrain)

    ############################# read in and plot other variables #############################

    # read in pressure used for multiple variables
    pres = getvar(wrf_ncfile, 'pressure')
    
    if plot_SALLJ_strength == True:
        
        SALLJ_label = '_relaxed_SALLJ_strength'
        
        hght = getvar(wrf_ncfile, 'height', msl=False, units='m') # in m NOTE: this is AGL not MSL!!
        speed, drct = getvar(wrf_ncfile, 'wspd_wdir', units='kt') # in kts

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

        ########### plot SALLJ presence/strength ###########

        SALLJ_strength = axs.contourf(lons, lats, max_wind_SALLJs, levels=np.arange(21,45,2), extend='max', cmap='Purples', zorder=2, transform=crs.PlateCarree())

        fig.subplots_adjust(right=0.8)
        cbar_ax3 = fig.add_axes([0.15, 0.85, 0.7, 0.05])
        fig.colorbar(SALLJ_strength,cax=cbar_ax3, label='SALLJ peak wind (kts)', orientation='horizontal')
        
    else:
        SALLJ_label = ''

    # ------------------------- calculate and plot wind shear ------------------------------
    
    if plot_shear_barbs == True:
        
        shear_barb_label = '_%d_shear_barbs_smooth_level_barbs_and_smooth_sfc_barbs_blue_few' %(level)
        
        u = getvar(wrf_ncfile, 'ua', units='kt') # in kts
        v = getvar(wrf_ncfile, 'va', units='kt') # in kts

        level = 700

        u_level = np.squeeze(interplevel(u, pres, level))
        v_level = np.squeeze(interplevel(v, pres, level))

        u_sfc = np.squeeze(u[3, :, :])
        v_sfc = np.squeeze(v[3, :, :])

        u_diff = u_level - u_sfc
        v_diff = v_level - v_sfc

        mag_diff = np.sqrt(u_diff**2 + v_diff**2)

        cmap=plt.cm.jet
        bounds = [0, 10, 20, 30, 40, 50]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    #                bounds = [0, 10, 20, 30, 40, 50]
    #                norm = matplotlib.colors.BoundaryNorm(bounds, new_cmap2.N)

        smooth_u_level = smooth2d(u_level, 10)
        smooth_v_level = smooth2d(v_level, 10)

        smooth_u_sfc = smooth2d(u_sfc, 10)
        smooth_v_sfc = smooth2d(v_sfc, 10)

        ########### plot wind and wind shear barbs ###########

        shear_barbs = axs.barbs(to_np(lons[::barb_interval, ::barb_interval]), to_np(lats[::barb_interval, ::barb_interval]),  to_np(u_diff[::barb_interval, ::barb_interval]), to_np(v_diff[::barb_interval, ::barb_interval]), to_np(mag_diff[::barb_interval, ::barb_interval]), transform=crs.PlateCarree(), cmap=cmap, norm=norm, length=10, linewidth=4.5, zorder=7)

        level_barbs = axs.barbs(to_np(lons[::barb_interval, ::barb_interval]), to_np(lats[::barb_interval, ::barb_interval]),  to_np(smooth_u_level[::barb_interval, ::barb_interval]), to_np(smooth_v_level[::barb_interval, ::barb_interval]), transform=crs.PlateCarree(), length=10, linewidth=4.5, zorder=6, color='lightslategray')

        sfc_barbs = axs.barbs(to_np(lons[::barb_interval, ::barb_interval]), to_np(lats[::barb_interval, ::barb_interval]),  to_np(smooth_u_sfc[::barb_interval, ::barb_interval]), to_np(smooth_v_sfc[::barb_interval, ::barb_interval]), transform=crs.PlateCarree(), length=10, linewidth=4.5, zorder=6, color='black')
        
    else:
        shear_barb_label = ''

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
    
    
    # plot MCS initation box
    axs.plot([lon_min, lon_max], [lat_max, lat_max], 'k-', lw=6, transform=crs.PlateCarree(), zorder=6)
    axs.plot([lon_min, lon_min], [lat_min, lat_max], 'k-', lw=6, transform=crs.PlateCarree(), zorder=6)
    axs.plot([lon_min, lon_max], [lat_min, lat_min], 'k-', lw=6, transform=crs.PlateCarree(), zorder=6)
    axs.plot([lon_max, lon_max], [lat_min, lat_max], 'k-', lw=6, transform=crs.PlateCarree(), zorder=6)
    
    # plot MCS centroid
    axs.scatter(MCS_center_lon, MCS_center_lat, transform=crs.PlateCarree(), color='black', marker='*', zorder=7, s=1600)
    
    # plot SALLJ search area 
    if plot_SALLJ_search_area == True:
        
        SALLJ_search_text = '__SALLJarea_%s' %(SALLJ_search_area)
    
        if SALLJ_search_area == '1deg3degBottomCentroid':

            lat_bottom_left_SALLJ = MCS_center_lat
            lon_bottom_left_SALLJ = MCS_center_lon - 1.50
            lat_top_right_SALLJ = MCS_center_lat + 1.00 
            lon_top_right_SALLJ = MCS_center_lon + 1.50

        elif SALLJ_search_area == '2deg4degOffset1degNFromCentroid':

            lat_bottom_left_SALLJ = MCS_center_lat + 1.00
            lon_bottom_left_SALLJ = MCS_center_lon - 2.00
            lat_top_right_SALLJ = MCS_center_lat + 3.00 
            lon_top_right_SALLJ = MCS_center_lon + 2.00
            
        elif SALLJ_search_area == '60-65W28-30SFixed':
            
            lat_bottom_left_SALLJ = -30
            lon_bottom_left_SALLJ = -65
            lat_top_right_SALLJ = -28
            lon_top_right_SALLJ = -60
            
        else:
            print('Please enter valid SALLJ search area') # will throw error
        
        
        # plot MCS initation box
        axs.plot([lon_bottom_left_SALLJ, lon_top_right_SALLJ], [lat_top_right_SALLJ, lat_top_right_SALLJ], 'k--', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_bottom_left_SALLJ, lon_bottom_left_SALLJ], [lat_bottom_left_SALLJ, lat_top_right_SALLJ], 'k--', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_bottom_left_SALLJ, lon_top_right_SALLJ], [lat_bottom_left_SALLJ, lat_bottom_left_SALLJ], 'k--', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_top_right_SALLJ, lon_top_right_SALLJ], [lat_bottom_left_SALLJ, lat_top_right_SALLJ], 'k--', lw=6, transform=crs.PlateCarree(), zorder=6)
            
    else:
        
        SALLJ_search_text = ''
        
        
        
    # plot environmental search area    
    if plot_env_search_area == True:
        
        env_search_text = '__EnvArea_%s' %(env_search_area)
    
        if env_search_area == '0.75fromMCScentroid':

            # get lats/lons of region based on centroid
            lat_bottom_left_env = MCS_center_lat - 0.75 
            lon_bottom_left_env = MCS_center_lon - 0.75

            lat_top_right_env = MCS_center_lat + 0.75 
            lon_top_right_env = MCS_center_lon + 0.75
            
        if env_search_area == '2.00fromMCScentroid':

            # get lats/lons of region based on centroid
            lat_bottom_left_env = MCS_center_lat - 2.00 
            lon_bottom_left_env = MCS_center_lon - 2.00

            lat_top_right_env = MCS_center_lat + 2.00 
            lon_top_right_env = MCS_center_lon + 2.00

        else:
            print('Please enter valid env search area') # will throw error
        
        
        # plot MCS initation box
        axs.plot([lon_bottom_left_env, lon_top_right_env], [lat_top_right_env, lat_top_right_env], 'c--', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_bottom_left_env, lon_bottom_left_env], [lat_bottom_left_env, lat_top_right_env], 'c--', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_bottom_left_env, lon_top_right_env], [lat_bottom_left_env, lat_bottom_left_env], 'c--', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_top_right_env, lon_top_right_env], [lat_bottom_left_env, lat_top_right_env], 'c--', lw=6, transform=crs.PlateCarree(), zorder=6)
            
    else:
        
        env_search_text = ''

    ########### plot relfectivity ###########
    
    if plot_radar_refl == True:
        radar_refl_max = getvar(wrf_ncfile, 'REFD_MAX')

        # smooth relfectivity variable

        smooth_radar_refl_max = smooth2d(radar_refl_max, 10)

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
    
    if plot_ctt == True:
            
            ctt_label = '_cloud_top_temp_241_225'
        
            ctt = getvar(wrf_ncfile, 'ctt', units='K')
            
            cloud_top_temp_241 = axs.contour(to_np(lons), to_np(lats), to_np(ctt), levels=[241], linestyles='solid', transform=crs.PlateCarree(), colors=['blue'], linewidths=np.arange(5,1,-1), zorder=3)
            
            cloud_top_temp_225 = axs.contour(to_np(lons), to_np(lats), to_np(ctt), levels=[225], linestyles='solid', transform=crs.PlateCarree(), colors=['darkblue'], linewidths=np.arange(5,1,-1), zorder=3)
            
    else:
        ctt_label = ''
        
    if plot_q_850 == True:
        
        q_850_label = '_q_850'
        
        ######### calculate q at 850 hPa #########
        level = 850

        # read in the variables
        temp_c_subset = getvar(wrf_ncfile, 'temp', units='degC')
        RH_subset = getvar(wrf_ncfile, 'rh') # in kg/kg

        # interpolate to the set level
        temp_c_subset_level = np.squeeze(interplevel(temp_c_subset, pres, level))

        RH_subset_level = np.squeeze(interplevel(RH_subset, pres, level))

        e_sat_subset_level = SatVap(temp_c_subset_level)

        e_subset_level = (RH_subset_level/100.) * e_sat_subset_level

        w_subset_level = MixRatio(e_subset_level,level*100.)
        #w = (e*Rs_da) / (Rs_v*(pres*100. - e))
        q_subset_level = (w_subset_level / (w_subset_level+1))*1000 #g/kg
        
        q_850 = axs.contourf(lons, lats, q_subset_level, extend='max', cmap='Greens', alpha=0.7, zorder=2, transform=crs.PlateCarree())
        
        fig.subplots_adjust(right=0.8)
        cbar_ax3 = fig.add_axes([0.15, 0.85, 0.7, 0.05])
        fig.colorbar(q_850,cax=cbar_ax3, label='850hPa q (g/kg)', orientation='horizontal')
        
    else:
        q_850_label = ''

    if plot_MUCAPE == True:
        
        MUCAPE_label = '_MUCAPE'

        MCAPE_original = getvar(wrf_ncfile, 'cape_2d')[0][int(centroid_xy[1]),int(centroid_xy[0])].values

        MCAPE = axs.contourf(lons, lats, MCAPE_original, extend='max', cmap='Reds', alpha=0.7, zorder=2, transform=crs.PlateCarree())

    else:
        MUCAPE_label = ''
        
    if plot_heights == True:
        
        level = 500
        
        hght = getvar(wrf_ncfile, 'height', units='m') # in m NOTE: this is MSL not AGL!!!
        pres = getvar(wrf_ncfile, 'pressure')
        
        hght_level = np.squeeze(interplevel(hght, pres, level))
        
        smooth_hght_level = smooth2d(hght_level, 10)
        
        heights = axs.contour(to_np(lons), to_np(lats), to_np(smooth_hght_level), levels=np.arange(5340,6000,60), extend='both', transform=crs.PlateCarree(), linewidths=6, cmap=new_cmap2, zorder=3)
        
        axs.clabel(heights, inline=1, fontsize=22, fmt='%.0f') 

    #######################################################

    #plt.title("%s" % (area))

    axs.scatter(-64.212, -31.298, s=240, color='red', transform=crs.PlateCarree(), zorder=5)
    axs.scatter(-63.726, -29.906, s=240, color='red', transform=crs.PlateCarree(), zorder=5)

    #print 'saving'

    #plt.suptitle("Max Reflectivity and sfc-%d wind shear: %sZ" % (level, file_name[11:24]))

    #plt.savefig('/home/disk/meso-home/crs326/Documents/Research/WRF_Paper/PNNL_WRF/WRF_MAPS/SALLJ_strength_and_shear_maps/PNNL_WRF_%sZ_relaxed_SALLJ_strength_%d_shear_barbs_smooth_level_barbs_and_smooth_sfc_barbs_%s.png' %(file_name[11:24], level, 'many_area'), dpi=200)
    
    general_outpath = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/firstStorms_CaseMaps_MCS_env/Figures_'
    
    specific_outpath = '%sarea_%s%s%s' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)
    
    plt.savefig(general_outpath + specific_outpath + '/%s_%s_init_area_WRF_%sZ%s%s%s%s%s_%s_larger_dots.png' %(MCS_init_area, MCS_file_label, conditions_datetime_dt.strftime('%Y%m%d%H'), ctt_label, q_850_label, MUCAPE_label, SALLJ_label, shear_barb_label, '1_panel'), dpi=600)

    print('saved')

    plt.close()
    wrf_ncfile.close()

