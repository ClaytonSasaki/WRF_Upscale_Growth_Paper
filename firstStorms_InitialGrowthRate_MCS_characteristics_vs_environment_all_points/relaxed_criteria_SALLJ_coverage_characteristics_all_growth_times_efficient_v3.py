#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Oct 8 2024

@author: crs326

Searches for SALLJ anytime during the growth period (first storms to MCS initation)

-----

Gets MCS characteristics (e.g., duration, area) and times for MCSs in a defined region.

Calculates environmental characteristics (shear, SALLJ) for *all points* in a region around MCS initiation location centroid. Saves MCS characteristics with each point at the MCS initation time

If percentage area coverage of the SALLJ meets a criteria, a SALLJ is found and the average SALLJ characteristics are used. Saves SALLJ characteristics from the whole area used as a characteristic for each point at the MCS initation time.

"""

# Use non-interactive backend for environments without display (e.g., remote servers)
import matplotlib
matplotlib.use('Agg') 

from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from numpy import exp,where,ma,cos,sin,pi,amax,amin
import pickle
import os
import xarray

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, vinterp, ll_to_xy, get_basemap, xy_to_ll, CoordPair, GeoBounds, srhel, g_cape, cape_2d, cape_3d)

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

# ----------------- USER-DEFINED SETTINGS ------------------

MCS_type = 'all_MCS'               # Options: 'all_MCS', 'robustMCS'
MCS_init_area = 'large_area1'      # Choose region for MCS initiation
MCS_init_location_filter = True
MCS_start_type_filter = True
offset_MCS_and_conditions = False
hours_offset = -1

SALLJ_search_area = '2deg4degOffset1degNFromCentroid' # '2deg4degOffset1degNFromCentroid', '1deg3degBottomCentroid', '60-65W28-30SFixed'

env_search_area = '2.00fromMCScentroid' # '0.75fromMCScentroid', '2.00fromMCScentroid'

# ----------------------------------------------------------

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

MCS_status = np.array(MCS_tracks_ncfile.variables['mcs_status'])

### find indicies of when MCS is first found (when MCS_status is first 1 or 2) ###
    
# Initialize a boolean array of the same shape, filled with False
first_MCS_bol_array = np.zeros_like(MCS_status, dtype=bool)

# Initialize a boolean array of the same shape, filled with False
time_0_bol_array = np.zeros_like(MCS_status, dtype=bool)

# Initalize array with the length of the number of tracks
MCS_growth_stage_time_length = np.full(MCS_status.shape[0], np.nan)

# Find the first occurrence of `1` in each row
for i in range(MCS_status.shape[0]):
    
    ## Check if `1` or `2` (MCS or squall line) is in the row ##
    
    contains = np.isin(MCS_status[i], [1, 2]) # create condition if 1 or 2 is the in the row (MCS track)
    
    if contains.any(): # check if 1 or 2 is the in the row (MCS track)
        
        # If MCS/Squall line found, get first tracked time (by setting first index to True)
        time_0_bol_array[i, 0] = True
        
        # Get the index where MCS/Squall line and set that index to True
        first_one_index = np.where(contains)[0][0]
        first_MCS_bol_array[i, first_one_index] = True
        
        # The index is uses to calculate the duration as each time step is 30 min
        MCS_growth_stage_time_length[i] = first_one_index/2
        
## Possible alternative way without loops
#mask = (MCS_status == 1) | (MCS_status == 2)
#
## Find the index of the first occurrence of 1 or 2 in each row
#first_occurrence_indices = np.argmax(mask, axis=1)
#
## Create an array of False with the same shape as arr
#bool_array = np.zeros_like(MCS_status, dtype=bool)
#
## Set the first occurrence of 1 or 2 in each row to True using advanced indexing
#bool_array[np.arange(MCS_status.shape[0]), first_occurrence_indices] = True
#
## Make sure True is only set where there was a 1 or 2 in the original array
#bool_array &= mask
#
#print(bool_array)

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

ccs_area_init = np.array(MCS_tracks_ncfile.variables['ccs_area'])[:,0]
MCS_ccs_area_init = np.array(MCS_tracks_ncfile.variables['ccs_area'])[first_MCS_bol_array]

MCS_ccs_area_growth = MCS_ccs_area_init - ccs_area_init
MCS_ccs_area_growth_filtered = MCS_ccs_area_growth[mask]

MCS_growth_stage_time_length_filtered = MCS_growth_stage_time_length[mask]

print('median', np.nanmedian(MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered))

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
MCS_tracks_ncfile.close

# create arrays that will hold MCS characterstics reshaped to hold data repeated for all points in area
MCS_duration_all_points = []
MCS_majoraxislength_growth_all_points = []
MCS_ccs_area_growth_all_points = []
MCS_growth_stage_time_length_all_points = []

SALLJ_count = 0
low_cov_SALLJ_count = 0
med_cov_SALLJ_count = 0
high_cov_SALLJ_count = 0

prop_area_SALLJ_all_growth_times = np.full((12, len(MCS_center_lons_initiation_filtered)), np.nan)


# go through list of times/centroids for each MCS to get corresponding environmental conditions
for count, (MCS_datetime, MCS_center_lon, MCS_center_lat, MCS_growth_time) in enumerate(zip(MCS_datetime_initiation_str_filtered, MCS_center_lons_initiation_filtered, MCS_center_lats_initiation_filtered, MCS_growth_stage_time_length)):

    # get the file necessary
    path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'                  
                                   
    MCS_datetime_dt = datetime.strptime(MCS_datetime, '%Y-%m-%d_%H:%M:')
    
    if offset_MCS_and_conditions == True:
        
        conditions_datetime_dt = MCS_datetime_dt + timedelta(hours = hours_offset)
        
    else: #offset_MCS_and_conditions == False
        
        conditions_datetime_dt = MCS_datetime_dt
        
    print('MCS_growth_time', MCS_growth_time)
        
    for growth_time_step in range(int(MCS_growth_time*2+1)):
        
        conditions_datetime_dt = MCS_datetime_dt + timedelta(hours = growth_time_step/2)

        MCS_datetime_str = conditions_datetime_dt.strftime('%Y%m%d')

        wrf_file_name = '/wrfout_d01_%s-%s-%s_%s:00:00' %(conditions_datetime_dt.strftime('%Y'), conditions_datetime_dt.strftime('%m'), conditions_datetime_dt.strftime('%d'), conditions_datetime_dt.strftime('%H'))

        wrf_file_path = path_wrf + MCS_datetime_str + wrf_file_name

        print('path of MCS within region: ', wrf_file_path)
        print('growth time: ', growth_time_step/2)

        # get the netCDF
        wrf_ncfile = Dataset(wrf_file_path,'r')

        ########## find if SALLJ is present and if so calculate characteristics #########

        # for finding LLJs #

        # 2: high - allows for jets centered higher in the atmosphere, and a bit weaker compared to mean flow
        crit = [2 ,3200, 5700, 19.4384, 7.77538] # relaxed to 10 m/s max and 4 m/s decrease

        crit_num = crit[0]
        max_search_hgt = crit[1]
        min_search_hgt = crit[2]
        max_wind_threshold = crit[3]
        decrease_to_min_threshold = crit[4]

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

        # get xy for regions
        bottom_left_xy = ll_to_xy(wrf_ncfile, lat_bottom_left_SALLJ, lon_bottom_left_SALLJ)
        top_right_xy = ll_to_xy(wrf_ncfile, lat_top_right_SALLJ, lon_top_right_SALLJ) 
    
        ############################# read in some variables #############################
        #print('reading in variables')
        # read in the variables subsetted to a specified area
        pres_subset = getvar(wrf_ncfile, 'pressure')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])]
        hght_subset = getvar(wrf_ncfile, 'height', msl=False, units='m')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in m NOTE: this is AGL not MSL!!
        radar_refl_max = getvar(wrf_ncfile, 'REFD_MAX')[int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])]
    
        speed, drct = getvar(wrf_ncfile, 'wspd_wdir', units='kt') # in kts
    
        speed_subset = speed[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in kts
    
        drct_subset = drct[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])]
    
    #                # get wind speed and direction at chosen height 1 and calculate wind shear
    #                wspd_chosen1 = interplevel(speed_subset, hght_subset, chosen_height1)
    #                wdir_chosen1 = interplevel(drct_subset, hght_subset, chosen_height1)
    #                
    #                udiff = wspd_chosen1*np.cos(np.radians(270-wdir_chosen1)) - speed_subset[3, :, :]*np.cos(np.radians(270-drct_subset[3, :, :]))
    #                vdiff = wspd_chosen1*np.sin(np.radians(270-wdir_chosen1)) - speed_subset[3, :, :]*np.sin(np.radians(270-drct_subset[3, :, :]))
    #
    #                shear_chosen1 = np.sqrt(udiff**2 + vdiff**2)
    #                
    #                print('wspd_chosen1', wspd_chosen1[30, 30].values)
    #                print('wdir_chosen1', wdir_chosen1[30, 30].values)
    #                print('wspd_sfc', speed_subset[3, 30, 30].values)
    #                print('wdir_sfc', drct_subset[3, 30, 30].values)
    #                print('shear_chosen1', shear_chosen1[30, 30].values)
    #                
    #                # get wind speed and direction at chosen height 2 and calculate wind shear
    #                wspd_chosen2 = interplevel(speed_subset, hght_subset, chosen_height2)
    #                wdir_chosen2 = interplevel(drct_subset, hght_subset, chosen_height2)
    #                
    #                udiff = wspd_chosen2*np.cos(np.radians(270-wdir_chosen2)) - speed_subset[3, :, :]*np.cos(np.radians(270-drct_subset[3, :, :]))
    #                vdiff = wspd_chosen2*np.sin(np.radians(270-wdir_chosen2)) - speed_subset[3, :, :]*np.sin(np.radians(270-drct_subset[3, :, :]))
    #
    #                shear_chosen2 = np.sqrt(udiff**2 + vdiff**2)
    #                
    #                # get wind speed and direction at chosen height 3 and calculate wind shear
    #                wspd_chosen3 = interplevel(speed_subset, hght_subset, chosen_height3)
    #                wdir_chosen3 = interplevel(drct_subset, hght_subset, chosen_height3)
    #                
    #                udiff = wspd_chosen3*np.cos(np.radians(270-wdir_chosen3)) - speed_subset[3, :, :]*np.cos(np.radians(270-drct_subset[3, :, :]))
    #                vdiff = wspd_chosen3*np.sin(np.radians(270-wdir_chosen3)) - speed_subset[3, :, :]*np.sin(np.radians(270-drct_subset[3, :, :]))
    #
    #                shear_chosen3 = np.sqrt(udiff**2 + vdiff**2)
    #                
    #                # get wind speed and direction at chosen height 6 and calculate wind shear
    #                wspd_chosen6 = interplevel(speed_subset, hght_subset, chosen_height6)
    #                wdir_chosen6 = interplevel(drct_subset, hght_subset, chosen_height6)
    #                
    #                udiff = wspd_chosen6*np.cos(np.radians(270-wdir_chosen6)) - speed_subset[3, :, :]*np.cos(np.radians(270-drct_subset[3, :, :]))
    #                vdiff = wspd_chosen6*np.sin(np.radians(270-wdir_chosen6)) - speed_subset[3, :, :]*np.sin(np.radians(270-drct_subset[3, :, :]))
    #
    #                shear_chosen6 = np.sqrt(udiff**2 + vdiff**2)
    #                
    #                # get wind shear 2-6 km
    #                
    #                udiff = wspd_chosen6*np.cos(np.radians(270-wdir_chosen6)) - wspd_chosen2*np.cos(np.radians(270-wdir_chosen2))
    #                vdiff = wspd_chosen6*np.sin(np.radians(270-wdir_chosen6)) - wspd_chosen2*np.sin(np.radians(270-wdir_chosen2))
    #
    #                shear_chosen2_6 = np.sqrt(udiff**2 + vdiff**2)
    
        ############################# find SALLJ #############################
    
        # interpolate variables to 250 meter vertical spacing up to the search heigth for the wind minimum 
        interp_levels_min = np.arange(0,min_search_hgt+250,250)
    
        pres_below_min = interplevel(pres_subset, hght_subset, interp_levels_min)
        hght_below_min = interplevel(hght_subset, hght_subset, interp_levels_min)
        speed_below_min = interplevel(speed_subset, hght_subset, interp_levels_min)
        drct_below_min = interplevel(drct_subset, hght_subset, interp_levels_min)
    
        # interpolate variables to 250 meter vertical spacing up to the search heigth for the wind maximum
        interp_levels_max = np.arange(0,max_search_hgt+250,250)
    
        pres_below_max = interplevel(pres_subset, hght_subset, interp_levels_max)
        hght_below_max = interplevel(hght_subset, hght_subset, interp_levels_max)
        speed_below_max = interplevel(speed_subset, hght_subset, interp_levels_max)
        drct_below_max = interplevel(drct_subset, hght_subset, interp_levels_max)
    
        #print(pres_below_max.shape)
        #print(hght_below_max.shape)
    
        #print('pres_below_max', pres_below_max.values)
        #print('hght_below_max', hght_below_max.values)
        #print('speed_below_max', speed_below_max.values)
    
        ################ with xarray ###################
    
        # get max wind
        max_wind = speed_below_max.max(dim='level')
    
        # get height of max wind
        level_max_wind = hght_below_max.isel(level=speed_below_max.argmax('level'))
    #                level_max_wind = speed_below_max.idxmax('level')
    #                level_max_wind = speed_below_max[np.where(np.equal(speed_below_max,max_wind))].get_index('level')
        #level_max_wind2 = speed_below_max.idxmax(dim='level')
    
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
        hght_below_min_masked_below_max = hght_below_min.where(hght_minus_level_max_wind > 0., np.nan)
    
        speed_below_min_masked_above_max = speed_below_min.where(hght_minus_level_max_wind < 0., np.nan)
        hght_below_min_masked_above_max = hght_below_min.where(hght_minus_level_max_wind < 0., np.nan)
    
        #print('speed_below_min_masked_below_max', speed_below_min_masked_below_max)
    
        # get min wind above max
        min_wind = speed_below_min_masked_below_max.min(dim='level')
    
        min_wind_below_max = speed_below_min_masked_above_max.min(dim='level')
    
        # get index of min above max
    
        #min_wind_index = speed_below_min_masked_below_max.idxmin(dim='level')
    
        #min_wind_index_below_max = speed_below_min_masked_above_max.idxmin(dim='level')
    
        level_min_wind = hght_below_min_masked_below_max.isel(level=speed_below_min_masked_below_max.argmin('level'))
    
        #level_min_wind_below_max = hght_below_min_masked_above_max.isel(level=speed_below_min_masked_above_max.argmin('level'))
    
        #print('min_wind', min_wind.values)
    
        # checks if max_wind meets threshold and keeps value if it does meet the threshold
        max_wind_meeting_threshold = max_wind.where(max_wind > max_wind_threshold, np.nan)
    
        #print('max_wind_meeting_threshold', max_wind_meeting_threshold.values)
    
        #print('max_wind_meeting_threshold', max_wind_meeting_threshold)
    
        # calculates decrease to min wind
        decrease_to_min = max_wind_meeting_threshold - min_wind
    
        decrease_to_min_below_max = max_wind_meeting_threshold - min_wind_below_max
    
        #print('decrease_to_min', decrease_to_min)
    
        # checks if decrease_to_min meets threshold and keeps value if it does meet the threshold
        decrease_to_min_meeting_threshold = decrease_to_min.where(decrease_to_min > decrease_to_min_threshold, np.nan)
    
        #print('decrease_to_min_meeting_threshold', decrease_to_min_meeting_threshold)
    
        # checks to see if the values met the other criteria (was not replaced with nan) and if it does leave the value
        drct_at_max_wind_meeting_threshold = drct_at_max_wind.where(np.isnan(decrease_to_min_meeting_threshold) == False, np.nan)
    
        #print('drct_at_max_wind_meeting_threshold', drct_at_max_wind_meeting_threshold.values)
    
        # checks to see if wind at max_wind is from a northerly direction and keeps value if it is
        drct_at_max_wind_meeting_threshold = drct_at_max_wind_meeting_threshold.where((drct_at_max_wind_meeting_threshold <= 45) | (drct_at_max_wind_meeting_threshold >= 315), np.nan)
    
        #print('drct_at_max_wind_meeting_threshold', drct_at_max_wind_meeting_threshold)
        
        num_points_SALLJ = np.count_nonzero(~np.isnan(drct_at_max_wind_meeting_threshold))
        num_total_points = drct_at_max_wind_meeting_threshold.size
        
        #print('num points', num_total_points)
        
        proporation_SALLJ = num_points_SALLJ/num_total_points
        
        print('proporation coverage: ', proporation_SALLJ)
        
        prop_area_SALLJ_all_growth_times[growth_time_step, count] = proporation_SALLJ
        
#        proporation_SALLJ_one_time = np.empty_like(q_subset_level)
#        proporation_SALLJ_one_time[:] = proporation_SALLJ
#        prop_area_SALLJ_all_points.append(proporation_SALLJ_one_time)
#        
#        print('proporation_SALLJ_one_time.shape', proporation_SALLJ_one_time.shape)
#        print('proporation_SALLJ_one_time', proporation_SALLJ_one_time)
#        
#        if proporation_SALLJ >= 0.12:
#            
#            print('found SALLJ')
#            
#            SALLJ_count = SALLJ_count + 1
#            
#            ### get SALLJ characteristics for points where the SALLJ criteria is met
#    
#            # get max_wind for points that meet SALLJ criteria
#            max_wind_SALLJs = max_wind.where(np.isnan(drct_at_max_wind_meeting_threshold) == False, np.nan)
#            
#            median_SALLJ_max_wind = np.nanmedian(max_wind_SALLJs)
#            
#            median_SALLJ_max_wind_one_time = np.empty_like(q_subset_level)
#            median_SALLJ_max_wind_one_time[:] = median_SALLJ_max_wind
#            median_SALLJ_max_wind_all_points.append(median_SALLJ_max_wind_one_time)
#            
#            print('median_SALLJ_max_wind_one_time.shape', median_SALLJ_max_wind_one_time.shape)
#            print('median_SALLJ_max_wind_one_time', median_SALLJ_max_wind_one_time)
#            
#            # get level of max_wind for points that meet SALLJ criteria
#            level_max_wind_SALLJs = level_max_wind.where(np.isnan(drct_at_max_wind_meeting_threshold) == False, np.nan)
#            
#            median_SALLJ_height = np.nanmedian(level_max_wind_SALLJs)
#            
#            median_SALLJ_height_one_time = np.empty_like(q_subset_level)
#            median_SALLJ_height_one_time[:] = median_SALLJ_height
#            median_SALLJ_height_all_points.append(median_SALLJ_height_one_time)
#            
#            print('median_SALLJ_height_one_time.shape', median_SALLJ_height_one_time.shape)
#            print('median_SALLJ_height_one_time', median_SALLJ_height_one_time)
#            
#            #print('Median SALLJ wind: %s   Median SALLJ level: %s' %(median_SALLJ_max_wind, median_SALLJ_level))
#            
#            if proporation_SALLJ < 0.3:
#                
#                low_cov_SALLJ_count = low_cov_SALLJ_count + 1
#                
#                print('low_coverage_SALLJ')
#            elif proporation_SALLJ <= 0.3 and proporation_SALLJ < 0.6:
#                
#                print('med_coverage_SALLJ')
#                
#                med_cov_SALLJ_count = med_cov_SALLJ_count + 1
#                
#            elif proporation_SALLJ >= 0.6:
#                
#                print('high_coverage_SALLJ')
#                
#                high_cov_SALLJ_count = high_cov_SALLJ_count + 1
#            
#        else:
#            print('no SALLJ found')
#            
#            # fill SALLJ characteristics with nans if SALLJ coverage threshold not met
#            
#            median_SALLJ_max_wind_one_time = np.empty_like(q_subset_level)
#            median_SALLJ_max_wind_one_time[:] = np.nan
#            median_SALLJ_max_wind_all_points.append(median_SALLJ_max_wind_one_time)
#            
#            print('median_SALLJ_max_wind_one_time.shape', median_SALLJ_max_wind_one_time.shape)
#            print('median_SALLJ_max_wind_one_time', median_SALLJ_max_wind_one_time)
#            
#            median_SALLJ_height_one_time = np.empty_like(q_subset_level)
#            median_SALLJ_height_one_time[:] = np.nan
#            median_SALLJ_height_all_points.append(median_SALLJ_height_one_time)
#            
#            print('median_SALLJ_height_one_time.shape', median_SALLJ_height_one_time.shape)
#            print('median_SALLJ_height_one_time', median_SALLJ_height_one_time)
#            
#    #        SALLJ_status = 'no_SALLJ'
#            #pass
        
    
        wrf_ncfile.close()
    
#print('SALLJ_count', SALLJ_count)
#print('low_cov_SALLJ_count', low_cov_SALLJ_count)
#print('med_cov_SALLJ_count', med_cov_SALLJ_count)
#print('high_cov_SALLJ_count', high_cov_SALLJ_count)
    
general_outpath = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/firstStorms_InitialGrowthRate_MCS_characteristics_vs_environment_all_points//'

specific_outpath = '%sarea_%s%s/data/' %(MCS_file_label, MCS_init_area, SALLJ_search_text)

pickle.dump(prop_area_SALLJ_all_growth_times, open(general_outpath + specific_outpath + "%s_prop_area_SALLJ_all_growth_times%s_%s_%s_.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area), "wb"))

    