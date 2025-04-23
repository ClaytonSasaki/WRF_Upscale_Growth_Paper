#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Oct 22 2024

@author: crs326

This script identifies MCS initiation times and locations,
and calculates environmental characteristics (like vertical wind shear) *only* at *MCS initiation centroids*.

It works with tracked MCS data and WRF model output, filtering based on defined criteria,
and saves results as pickled files for further analysis.

"""

# Use non-interactive backend for environments without display (e.g., remote servers)
import matplotlib
matplotlib.use('Agg') 

from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
import pickle
import os
from scipy.interpolate import interp1d
from wrf import (getvar, ll_to_xy)

def get_wind_direction(u, v):
    """Calculates the wind direction in degrees from u and v components.

    Args:
        u (float): East-West wind component (positive is eastward)
        v (float): North-South wind component (positive is northward)

    Returns:
        float: Wind direction in degrees (0 is North, 90 is East, 180 is South, 270 is West)
    """

    # calculates the angle in radians between the positive x-axis and the line connecting the origin to the point (u, v), then converts the angle from radians to degrees, and finally adjusts the angle to represent the direction from which the wind is blowing (meteorological convention)
    direction = np.degrees(np.arctan2(u, v)) + 180
    
    # ensures the direction is within the range of 0 to 360 degrees
    return direction % 360

# ----------------- USER-DEFINED SETTINGS ------------------

MCS_type = 'all_MCS'               # Options: 'all_MCS', 'robustMCS'
MCS_init_area = 'large_area1'  # Options: 'large_area1', 'large_area2', 'Zhang_area', 'SDC_area1'
offset_MCS_and_conditions = False            # If True, offset environment data from MCS start time
hours_offset = 1 
MCS_init_location_filter = True
MCS_start_type_filter = True

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
    
    #print('MCS_status[i]', MCS_status[i])
    
    contains = np.isin(MCS_status[i], [1, 2]) # create condition if 1 or 2 is the in the row (MCS track)
    
    #print('contains', contains)
    
    if contains.any(): # check if 1 or 2 is the in the row (MCS track)
        
        # If MCS/Squall line found, get first tracked time (by setting first index to True)
        time_0_bol_array[i, 0] = True
        
        # Get the index where MCS/Squall line and set that index to True
        first_one_index = np.where(contains)[0][0]
        #print('np.where(contains)', np.where(contains))
        #print('first_one_index', first_one_index)
        first_MCS_bol_array[i, first_one_index] = True
        #print('first_MCS_bol_array[i, :]', first_MCS_bol_array[i, :])
        
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

first_MCS_bol_array_3d = np.repeat(first_MCS_bol_array[:, :, np.newaxis], repeats=17, axis=2)

# gets MCS iniation times still in byte format
MCS_datetime_initiation_bytes = np.array(MCS_datetime_bytes)[first_MCS_bol_array_3d]

# boolean mask flattened array so reshaping to orignal shape
MCS_datetime_initiation_bytes = MCS_datetime_initiation_bytes.reshape(697, 17)

print(MCS_datetime_initiation_bytes.shape)

# concatenates the bytes that make up one MCS initation time into a single byte string and decodes to a string
MCS_datetime_initiation_str = np.array([bytes(bytes_for_one_time).decode('UTF-8') for bytes_for_one_time in MCS_datetime_initiation_bytes])

print(MCS_datetime_initiation_str.shape)

# another way to do the same as the line above with the join function instead
#print([b''.join(row).decode('UTF-8') for row in datetime_bytes])

# gets the MCS center lons and lats (nmaxpf chosen to be 0, NOTE: not sure what nmaxpf is?????????)
MCS_center_lons = MCS_tracks_ncfile.variables['meanlon']
MCS_center_lats = MCS_tracks_ncfile.variables['meanlat']

# get MCS initation center lons and lats
MCSfirstStorms_center_lons_initiation = MCS_center_lons[:,0]
MCSfirstStorms_center_lats_initiation = MCS_center_lats[:,0]

# get MCS initation center lons and lats
MCS_center_lons_initiation = np.array(MCS_center_lons)[first_MCS_bol_array]
MCS_center_lats_initiation = np.array(MCS_center_lats)[first_MCS_bol_array]

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
        
condition_lat_min = lat_min <= MCSfirstStorms_center_lats_initiation
condition_lat_max = MCSfirstStorms_center_lats_initiation <= lat_max
condition_lon_min = lon_min <= MCSfirstStorms_center_lons_initiation
condition_lon_max = MCSfirstStorms_center_lons_initiation <= lon_max

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

print(mask)

print(MCS_datetime_initiation_str.shape)
print(MCS_datetime_initiation_str_filtered.shape)

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

print('ccs_area_init', ccs_area_init)
print('MCS_ccs_area_init', MCS_ccs_area_init)

MCS_ccs_area_growth = MCS_ccs_area_init - ccs_area_init

MCS_ccs_area_growth_filtered = MCS_ccs_area_growth[mask]

print('MCS_ccs_area_growth_filtered', MCS_ccs_area_growth_filtered)

MCS_growth_stage_time_length_filtered = MCS_growth_stage_time_length[mask]

print('MCS_growth_stage_time_length_filtered', MCS_growth_stage_time_length_filtered)

print('rate', MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered)
print('median', np.nanmedian(MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered))

# close the MCS tracks file
MCS_tracks_ncfile.close

shearDir_0_3km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
shearDir_0_6km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
shearDir_2_6km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
MCS_dt = []
MCS_centroid_elevation = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)

# go through list of times/centroids for each MCS to get corresponding environmental conditions
for count, (MCS_datetime, MCS_center_lon, MCS_center_lat) in enumerate(zip(MCS_datetime_initiation_str_filtered, MCS_center_lons_initiation_filtered, MCS_center_lats_initiation_filtered)):

    # get the file necessary
    path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'                  
                                   
    MCS_datetime_dt = datetime.strptime(MCS_datetime, '%Y-%m-%d_%H:%M:')
    
    MCS_dt.append(MCS_datetime_dt)
    
    if offset_MCS_and_conditions == True:
        
        conditions_datetime_dt = MCS_datetime_dt + timedelta(hours = hours_offset)
        
        if hours_offset < 0:
            offset_label = '_%dhrPrior' %(abs(hours_offset))
        elif hours_offset > 0:
            offset_label = '_%dhrAfter' %(abs(hours_offset))
        else:
            print('Enter valid hours_offset')
        
    else: #offset_MCS_and_conditions == False
        
        conditions_datetime_dt = MCS_datetime_dt
        offset_label = ''
    
    MCS_datetime_str = conditions_datetime_dt.strftime('%Y%m%d')

    wrf_file_name = '/wrfout_d01_%s-%s-%s_%s:00:00' %(conditions_datetime_dt.strftime('%Y'), conditions_datetime_dt.strftime('%m'), conditions_datetime_dt.strftime('%d'), conditions_datetime_dt.strftime('%H'))

    wrf_file_path = path_wrf + MCS_datetime_str + wrf_file_name
    
    print('path of MCS within region: ', wrf_file_path)

    # get the netCDF
    wrf_ncfile = Dataset(wrf_file_path,'r')
    
    # get xy for centroid
    centroid_xy = ll_to_xy(wrf_ncfile, MCS_center_lat+0.5, MCS_center_lon)
    
    # Get the terrain height at MCS centroid
    terrain_centroid = getvar(wrf_ncfile, 'HGT')[int(centroid_xy[1]),int(centroid_xy[0])].values
    
    MCS_centroid_elevation[count] = terrain_centroid

    # read in the variables 
    hght_original = getvar(wrf_ncfile, 'height', msl=False, units='m')[:,int(centroid_xy[1]),int(centroid_xy[0])].values # in m NOTE: this is AGL not MSL!!
    u_original = getvar(wrf_ncfile, 'ua', units='kt')[:,int(centroid_xy[1]),int(centroid_xy[0])].values  # in kts
    v_original = getvar(wrf_ncfile, 'va', units='kt')[:,int(centroid_xy[1]),int(centroid_xy[0])].values # in kts
    
    u2km = interp1d(u_original, hght_original, np.asarray([2000.])).values
    v2km = interp1d(v_original, hght_original, np.asarray([2000.])).values

    u3km = interp1d(u_original, hght_original, np.asarray([3000.])).values
    v3km = interp1d(v_original, hght_original, np.asarray([3000.])).values
    
    u6km = interp1d(u_original, hght_original, np.asarray([6000.])).values
    v6km = interp1d(v_original, hght_original, np.asarray([6000.])).values
    
    u_sfc = u_original[3]
    v_sfc = v_original[3]
    
    #print('v_sfc', v_sfc)
    
    # calculate 0-3 km shear at all points
    udiff_0_3 = u3km - u_sfc
    vdiff_0_3 = v3km - v_sfc
    
    dir_shear_0_3km = get_wind_direction(udiff_0_3, vdiff_0_3)
    
    shearDir_0_3km[count] = dir_shear_0_3km[0]
    
    # calculate 0-6 km shear at all points
    udiff_0_6 = u6km - u_sfc
    vdiff_0_6 = v6km - v_sfc
    
    dir_shear_0_6km = get_wind_direction(udiff_0_6, vdiff_0_6)
    
    shearDir_0_6km[count] = dir_shear_0_6km[0]
    
    # calculate 2-6 km shear at all points
    udiff_2_6 = u6km - u2km
    vdiff_2_6 = v6km - v2km
    
    dir_shear_2_6km = get_wind_direction(udiff_2_6, vdiff_2_6)
    
    shearDir_2_6km[count] = dir_shear_2_6km[0]

    wrf_ncfile.close()
    
general_outpath = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/MCSinit_InitialGrowthRate_shear_direction_MCS_centroids_plus0.5deg'

specific_outpath = '/data/'

pickle.dump(MCS_dt, open(general_outpath + specific_outpath + "%s%s_MCS_dt%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
pickle.dump(MCS_centroid_elevation, open(general_outpath + specific_outpath + "%s%s_MCS_centroid_elevation%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
pickle.dump(MCS_center_lons_initiation_filtered, open(general_outpath + specific_outpath + "%s%s_MCS_center_lons%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
pickle.dump(MCS_center_lats_initiation_filtered, open(general_outpath + specific_outpath + "%s%s_MCS_center_lats%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
pickle.dump(MCS_ccs_area_growth_filtered, open(general_outpath + specific_outpath + "%s%s_ccs_area_growth%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
pickle.dump(MCS_growth_stage_time_length_filtered, open(general_outpath + specific_outpath + "%s%s_growth_stage_time_length%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
pickle.dump(shearDir_0_3km, open(general_outpath + specific_outpath + "%s%s_shearDir_0_3km_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
pickle.dump(shearDir_0_6km, open(general_outpath + specific_outpath + "%s%s_shearDir_0_6km_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
pickle.dump(shearDir_2_6km, open(general_outpath + specific_outpath + "%s%s_shearDir_2_6km_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))