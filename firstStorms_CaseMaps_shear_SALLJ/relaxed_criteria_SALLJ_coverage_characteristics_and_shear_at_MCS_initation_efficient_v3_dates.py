#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sep 18 2023

@author: crs326

Gets MCS characteristics (e.g., duration, area) and times for MCSs in a defined region.

Finds average environmental characteristics (shear, SALLJ) in a region around MCS initiation location centroid. 

If percentage area coverage of the SALLJ meets a criteria, a SALLJ is found and the average SALLJ characteristics are used.


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
from scipy.interpolate import interp1d
import xarray

from wrf import (getvar, interplevel, ll_to_xy)

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

for index in range(len(MCS_datetime_initiation_str_filtered)):

    print(index, MCS_datetime_initiation_str_filtered[index])

#MCS_duration = np.array(MCS_tracks_ncfile.variables['length'])
#
#print('MCS_duration', MCS_duration)
#print('len(MCS_duration)', len(MCS_duration))
#
#MCS_duration_filtered = MCS_duration[mask]
#
#MCS_majoraxislength_init = np.array(MCS_tracks_ncfile.variables['majoraxislength'])[:,0]
#
#MCS_majoraxislength_init_2 = np.array(MCS_tracks_ncfile.variables['majoraxislength'])[:,1]
#
#MCS_majoraxislength_growth = MCS_majoraxislength_init_2 - MCS_majoraxislength_init
#
#MCS_majoraxislength_growth_filtered = MCS_majoraxislength_growth[mask]
#
#ccs_area_init = np.array(MCS_tracks_ncfile.variables['ccs_area'])[:,0]
#MCS_ccs_area_init = np.array(MCS_tracks_ncfile.variables['ccs_area'])[first_MCS_bol_array]
#
#print('ccs_area_init', ccs_area_init)
#print('MCS_ccs_area_init', MCS_ccs_area_init)
#
##print('first track area', np.array(MCS_tracks_ncfile.variables['ccs_area'])[0,:])
##print('first track status', np.array(MCS_tracks_ncfile.variables['mcs_status'])[0,:])
#
#MCS_ccs_area_growth = MCS_ccs_area_init - ccs_area_init
#
#MCS_ccs_area_growth_filtered = MCS_ccs_area_growth[mask]
#
#print('MCS_ccs_area_growth_filtered', MCS_ccs_area_growth_filtered)
#
#ccs_area_init_filtered = ccs_area_init[mask]
#
#print('ccs_area_init_filtered', ccs_area_init_filtered)
#
#MCS_growth_stage_time_length_filtered = MCS_growth_stage_time_length[mask]
#
#print('MCS_growth_stage_time_length_filtered', MCS_growth_stage_time_length_filtered)
#
#print('rate', MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered)
#print('median', np.nanmedian(MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered))
#
#condition_growth_time = 0.4 >= MCS_growth_stage_time_length_filtered
#
#print('dates', MCS_datetime_initiation_str_filtered[condition_growth_time])

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

########## Get area average environmental data centered on each MCS track chosen ######### 
#
#mean_bulk_shear_0_3km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#mean_bulk_shear_0_6km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#mean_bulk_shear_2_6km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#q_850 = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#MUCAPE = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#wv_flux_850_v = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#q_flux_850 = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#q_flux_850_v = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#prop_area_SALLJ = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#median_SALLJ_max_wind = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#median_SALLJ_height = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#
#SALLJ_count = 0
#low_cov_SALLJ_count = 0
#med_cov_SALLJ_count = 0
#high_cov_SALLJ_count = 0
#
## go through list of times/centroids for each MCS to get corresponding environmental conditions
#for count, (MCS_datetime, MCS_center_lon, MCS_center_lat) in enumerate(zip(MCS_datetime_initiation_str_filtered, MCS_center_lons_initiation_filtered, MCS_center_lats_initiation_filtered)):
#
#    # get the file necessary
#    path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'                  
#                                   
#    MCS_datetime_dt = datetime.strptime(MCS_datetime, '%Y-%m-%d_%H:%M:')
#    
#    if offset_MCS_and_conditions == True:
#        
#        conditions_datetime_dt = MCS_datetime_dt + timedelta(hours = hours_offset)
#        
#    else: #offset_MCS_and_conditions == False
#        
#        conditions_datetime_dt = MCS_datetime_dt
#    
#    MCS_datetime_str = conditions_datetime_dt.strftime('%Y%m%d')
#
#    wrf_file_name = '/wrfout_d01_%s-%s-%s_%s:00:00' %(conditions_datetime_dt.strftime('%Y'), conditions_datetime_dt.strftime('%m'), conditions_datetime_dt.strftime('%d'), conditions_datetime_dt.strftime('%H'))
#
#    wrf_file_path = path_wrf + MCS_datetime_str + wrf_file_name
#    
#    print('path of MCS within region: ', wrf_file_path)
#
#    # get the netCDF
#    wrf_ncfile = Dataset(wrf_file_path,'r')
#    
#    ####### get environmental profiles #######
#    
#    env_search_text = '__EnvArea_%s' %(env_search_area)
#
#    if env_search_area == '0.75fromMCScentroid':
#
#        # get lats/lons of region based on centroid
#        lat_bottom_left_env = MCS_center_lat - 0.75 
#        lon_bottom_left_env = MCS_center_lon - 0.75
#
#        lat_top_right_env = MCS_center_lat + 0.75 
#        lon_top_right_env = MCS_center_lon + 0.75
#            
#    elif env_search_area == '2.00fromMCScentroid':
#
#        # get lats/lons of region based on centroid
#        lat_bottom_left_env = MCS_center_lat - 2.00 
#        lon_bottom_left_env = MCS_center_lon - 2.00
#
#        lat_top_right_env = MCS_center_lat + 2.00 
#        lon_top_right_env = MCS_center_lon + 2.00
#        
#    else:
#        print('Please enter valid env search area') # will throw error
#    
#    #print(MCS_center_lat, MCS_center_lon)
#    #print(lat_bottom_left_env, lon_bottom_left_env)
#    #print(lat_top_right_env, lon_top_right_env)
#    
#    # get xy for regions
#    bottom_left_xy = ll_to_xy(wrf_ncfile, lat_bottom_left_env, lon_bottom_left_env)
#    top_right_xy = ll_to_xy(wrf_ncfile, lat_top_right_env, lon_top_right_env)  
#
#    # read in the variables 
#    hght_original = getvar(wrf_ncfile, 'height', msl=False, units='m')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in m NOTE: this is AGL not MSL!!
#    u_original = getvar(wrf_ncfile, 'ua', units='kt')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in kts
#    v_original = getvar(wrf_ncfile, 'va', units='kt')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in kts
#    
#    #print(u_original)
#    #print(hght_original)
#    
#    u2km = interplevel(u_original, hght_original, 2000)
#    v2km = interplevel(v_original, hght_original, 2000)
#
#    u3km = interplevel(u_original, hght_original, 3000)
#    v3km = interplevel(v_original, hght_original, 3000)
#    
#    u6km = interplevel(u_original, hght_original, 6000)
#    v6km = interplevel(v_original, hght_original, 6000)
#    
#    u_sfc = u_original[3, :, :]
#    v_sfc = v_original[3, :, :]
#
#    # calculate mean values
#    mean_u2km = np.nanmean(u2km)
#    mean_v2km = np.nanmean(v2km)
#    
#    mean_u3km = np.nanmean(u3km)
#    mean_v3km = np.nanmean(v3km)
#    
#    mean_u6km = np.nanmean(u6km)
#    mean_v6km = np.nanmean(v6km)
#    
#    mean_u_sfc = np.nanmean(u_sfc)
#    mean_v_sfc = np.nanmean(v_sfc)
#    
##    mean_u_sfc = u_mean_profile[3]
##    mean_v_sfc = v_mean_profile[3]
##    
##    mean_u2km = interplevel(u_mean_profile, hght_mean_profile, 2000)
##    mean_v2km = interplevel(v_mean_profile, hght_mean_profile, 2000)
##
##    mean_u3km = interplevel(u_mean_profile, hght_mean_profile, 3000)
##    mean_v3km = interplevel(v_mean_profile, hght_mean_profile, 3000)
##    
##    mean_u6km = interplevel(u_mean_profile, hght_mean_profile, 6000)
##    mean_v6km = interplevel(v_mean_profile, hght_mean_profile, 6000)
#    
#    # calcualte mean shear
#    mean_udiff_0_3 = mean_u3km - mean_u_sfc
#    mean_vdiff_0_3 = mean_v3km - mean_v_sfc
#    
#    mean_shear3km = np.sqrt(mean_udiff_0_3**2 + mean_vdiff_0_3**2)
#
#    mean_bulk_shear_0_3km[count] = mean_shear3km
#    
#    mean_udiff_0_6 = mean_u6km - mean_u_sfc
#    mean_vdiff_0_6 = mean_v6km - mean_v_sfc
#    
#    mean_shear6km = np.sqrt(mean_udiff_0_6**2 + mean_vdiff_0_6**2)
#
#    mean_bulk_shear_0_6km[count] = mean_shear6km
#    
#    mean_udiff_2_6 = mean_u6km - mean_u2km
#    mean_vdiff_2_6 = mean_v6km - mean_v2km
#    
#    mean_shear2_6km = np.sqrt(mean_udiff_2_6**2 + mean_vdiff_2_6**2)
#
#    mean_bulk_shear_2_6km[count] = mean_shear2_6km
#    
#    ######### calculate q at 850 hPa #########
#    level = 850
#    
#    # read in the variables for the subsetted area 
#    pres = getvar(wrf_ncfile, 'pressure')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])]
#
#    temp_c_subset = getvar(wrf_ncfile, 'temp', units='degC')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])]
#
#    RH_subset = getvar(wrf_ncfile, 'rh')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in kg/kg
#
#    # interpolate to the set level
#    temp_c_subset_level = np.squeeze(interplevel(temp_c_subset, pres, level))
#
#    RH_subset_level = np.squeeze(interplevel(RH_subset, pres, level))
#
#    e_sat_subset_level = SatVap(temp_c_subset_level)
#
#    e_subset_level = (RH_subset_level/100.) * e_sat_subset_level
#
#    w_subset_level = MixRatio(e_subset_level,level*100.)
#    #w = (e*Rs_da) / (Rs_v*(pres*100. - e))
#    q_subset_level = (w_subset_level / (w_subset_level+1))*1000 #g/kg
#
#    q_850[count] = np.nanmean(q_subset_level)
#    
#    ######### calculate MUCAPE #########
#    
#    MCAPE_original = getvar(wrf_ncfile, 'cape_2d')[0][int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])]
#    
#    MUCAPE[count] = np.nanmean(MCAPE_original)
#    
#    ######### calculate water vapor mass flux and q flux at 850 hPa #########
#    
#    v_850 = interplevel(v_original, pres, 850)
#    v_850_ms = v_850 * 0.5144 # convert from kts to m/s
#    
#    u_850 = interplevel(u_original, pres, 850)
#    u_850_ms = u_850 * 0.5144 # convert from kts to m/s
#    
#    
#    water_vapor_mass_flux_850_v = w_subset_level*v_850_ms # (g/kg)(m/s)
#    wv_flux_850_v[count] = np.nanmean(water_vapor_mass_flux_850_v)
#    
#    q_subset_level_kg_kg = q_subset_level/1000 # convert g/kg back to kg/kg
#    
#    qu_flux_850 = q_subset_level_kg_kg*u_850_ms
#    qv_flux_850 = q_subset_level_kg_kg*v_850_ms
#    q_total_flux_850 = np.sqrt((qu_flux_850**2)+(qv_flux_850**2))
#    
#    q_flux_850[count] = np.nanmean(q_total_flux_850)
#    q_flux_850_v[count] = np.nanmean(qv_flux_850)
#    
#    
#    ########## find if SALLJ is present and if so calculate characteristics #########
#    
#    # for finding LLJs #
#
#    # 2: high - allows for jets centered higher in the atmosphere, and a bit weaker compared to mean flow
#    crit = [2 ,3200, 5700, 19.4384, 7.77538] # relaxed to 10 m/s max and 4 m/s decrease
#
#    crit_num = crit[0]
#    max_search_hgt = crit[1]
#    min_search_hgt = crit[2]
#    max_wind_threshold = crit[3]
#    decrease_to_min_threshold = crit[4]
#    
#    SALLJ_search_text = '__SALLJarea_%s' %(SALLJ_search_area)
#
#    if SALLJ_search_area == '1deg3degBottomCentroid':
#
#        lat_bottom_left_SALLJ = MCS_center_lat
#        lon_bottom_left_SALLJ = MCS_center_lon - 1.50
#        lat_top_right_SALLJ = MCS_center_lat + 1.00 
#        lon_top_right_SALLJ = MCS_center_lon + 1.50
#
#    elif SALLJ_search_area == '2deg4degOffset1degNFromCentroid':
#
#        lat_bottom_left_SALLJ = MCS_center_lat + 1.00
#        lon_bottom_left_SALLJ = MCS_center_lon - 2.00
#        lat_top_right_SALLJ = MCS_center_lat + 3.00 
#        lon_top_right_SALLJ = MCS_center_lon + 2.00
#
#    elif SALLJ_search_area == '60-65W28-30SFixed':
#
#        lat_bottom_left_SALLJ = -30
#        lon_bottom_left_SALLJ = -65
#        lat_top_right_SALLJ = -28
#        lon_top_right_SALLJ = -60
#
#    else:
#        print('Please enter valid SALLJ search area') # will throw error
#    
#    # get xy for regions
#    bottom_left_xy = ll_to_xy(wrf_ncfile, lat_bottom_left_SALLJ, lon_bottom_left_SALLJ)
#    top_right_xy = ll_to_xy(wrf_ncfile, lat_top_right_SALLJ, lon_top_right_SALLJ) 
#
#    ############################# read in some variables #############################
#    #print('reading in variables')
#    # read in the variables subsetted to a specified area
#    pres_subset = getvar(wrf_ncfile, 'pressure')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])]
#    hght_subset = getvar(wrf_ncfile, 'height', msl=False, units='m')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in m NOTE: this is AGL not MSL!!
#    radar_refl_max = getvar(wrf_ncfile, 'REFD_MAX')[int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])]
#
#    speed, drct = getvar(wrf_ncfile, 'wspd_wdir', units='kt') # in kts
#
#    speed_subset = speed[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in kts
#
#    drct_subset = drct[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])]
#
##                # get wind speed and direction at chosen height 1 and calculate wind shear
##                wspd_chosen1 = interplevel(speed_subset, hght_subset, chosen_height1)
##                wdir_chosen1 = interplevel(drct_subset, hght_subset, chosen_height1)
##                
##                udiff = wspd_chosen1*np.cos(np.radians(270-wdir_chosen1)) - speed_subset[3, :, :]*np.cos(np.radians(270-drct_subset[3, :, :]))
##                vdiff = wspd_chosen1*np.sin(np.radians(270-wdir_chosen1)) - speed_subset[3, :, :]*np.sin(np.radians(270-drct_subset[3, :, :]))
##
##                shear_chosen1 = np.sqrt(udiff**2 + vdiff**2)
##                
##                print('wspd_chosen1', wspd_chosen1[30, 30].values)
##                print('wdir_chosen1', wdir_chosen1[30, 30].values)
##                print('wspd_sfc', speed_subset[3, 30, 30].values)
##                print('wdir_sfc', drct_subset[3, 30, 30].values)
##                print('shear_chosen1', shear_chosen1[30, 30].values)
##                
##                # get wind speed and direction at chosen height 2 and calculate wind shear
##                wspd_chosen2 = interplevel(speed_subset, hght_subset, chosen_height2)
##                wdir_chosen2 = interplevel(drct_subset, hght_subset, chosen_height2)
##                
##                udiff = wspd_chosen2*np.cos(np.radians(270-wdir_chosen2)) - speed_subset[3, :, :]*np.cos(np.radians(270-drct_subset[3, :, :]))
##                vdiff = wspd_chosen2*np.sin(np.radians(270-wdir_chosen2)) - speed_subset[3, :, :]*np.sin(np.radians(270-drct_subset[3, :, :]))
##
##                shear_chosen2 = np.sqrt(udiff**2 + vdiff**2)
##                
##                # get wind speed and direction at chosen height 3 and calculate wind shear
##                wspd_chosen3 = interplevel(speed_subset, hght_subset, chosen_height3)
##                wdir_chosen3 = interplevel(drct_subset, hght_subset, chosen_height3)
##                
##                udiff = wspd_chosen3*np.cos(np.radians(270-wdir_chosen3)) - speed_subset[3, :, :]*np.cos(np.radians(270-drct_subset[3, :, :]))
##                vdiff = wspd_chosen3*np.sin(np.radians(270-wdir_chosen3)) - speed_subset[3, :, :]*np.sin(np.radians(270-drct_subset[3, :, :]))
##
##                shear_chosen3 = np.sqrt(udiff**2 + vdiff**2)
##                
##                # get wind speed and direction at chosen height 6 and calculate wind shear
##                wspd_chosen6 = interplevel(speed_subset, hght_subset, chosen_height6)
##                wdir_chosen6 = interplevel(drct_subset, hght_subset, chosen_height6)
##                
##                udiff = wspd_chosen6*np.cos(np.radians(270-wdir_chosen6)) - speed_subset[3, :, :]*np.cos(np.radians(270-drct_subset[3, :, :]))
##                vdiff = wspd_chosen6*np.sin(np.radians(270-wdir_chosen6)) - speed_subset[3, :, :]*np.sin(np.radians(270-drct_subset[3, :, :]))
##
##                shear_chosen6 = np.sqrt(udiff**2 + vdiff**2)
##                
##                # get wind shear 2-6 km
##                
##                udiff = wspd_chosen6*np.cos(np.radians(270-wdir_chosen6)) - wspd_chosen2*np.cos(np.radians(270-wdir_chosen2))
##                vdiff = wspd_chosen6*np.sin(np.radians(270-wdir_chosen6)) - wspd_chosen2*np.sin(np.radians(270-wdir_chosen2))
##
##                shear_chosen2_6 = np.sqrt(udiff**2 + vdiff**2)
#
#    ############################# find SALLJ #############################
#
#    # interpolate variables to 250 meter vertical spacing up to the search heigth for the wind minimum 
#    interp_levels_min = np.arange(0,min_search_hgt+250,250)
#
#    pres_below_min = interplevel(pres_subset, hght_subset, interp_levels_min)
#    hght_below_min = interplevel(hght_subset, hght_subset, interp_levels_min)
#    speed_below_min = interplevel(speed_subset, hght_subset, interp_levels_min)
#    drct_below_min = interplevel(drct_subset, hght_subset, interp_levels_min)
#
#    # interpolate variables to 250 meter vertical spacing up to the search heigth for the wind maximum
#    interp_levels_max = np.arange(0,max_search_hgt+250,250)
#
#    pres_below_max = interplevel(pres_subset, hght_subset, interp_levels_max)
#    hght_below_max = interplevel(hght_subset, hght_subset, interp_levels_max)
#    speed_below_max = interplevel(speed_subset, hght_subset, interp_levels_max)
#    drct_below_max = interplevel(drct_subset, hght_subset, interp_levels_max)
#
#    #print(pres_below_max.shape)
#    #print(hght_below_max.shape)
#
#    #print('pres_below_max', pres_below_max.values)
#    #print('hght_below_max', hght_below_max.values)
#    #print('speed_below_max', speed_below_max.values)
#
#    ################ with xarray ###################
#
#    # get max wind
#    max_wind = speed_below_max.max(dim='level')
#
#    # get height of max wind
#    level_max_wind = hght_below_max.isel(level=speed_below_max.argmax('level'))
##                level_max_wind = speed_below_max.idxmax('level')
##                level_max_wind = speed_below_max[np.where(np.equal(speed_below_max,max_wind))].get_index('level')
#    #level_max_wind2 = speed_below_max.idxmax(dim='level')
#
#    # get pressure at max wind
#    pres_max_wind = pres_below_max.isel(level=speed_below_max.argmax('level'))
#
#    # get direction at max wind
#    drct_at_max_wind = drct_below_max.isel(level=speed_below_max.argmax('level'))
#
#    #print('drct_at_max_wind', drct_at_max_wind)
#
#    # set dim of level of max wind array to 'level' (height)
#    level_max_wind_subtract = xarray.DataArray(level_max_wind, dims=['south_north', 'west_east'])
#
#    #print('level_max_wind_subtract', level_max_wind_subtract)
#
#    #print('hght_below_max', hght_below_max.values)
#
#    # subtract the heights of max wind from the heights
#    hght_minus_level_max_wind = hght_below_min - level_max_wind_subtract
#
#    #print('hght_minus_level_max_wind', hght_minus_level_max_wind)
#
#    speed_below_min_masked_below_max = speed_below_min.where(hght_minus_level_max_wind > 0., np.nan)
#    hght_below_min_masked_below_max = hght_below_min.where(hght_minus_level_max_wind > 0., np.nan)
#
#    speed_below_min_masked_above_max = speed_below_min.where(hght_minus_level_max_wind < 0., np.nan)
#    hght_below_min_masked_above_max = hght_below_min.where(hght_minus_level_max_wind < 0., np.nan)
#
#    #print('speed_below_min_masked_below_max', speed_below_min_masked_below_max)
#
#    # get min wind above max
#    min_wind = speed_below_min_masked_below_max.min(dim='level')
#
#    min_wind_below_max = speed_below_min_masked_above_max.min(dim='level')
#
#    # get index of min above max
#
#    #min_wind_index = speed_below_min_masked_below_max.idxmin(dim='level')
#
#    #min_wind_index_below_max = speed_below_min_masked_above_max.idxmin(dim='level')
#
#    level_min_wind = hght_below_min_masked_below_max.isel(level=speed_below_min_masked_below_max.argmin('level'))
#
#    #level_min_wind_below_max = hght_below_min_masked_above_max.isel(level=speed_below_min_masked_above_max.argmin('level'))
#
#    #print('min_wind', min_wind.values)
#
#    # checks if max_wind meets threshold and keeps value if it does meet the threshold
#    max_wind_meeting_threshold = max_wind.where(max_wind > max_wind_threshold, np.nan)
#
#    #print('max_wind_meeting_threshold', max_wind_meeting_threshold.values)
#
#    #print('max_wind_meeting_threshold', max_wind_meeting_threshold)
#
#    # calculates decrease to min wind
#    decrease_to_min = max_wind_meeting_threshold - min_wind
#
#    decrease_to_min_below_max = max_wind_meeting_threshold - min_wind_below_max
#
#    #print('decrease_to_min', decrease_to_min)
#
#    # checks if decrease_to_min meets threshold and keeps value if it does meet the threshold
#    decrease_to_min_meeting_threshold = decrease_to_min.where(decrease_to_min > decrease_to_min_threshold, np.nan)
#
#    #print('decrease_to_min_meeting_threshold', decrease_to_min_meeting_threshold)
#
#    # checks to see if the values met the other criteria (was not replaced with nan) and if it does leave the value
#    drct_at_max_wind_meeting_threshold = drct_at_max_wind.where(np.isnan(decrease_to_min_meeting_threshold) == False, np.nan)
#
#    #print('drct_at_max_wind_meeting_threshold', drct_at_max_wind_meeting_threshold.values)
#
#    # checks to see if wind at max_wind is from a northerly direction and keeps value if it is
#    drct_at_max_wind_meeting_threshold = drct_at_max_wind_meeting_threshold.where((drct_at_max_wind_meeting_threshold <= 45) | (drct_at_max_wind_meeting_threshold >= 315), np.nan)
#
#    #print('drct_at_max_wind_meeting_threshold', drct_at_max_wind_meeting_threshold)
#    
#    num_points_SALLJ = np.count_nonzero(~np.isnan(drct_at_max_wind_meeting_threshold))
#    num_total_points = drct_at_max_wind_meeting_threshold.size
#    
#    #print('num points', num_total_points)
#    
#    proporation_SALLJ = num_points_SALLJ/num_total_points
#    
#    #print('proporation coverage: ', proporation_SALLJ)
#    
#    prop_area_SALLJ[count] = proporation_SALLJ
#    
#    if proporation_SALLJ >= 0.12:
#        
#        print('found SALLJ')
#        
#        SALLJ_count = SALLJ_count + 1
#        
#        ### get SALLJ characteristics for points where the SALLJ criteria is met
#
#        # get max_wind for points that meet SALLJ criteria
#        max_wind_SALLJs = max_wind.where(np.isnan(drct_at_max_wind_meeting_threshold) == False, np.nan)
#        
#        median_SALLJ_max_wind[count] = np.nanmedian(max_wind_SALLJs)
#        
#        # get level of max_wind for points that meet SALLJ criteria
#        level_max_wind_SALLJs = level_max_wind.where(np.isnan(drct_at_max_wind_meeting_threshold) == False, np.nan)
#        
#        median_SALLJ_height[count] = np.nanmedian(level_max_wind_SALLJs)
#        
#        #print('Median SALLJ wind: %s   Median SALLJ level: %s' %(median_SALLJ_max_wind, median_SALLJ_level))
#        
#        if proporation_SALLJ < 0.3:
#            
#            low_cov_SALLJ_count = low_cov_SALLJ_count + 1
#            
#            print('low_coverage_SALLJ')
#        elif proporation_SALLJ <= 0.3 and proporation_SALLJ < 0.6:
#            
#            print('med_coverage_SALLJ')
#            
#            med_cov_SALLJ_count = med_cov_SALLJ_count + 1
#            
#        elif proporation_SALLJ >= 0.6:
#            
#            print('high_coverage_SALLJ')
#            
#            high_cov_SALLJ_count = high_cov_SALLJ_count + 1
#        
#    else:
#        print('no SALLJ found')
#    
#    wrf_ncfile.close()
#    
#print('SALLJ_count', SALLJ_count)
#print('low_cov_SALLJ_count', low_cov_SALLJ_count)
#print('med_cov_SALLJ_count', med_cov_SALLJ_count)
#print('high_cov_SALLJ_count', high_cov_SALLJ_count)
#    
#general_outpath = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/firstStorms_InitialGrowthRate_MCS_characteristics_vs_environment/'
#
##specific_outpath = '%sarea_%s%s%s/data/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)
#
### output arrays as pickles (.dat files)
##pickle.dump(MCS_duration_filtered, open(general_outpath + specific_outpath + "2%s_duration2%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
###pickle.dump(MCS_ccs_area_filtered, open(general_outpath + specific_outpath + "%s_ccs_area%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(MCS_majoraxislength_growth_filtered, open(general_outpath + specific_outpath + "2%s_majoraxislength_growth%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(MCS_ccs_area_growth_filtered, open(general_outpath + specific_outpath + "2%s_ccs_area_growth%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
#pickle.dump(ccs_area_init_filtered, open(general_outpath + "2%s_ccs_area_init%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(MCS_growth_stage_time_length_filtered, open(general_outpath + specific_outpath + "2%s_growth_stage_time_length%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(mean_bulk_shear_0_3km, open(general_outpath + specific_outpath + "2%s_mean_bulk_shear_0_3km%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(mean_bulk_shear_0_6km, open(general_outpath + specific_outpath + "2%s_mean_bulk_shear_0_6km%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(mean_bulk_shear_2_6km, open(general_outpath + specific_outpath + "2%s_mean_bulk_shear_2_6km%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(q_850, open(general_outpath + specific_outpath + "2%s_q_850%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(MUCAPE, open(general_outpath + specific_outpath + "2%s_MUCAPE%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(wv_flux_850_v, open(general_outpath + specific_outpath + "2%s_wv_flux_850_v%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(q_flux_850, open(general_outpath + specific_outpath + "2%s_q_flux_850%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(q_flux_850_v, open(general_outpath + specific_outpath + "2%s_q_flux_850_v%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(prop_area_SALLJ, open(general_outpath + specific_outpath + "2%s_prop_area_SALLJ%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(median_SALLJ_max_wind, open(general_outpath + specific_outpath + "2%s_median_SALLJ_max_wind%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
##pickle.dump(median_SALLJ_height, open(general_outpath + specific_outpath + "2%s_median_SALLJ_height%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "wb"))
### Note: this could be improved by using a dict object
##    