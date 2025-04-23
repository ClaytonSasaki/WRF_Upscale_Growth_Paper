#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on fri Sept 30 11:59:01 2022

@author: crs326

Compares SALLJ height (and speed) with 0-3km wind shear

"""

# Use non-interactive backend for environments without display (e.g., remote servers)
import matplotlib
matplotlib.use('Agg') 

import numpy as np
import pickle

# Set plot-wide font size
matplotlib.rcParams.update({'font.size': 18})

# ----------------- USER-DEFINED SETTINGS ------------------

filter_label = '_filtered_init_loc_and_start_type'
MCS_type = 'all_MCS'               # Options: 'all_MCS', 'robustMCS'
MCS_init_area = 'large_area1'      # Choose region for MCS initiation
SALLJ_search_area = '2deg4degOffset1degNFromCentroid' # '2deg4degOffset1degNFromCentroid', '1deg3degBottomCentroid', '60-65W28-30SFixed'
env_search_area = '2.00fromMCScentroid' # '0.75fromMCScentroid', '2.00fromMCScentroid'

add_SALLJ = False
add_statistics = False

events_removed = True

# NOTE: some variables may not be available script that created them was before they were added. Can either rerun that script or comment out certain varibale 
# ----------------------------------------------------------

# get corresponding file labels using chosen inputs
SALLJ_search_text = '__SALLJarea_%s' %(SALLJ_search_area)
env_search_text = '__EnvArea_%s' %(env_search_area)

# Set file label for type of MCS
if MCS_type == 'all_MCS':
    MCS_file_label = 'MCS'
elif MCS_type == 'robustMCS':
    MCS_file_label = 'robustMCS'
else:
    print('MUST choose valid MCS_type')

# get input files paths
general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/firstStorms_InitialGrowthRate_MCS_characteristics_vs_environment_all_points/'

specific_inpath = '%sarea_%s%s%s/data/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

# read in files
MCS_prop_area_SALLJ_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_prop_area_SALLJ_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
##MCS_duration_filtered = pickle.load(open("MCS_duration%s.dat" %(filter_label), "rb")) # NOTE: this duration from 'mcs_length' variable is all zeros
#MCS_duration_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_duration2_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_majoraxislength_growth_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_majoraxislength_growth_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_growth_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_growth_stage_time_length_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_growth_stage_time_length_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_ccs_area_growth_filtered_all_points_byEvent_2hr = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_2hr_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
bulk_shear_0_3km_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_bulk_shear_0_3km_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
bulk_shear_0_6km_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_bulk_shear_0_6km_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
bulk_shear_2_6km_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_bulk_shear_2_6km_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_max_wind_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_median_SALLJ_max_wind_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_height_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_median_SALLJ_height_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_q_850_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_q_850_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_MUCAPE_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_MUCAPE_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_q_flux_850_v_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_q_flux_850_v_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_refl_coverage_grt0_byEvent = pickle.load( open(general_path + specific_inpath + "%s_refl_coverage_grt0_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_refl_coverage_grt30_byEvent = pickle.load( open(general_path + specific_inpath + "%s_refl_coverage_grt30_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))

# Filter events
if events_removed == True:

    events_bool = np.full(len(bulk_shear_0_3km_all_points_byEvent), True)
    remove_events_nums = [1,7,10,154,72,73,81,84,88,93,95,99,136,141,142,153,57]
    events_bool[remove_events_nums] = False
    events_removed_label = '_events_removed'
    
else:
    
    events_bool = np.full(len(bulk_shear_0_3km_all_points_byEvent), True)
    events_removed_label = ''  

MCS_prop_area_SALLJ_all_points_byEvent = np.array(MCS_prop_area_SALLJ_all_points_byEvent, dtype="object")[events_bool]
median_SALLJ_max_wind_all_points_byEvent = np.array(median_SALLJ_max_wind_all_points_byEvent, dtype="object")[events_bool]
median_SALLJ_height_all_points_byEvent = np.array(median_SALLJ_height_all_points_byEvent, dtype="object")[events_bool]
MCS_ccs_area_growth_filtered_all_points_byEvent = np.array(MCS_ccs_area_growth_filtered_all_points_byEvent, dtype="object")[events_bool]
MCS_growth_stage_time_length_filtered_all_points_byEvent = np.array(MCS_growth_stage_time_length_filtered_all_points_byEvent, dtype="object")[events_bool]
bulk_shear_0_3km_all_points_byEvent = np.array(bulk_shear_0_3km_all_points_byEvent, dtype="object")[events_bool]
bulk_shear_0_6km_all_points_byEvent = np.array(bulk_shear_0_6km_all_points_byEvent, dtype="object")[events_bool]
bulk_shear_2_6km_all_points_byEvent = np.array(bulk_shear_2_6km_all_points_byEvent, dtype="object")[events_bool]
MCS_q_850_all_points_byEvent = np.array(MCS_q_850_all_points_byEvent, dtype="object")[events_bool]
MCS_MUCAPE_all_points_byEvent = np.array(MCS_MUCAPE_all_points_byEvent, dtype="object")[events_bool]
MCS_q_flux_850_v_all_points_byEvent = np.array(MCS_q_flux_850_v_all_points_byEvent, dtype="object")[events_bool]

# calculate mean/median in variables and save them into 1D array (for MCS and SALLJ characteristics, first value is taken as does not vary spatially)

MCS_prop_area_SALLJ = np.array([arr[0,0] for arr in MCS_prop_area_SALLJ_all_points_byEvent])
median_SALLJ_max_wind = np.array([arr[0,0] for arr in median_SALLJ_max_wind_all_points_byEvent])
median_SALLJ_height = np.array([arr[0,0] for arr in median_SALLJ_height_all_points_byEvent])
MCS_ccs_area_growth_filtered = np.array([arr[0,0] for arr in MCS_ccs_area_growth_filtered_all_points_byEvent])
MCS_growth_stage_time_length_filtered = np.array([arr[0,0] for arr in MCS_growth_stage_time_length_filtered_all_points_byEvent])
mean_bulk_shear_0_3km = np.array([np.nanmean(arr) for arr in bulk_shear_0_3km_all_points_byEvent])
mean_bulk_shear_0_6km = np.array([np.nanmean(arr) for arr in bulk_shear_0_6km_all_points_byEvent])
mean_bulk_shear_2_6km = np.array([np.nanmean(arr) for arr in bulk_shear_2_6km_all_points_byEvent])
mean_q_850 = np.array([np.nanmean(arr) for arr in MCS_q_850_all_points_byEvent])
mean_MUCAPE = np.array([np.nanmean(arr) for arr in MCS_MUCAPE_all_points_byEvent])
mean_q_flux_850_v = np.array([np.nanmean(arr) for arr in MCS_q_flux_850_v_all_points_byEvent])
median_bulk_shear_0_3km = np.array([np.nanmedian(arr) for arr in bulk_shear_0_3km_all_points_byEvent])
median_bulk_shear_0_6km = np.array([np.nanmedian(arr) for arr in bulk_shear_0_6km_all_points_byEvent])
median_bulk_shear_2_6km = np.array([np.nanmedian(arr) for arr in bulk_shear_2_6km_all_points_byEvent])
median_q_850 = np.array([np.nanmedian(arr) for arr in MCS_q_850_all_points_byEvent])
median_MUCAPE= np.array([np.nanmedian(arr) for arr in MCS_MUCAPE_all_points_byEvent])
q75_bulk_shear_0_3km = np.array([np.quantile(arr,0.75) for arr in bulk_shear_0_3km_all_points_byEvent])
q75_bulk_shear_0_6km = np.array([np.quantile(arr,0.75) for arr in bulk_shear_0_6km_all_points_byEvent])
q75_bulk_shear_2_6km = np.array([np.quantile(arr,0.75) for arr in bulk_shear_2_6km_all_points_byEvent])
q75_q_850 = np.array([np.quantile(arr,0.75) for arr in MCS_q_850_all_points_byEvent])
q75_MUCAPE= np.array([np.quantile(arr,0.75) for arr in MCS_MUCAPE_all_points_byEvent])

# calculate growth rate

MCS_ccs_area_growth_rate_filtered = MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered
print(MCS_ccs_area_growth_rate_filtered)
#print(mean_q_850)

# remove events where growth rate is nan (because MCS was found at time 0)

nan_indices = np.isnan(MCS_ccs_area_growth_rate_filtered)

MCS_prop_area_SALLJ = np.delete(MCS_prop_area_SALLJ, nan_indices)
median_SALLJ_max_wind = np.delete(median_SALLJ_max_wind, nan_indices)
median_SALLJ_height = np.delete(median_SALLJ_height, nan_indices)
MCS_ccs_area_growth_filtered = np.delete(MCS_ccs_area_growth_filtered, nan_indices)
MCS_growth_stage_time_length_filtered = np.delete(MCS_growth_stage_time_length_filtered, nan_indices)
MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, nan_indices)
mean_bulk_shear_0_3km = np.delete(mean_bulk_shear_0_3km, nan_indices)
mean_bulk_shear_0_6km = np.delete(mean_bulk_shear_0_6km, nan_indices)
mean_bulk_shear_2_6km = np.delete(mean_bulk_shear_2_6km, nan_indices)
mean_q_850 = np.delete(mean_q_850, nan_indices)
mean_MUCAPE = np.delete(mean_MUCAPE, nan_indices)
mean_q_flux_850_v = np.delete(mean_q_flux_850_v, nan_indices)
median_bulk_shear_0_3km = np.delete(median_bulk_shear_0_3km, nan_indices)
median_bulk_shear_0_6km = np.delete(median_bulk_shear_0_6km, nan_indices)
median_bulk_shear_2_6km = np.delete(median_bulk_shear_2_6km, nan_indices)
median_q_850 = np.delete(median_q_850, nan_indices)
median_MUCAPE = np.delete(median_MUCAPE, nan_indices)
q75_bulk_shear_0_3km = np.delete(q75_bulk_shear_0_3km, nan_indices)
q75_bulk_shear_0_6km = np.delete(q75_bulk_shear_0_6km, nan_indices)
q75_bulk_shear_2_6km = np.delete(q75_bulk_shear_2_6km, nan_indices)
q75_q_850 = np.delete(q75_q_850, nan_indices)
q75_MUCAPE = np.delete(q75_MUCAPE, nan_indices)

############ Change variable for seperation ############

var_seperation = MCS_ccs_area_growth_rate_filtered # MCS_ccs_area_growth_filtered, MCS_growth_stage_time_length_filtered
var_seperation_name = 'growth_rate' # 'growth_rate', 'growth_stage_time_length'
seperation_threshold = np.nanmedian(var_seperation)

# create masked array based on MCS css area growth rate

rapid_growth = var_seperation >= seperation_threshold

MCS_prop_area_SALLJ_rapid_growth = MCS_prop_area_SALLJ[rapid_growth]
median_SALLJ_max_wind_rapid_growth = median_SALLJ_max_wind[rapid_growth]
median_SALLJ_height_rapid_growth = median_SALLJ_height[rapid_growth]
#MCS_ccs_area_growth_rate_filtered_rapid_growth = MCS_ccs_area_growth_rate_filtered[rapid_growth]
mean_q_850_rapid_growth = mean_q_850[rapid_growth]
mean_MUCAPE_rapid_growth = mean_MUCAPE[rapid_growth]
mean_q_flux_850_v_rapid_growth = mean_q_flux_850_v[rapid_growth]
mean_bulk_shear_0_3km_rapid_growth = mean_bulk_shear_0_3km[rapid_growth]
mean_bulk_shear_0_6km_rapid_growth = mean_bulk_shear_0_6km[rapid_growth]
mean_bulk_shear_2_6km_rapid_growth = mean_bulk_shear_2_6km[rapid_growth]

slow_growth = var_seperation < seperation_threshold

MCS_prop_area_SALLJ_slow_growth = MCS_prop_area_SALLJ[slow_growth]
median_SALLJ_max_wind_slow_growth = median_SALLJ_max_wind[slow_growth]
median_SALLJ_height_slow_growth = median_SALLJ_height[slow_growth]
#MCS_ccs_area_growth_rate_filtered_slow_growth = MCS_ccs_area_growth_rate_filtered[slow_growth]
mean_q_850_slow_growth = mean_q_850[slow_growth]
mean_MUCAPE_slow_growth = mean_MUCAPE[slow_growth]
mean_q_flux_850_v_slow_growth = mean_q_flux_850_v[slow_growth]
mean_bulk_shear_0_3km_slow_growth = mean_bulk_shear_0_3km[slow_growth]
mean_bulk_shear_0_6km_slow_growth = mean_bulk_shear_0_6km[slow_growth]
mean_bulk_shear_2_6km_slow_growth = mean_bulk_shear_2_6km[slow_growth]

## calculate numbers in each growth rate group and statistics (medians, quartiles)
#
#num_MCS_all_growth = len(mean_q_850)
#
#all_growth_q_850_q75 = np.quantile(mean_q_850,0.75)
#all_growth_q_850_q50 = np.quantile(mean_q_850,0.50)
#all_growth_q_850_q25 = np.quantile(mean_q_850,0.25)
#
#num_MCS_growth_grt_10000 = len(mean_q_850_rapid_growth)
#
#growth_grt_10000_q_850_q75 = np.quantile(mean_q_850_rapid_growth,0.75)
#growth_grt_10000_q_850_q50 = np.quantile(mean_q_850_rapid_growth,0.50)
#growth_grt_10000_q_850_q25 = np.quantile(mean_q_850_rapid_growth,0.25)
#
#num_MCS_growth_less_10000 = len(mean_q_850_slow_growth)
#
#growth_less_10000_q_850_q75 = np.quantile(mean_q_850_slow_growth,0.75)
#growth_less_10000_q_850_q50 = np.quantile(mean_q_850_slow_growth,0.50)
#growth_less_10000_q_850_q25 = np.quantile(mean_q_850_slow_growth,0.25)
#
#print('MCS of ALL growth rate, n = %d' %(num_MCS_all_growth))
#print('850hPa q: 25pct - %.1f, 50pct - %.1f, 75pct - %.1f' %(all_growth_q_850_q25, all_growth_q_850_q50, all_growth_q_850_q75))
#print(sorted(mean_q_850))
#
#print('MCS growth rate >= %.0f km^2/hr, n = %d' %(seperation_threshold, num_MCS_growth_grt_10000))
#print('850hPa q: 25pct - %.1f, 50pct - %.1f, 75pct - %.1f' %(growth_grt_10000_q_850_q25, growth_grt_10000_q_850_q50, growth_grt_10000_q_850_q75))
#print(sorted(mean_q_850_rapid_growth))
#
#print('MCS growth rate < %.0f km^2/hr, n = %d' %(seperation_threshold, num_MCS_growth_less_10000))
#print('850hPa q: 25pct - %.1f, 50pct - %.1f, 75pct - %.1f' %(growth_less_10000_q_850_q25, growth_less_10000_q_850_q50, growth_less_10000_q_850_q75))
#print(sorted(mean_q_850_slow_growth))

MCS_prop_area_SALLJ_rapid_growth = MCS_prop_area_SALLJ_rapid_growth[MCS_prop_area_SALLJ_rapid_growth >= 0.30]
MCS_prop_area_SALLJ_slow_growth = MCS_prop_area_SALLJ_slow_growth[MCS_prop_area_SALLJ_slow_growth >= 0.30]

print('MCS_prop_area_SALLJ_rapid_growth', np.sort(MCS_prop_area_SALLJ_rapid_growth))
print('median_SALLJ_height_rapid_growth', np.sort(median_SALLJ_height_rapid_growth))
print('median_SALLJ_max_wind_rapid_growth', np.sort(median_SALLJ_max_wind_rapid_growth))
print('count rapid SALLJ: %s / %s' %(np.count_nonzero(~np.isnan(median_SALLJ_max_wind_rapid_growth)), len(median_SALLJ_max_wind_rapid_growth)))
print('Median SALLJ max wind: %s' %(np.nanmedian(median_SALLJ_max_wind_rapid_growth)))
print('Median SALLJ height: %s' %(np.nanmedian(median_SALLJ_height_rapid_growth)))
print('Median prop SALLJ: %s' %(np.nanmedian(MCS_prop_area_SALLJ_rapid_growth)))

print('MCS_prop_area_SALLJ_slow_growth', np.sort(MCS_prop_area_SALLJ_slow_growth))
print('median_SALLJ_height_slow_growth', np.sort(median_SALLJ_height_slow_growth))
print('median_SALLJ_max_wind_slow_growth', np.sort(median_SALLJ_max_wind_slow_growth))
print('count slow SALLJ: %s / %s' %(np.count_nonzero(~np.isnan(median_SALLJ_max_wind_slow_growth)), len(median_SALLJ_max_wind_slow_growth)))
print('Median SALLJ max wind: %s' %(np.nanmedian(median_SALLJ_max_wind_slow_growth)))
print('Median SALLJ height: %s' %(np.nanmedian(median_SALLJ_height_slow_growth)))
print('Median prop SALLJ: %s' %(np.nanmedian(MCS_prop_area_SALLJ_slow_growth)))


    