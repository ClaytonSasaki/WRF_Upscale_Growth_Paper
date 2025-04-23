#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Oct 21 2024

@author: crs326

Creates histograms of environmental variables seperated by MCS characteristics derived all points data. Can calculate statistic (e.g. spread) over the area at first storms times.

Specifically used for plotting the proportion of the area covered by reflectivity.

"""

# Use non-interactive backend for environments without display (e.g., remote servers)
import matplotlib
matplotlib.use('Agg') 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
from scipy import stats

matplotlib.rcParams.update({'font.size': 15})

#################### VARIABLES TO CHANGE ##########################

filter_label = '_filtered_init_loc_and_start_type'

MCS_type = 'all_MCS' # 'all_MCS' or 'robustMCS'
MCS_init_area = 'large_area1'
SALLJ_search_area = '2deg4degOffset1degNFromCentroid' # '2deg4degOffset1degNFromCentroid', '1deg3degBottomCentroid', '60-65W28-30SFixed'
env_search_area = '2.00fromMCScentroid' # '0.75fromMCScentroid', '2.00fromMCScentroid'

plot_rapid_growth_dist = False
plot_slow_growth_dist = False
plot_all_growth_dist = True

##############################################################################

# get corresponding file labels using chosen inputs
SALLJ_search_text = '__SALLJarea_%s' %(SALLJ_search_area)
env_search_text = '__EnvArea_%s' %(env_search_area)

if MCS_type == 'all_MCS':
    MCS_file_label = 'MCS'
elif MCS_type == 'robustMCS':
    MCS_file_label = 'robustMCS'
else:
    print('MUST choose valid MCS_type')

# Input paths
general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/MCS_characteristics_vs_environment_all_points/'

specific_inpath = '%sarea_%s%s%s/data/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

# Load first storm data
MCS_prop_area_SALLJ_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_prop_area_SALLJ_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_duration_filtered = pickle.load(open("MCS_duration%s.dat" %(filter_label), "rb")) # NOTE: this duration from 'mcs_length' variable is all zeros
MCS_duration_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_duration2_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_majoraxislength_growth_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_majoraxislength_growth_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_growth_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_growth_filtered_all_points_byEvent_2hr = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_2hr_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
bulk_shear_0_3km_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_bulk_shear_0_3km_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
bulk_shear_0_6km_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_bulk_shear_0_6km_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
bulk_shear_2_6km_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_bulk_shear_2_6km_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_max_wind_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_median_SALLJ_max_wind_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_height_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_median_SALLJ_height_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_q_850_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_q_850_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_MUCAPE_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_MUCAPE_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_wv_flux_850_v_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_wv_flux_850_v_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_refl_coverage_grt0_byEvent = pickle.load( open(general_path + specific_inpath + "%s_refl_coverage_grt0_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_refl_coverage_grt30_byEvent = pickle.load( open(general_path + specific_inpath + "%s_refl_coverage_grt30_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))

# calculate spread in variables and save them into 1D array (for MCS and SALLJ characteristics, first value is taken as does not vary spatially)

MCS_prop_area_SALLJ_spreadByEvent = np.array([arr[0,0] for arr in MCS_prop_area_SALLJ_all_points_byEvent])
MCS_duration_filtered_spreadByEvent = np.array([arr[0,0] for arr in MCS_duration_filtered_all_points_byEvent])
MCS_majoraxislength_growth_filtered_spreadByEvent = np.array([arr[0,0] for arr in MCS_majoraxislength_growth_filtered_all_points_byEvent])
MCS_ccs_area_growth_filtered_spreadByEvent = np.array([arr[0,0] for arr in MCS_ccs_area_growth_filtered_all_points_byEvent])
MCS_ccs_area_growth_filtered_spreadByEvent_2hr = np.array([arr[0,0] for arr in MCS_ccs_area_growth_filtered_all_points_byEvent_2hr])
print('MCS_ccs_area_growth_filtered_spreadByEvent', MCS_ccs_area_growth_filtered_spreadByEvent)
bulk_shear_0_3km_spreadByEvent = np.array([np.nanmax(arr)-np.nanmin(arr) for arr in bulk_shear_0_3km_all_points_byEvent])
bulk_shear_0_6km_spreadByEvent = np.array([np.nanmax(arr)-np.nanmin(arr) for arr in bulk_shear_0_6km_all_points_byEvent])
bulk_shear_2_6km_spreadByEvent = np.array([np.nanmax(arr)-np.nanmin(arr) for arr in bulk_shear_2_6km_all_points_byEvent])
q_850_spreadByEvent = np.array([np.nanmax(arr)-np.nanmin(arr) for arr in MCS_q_850_all_points_byEvent])
median_SALLJ_max_wind_spreadByEvent =  np.array([arr[0,0] for arr in median_SALLJ_max_wind_all_points_byEvent])
median_SALLJ_height_spreadByEvent =  np.array([arr[0,0] for arr in median_SALLJ_height_all_points_byEvent])

q_850_quart75ByEvent = np.array([np.percentile(arr, 75) for arr in MCS_q_850_all_points_byEvent])
q_850_quart90ByEvent = np.array([np.percentile(arr, 90) for arr in MCS_q_850_all_points_byEvent])
MCS_MUCAPE_quart90ByEvent = np.array([np.percentile(arr, 90) for arr in MCS_MUCAPE_all_points_byEvent])
MCS_wv_flux_quart90ByEvent = np.array([np.percentile(arr, 90) for arr in MCS_wv_flux_850_v_all_points_byEvent])

q_850_mean = np.array([np.nanmean(arr) for arr in MCS_q_850_all_points_byEvent])
MCS_MUCAPE_mean = np.array([np.nanmean(arr) for arr in MCS_MUCAPE_all_points_byEvent])
MCS_wv_flux_850_mean = np.array([np.nanmean(arr) for arr in MCS_wv_flux_850_v_all_points_byEvent])

median_MCS_ccs_area_growth_filtered_spreadByEvent = np.nanmedian(MCS_ccs_area_growth_filtered_spreadByEvent)
print(median_MCS_ccs_area_growth_filtered_spreadByEvent)
median_MCS_ccs_area_growth_filtered_spreadByEvent_2hr = np.nanmedian(MCS_ccs_area_growth_filtered_spreadByEvent_2hr)
print(median_MCS_ccs_area_growth_filtered_spreadByEvent_2hr)

############ Change variable for seperation ############

var_seperation = MCS_ccs_area_growth_filtered_spreadByEvent_2hr
var_seperation_name = 'ccs_area_growth_2hr'
seperation_threshold = median_MCS_ccs_area_growth_filtered_spreadByEvent_2hr

########################################################

mask1 = var_seperation >= seperation_threshold

MCS_prop_area_SALLJ_spreadByEvent_mask1 = MCS_prop_area_SALLJ_spreadByEvent[mask1]
MCS_duration_filtered_spreadByEvent_mask1 = MCS_duration_filtered_spreadByEvent[mask1]
MCS_majoraxislength_growth_filtered_spreadByEvent_mask1 = MCS_majoraxislength_growth_filtered_spreadByEvent[mask1]
MCS_ccs_area_growth_filtered_spreadByEvent_mask1 = MCS_ccs_area_growth_filtered_spreadByEvent[mask1]
MCS_ccs_area_growth_filtered_2hr_spreadByEvent_mask1 = MCS_ccs_area_growth_filtered_spreadByEvent_2hr[mask1]
bulk_shear_0_3km_spreadByEvent_mask1 = bulk_shear_0_3km_spreadByEvent[mask1]
bulk_shear_0_6km_spreadByEvent_mask1 = bulk_shear_0_6km_spreadByEvent[mask1]
bulk_shear_2_6km_spreadByEvent_mask1 = bulk_shear_2_6km_spreadByEvent[mask1]
median_SALLJ_max_wind_spreadByEvent_mask1 = median_SALLJ_max_wind_spreadByEvent[mask1]
median_SALLJ_height_spreadByEvent_mask1 = median_SALLJ_height_spreadByEvent[mask1]
q_850_spreadByEvent_mask1 = q_850_spreadByEvent[mask1]
q_850_quart75ByEvent_mask1 = q_850_quart75ByEvent[mask1]
q_850_quart90ByEvent_mask1 = q_850_quart90ByEvent[mask1]
MCS_MUCAPE_quart90ByEvent_mask1 = MCS_MUCAPE_quart90ByEvent[mask1]
MCS_wv_flux_quart90ByEvent_mask1 = MCS_wv_flux_quart90ByEvent[mask1]
q_850_mean_mask1 = q_850_mean[mask1]
MCS_MUCAPE_mean_mask1 = MCS_MUCAPE_mean[mask1]
MCS_wv_flux_850_mean_mask1 = MCS_wv_flux_850_mean[mask1]
#MCS_refl_coverage_grt0_byEvent_mask1 = MCS_refl_coverage_grt0_byEvent[mask1]

print('# rapid growth MCS', len(MCS_ccs_area_growth_filtered_2hr_spreadByEvent_mask1))

mask2 = var_seperation < seperation_threshold

MCS_prop_area_SALLJ_spreadByEvent_mask2 = MCS_prop_area_SALLJ_spreadByEvent[mask2]
MCS_duration_filtered_spreadByEvent_mask2 = MCS_duration_filtered_spreadByEvent[mask2]
MCS_majoraxislength_growth_filtered_spreadByEvent_mask2 = MCS_majoraxislength_growth_filtered_spreadByEvent[mask2]
MCS_ccs_area_growth_filtered_spreadByEvent_mask2 = MCS_ccs_area_growth_filtered_spreadByEvent[mask2]
MCS_ccs_area_growth_filtered_2hr_spreadByEvent_mask2 = MCS_ccs_area_growth_filtered_spreadByEvent_2hr[mask2]
bulk_shear_0_3km_spreadByEvent_mask2 = bulk_shear_0_3km_spreadByEvent[mask2]
bulk_shear_0_6km_spreadByEvent_mask2 = bulk_shear_0_6km_spreadByEvent[mask2]
bulk_shear_2_6km_spreadByEvent_mask2 = bulk_shear_2_6km_spreadByEvent[mask2]
median_SALLJ_max_wind_spreadByEvent_mask2 = median_SALLJ_max_wind_spreadByEvent[mask2]
median_SALLJ_height_spreadByEvent_mask2 = median_SALLJ_height_spreadByEvent[mask2]
q_850_spreadByEvent_mask2 = q_850_spreadByEvent[mask2]
q_850_quart75ByEvent_mask2 = q_850_quart75ByEvent[mask2]
q_850_quart90ByEvent_mask2 = q_850_quart90ByEvent[mask2]
MCS_MUCAPE_quart90ByEvent_mask2 = MCS_MUCAPE_quart90ByEvent[mask2]
MCS_wv_flux_quart90ByEvent_mask2 = MCS_wv_flux_quart90ByEvent[mask2]
q_850_mean_mask2 = q_850_mean[mask2]
MCS_MUCAPE_mean_mask2 = MCS_MUCAPE_mean[mask2]
MCS_wv_flux_850_mean_mask2 = MCS_wv_flux_850_mean[mask2]
#MCS_refl_coverage_grt0_byEvent_mask2 = MCS_refl_coverage_grt0_byEvent[mask2]


print('# slow growth MCS', len(MCS_ccs_area_growth_filtered_2hr_spreadByEvent_mask2))

############ Change variable plotting ############ 

data_all_growth = MCS_refl_coverage_grt30_byEvent
#data_rapid_growth = MCS_refl_coverage_grt0_byEvent_mask1
#data_slow_growth = MCS_refl_coverage_grt0_byEvent_mask2
variable_name = 'MCS_refl_coverage_grt30_byEvent' # q_850_spreadByEvent, bulk_shear_0_3km_spreadByEvent, MCS_MUCAPE_quart90ByEvent, MCS_wv_flux_quart90ByEvent, q_850_mean, MCS_MUCAPE_mean
x_label = 'refl coverage (dBZ > 30) by event' # 850hPa q spread by event

################################################## 

fig, ax = plt.subplots()

if plot_all_growth_dist == True:
    
    median_value = np.nanmedian(data_all_growth)
    dist_text_all_growth = '_all'   
    plt.hist(data_all_growth, bins=30, alpha=0.3, label='all MCS', color='gray')
    plt.axvline(median_value, color='gray', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.2f}')

else:
    dist_text_all_growth = ''
    
if plot_rapid_growth_dist == True:
    
    median_value = np.nanmedian(data_rapid_growth)
    dist_text_rapid_growth = '_rapid'   
    plt.hist(data_rapid_growth, bins=30, alpha=0.3, label='rapid growth MCS', color='blue')
    plt.axvline(median_value, color='blue', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.2f}')
    
else:
    dist_text_rapid_growth = ''
    
if plot_slow_growth_dist == True:
    
    median_value = np.nanmedian(data_slow_growth)
    dist_text_slow_growth = '_slow'
    plt.hist(data_slow_growth, bins=30, alpha=0.3, label='slow growth MCS', color='red')
    plt.axvline(median_value, color='red', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.2f}')
    
else:
    dist_text_slow_growth = ''

# Adding labels and title
plt.xlabel(x_label)
plt.ylabel('Frequency')

# Show legend
plt.legend()

plt.tight_layout()

print('saving')

specific_outpath = '%sarea_%s%s%s/plots/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

plt.savefig(general_path + specific_outpath + '%s_hist_by_%s_%d%s%s%s%s.png' %(variable_name, var_seperation_name, seperation_threshold, dist_text_all_growth, dist_text_rapid_growth, dist_text_slow_growth, filter_label), dpi=200)

print('saved')
    