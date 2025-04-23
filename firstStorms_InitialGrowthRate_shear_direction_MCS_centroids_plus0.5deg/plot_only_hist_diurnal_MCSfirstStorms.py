#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Oct 24 2024

@author: crs326

Plots hist of the diurnal timing of MCS firstStorms Additionally seperates by MCS growth rate.

"""

# Use non-interactive backend for environments without display (e.g., remote servers)
import matplotlib
matplotlib.use('Agg') 

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pickle

# Set plot-wide font size
matplotlib.rcParams.update({'font.size': 18})

# ----------------- USER-DEFINED SETTINGS ------------------

filter_label = '_filtered_init_loc_and_start_type'
MCS_type = 'all_MCS'               # Options: 'all_MCS', 'robustMCS'
MCS_init_area = 'large_area1'      # Choose region for MCS initiation
env_offset_MCS = False             # If True, offset environment data from MCS start time
hours_offset = -1                  # Hours offset for environmental sampling
seperation_plot = False            # Whether to plot MCS growth rate separation
growth_type = 'slow'          
events_removed = True              # Whether to exclude certain MCS events

# ----------------------------------------------------------

# Set file label for type of MCS
if MCS_type == 'all_MCS':
    MCS_file_label = 'MCS'
elif MCS_type == 'robustMCS':
    MCS_file_label = 'robustMCS'
else:
    print('MUST choose valid MCS_type')
    
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
    

# Input paths
general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/firstStorms_InitialGrowthRate_shear_direction_MCS_centroids_plus0.5deg'

specific_inpath = '/data/'

# Load MCS initiation data
MCS_dt = pickle.load(open(general_path + specific_inpath + "%s%s_MCS_dt%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_centroid_elevation = pickle.load(open(general_path + specific_inpath + "%s%s_MCS_centroid_elevation%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_center_lons_initiation_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_MCS_center_lons%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_center_lats_initiation_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_MCS_center_lats%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_ccs_area_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_ccs_area_growth%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_growth_stage_time_length_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_growth_stage_time_length%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))


# Filter events
if events_removed == True:

    events_bool = np.full(len(MCS_dt), True)
    remove_events_nums = [1,7,10,154,72,73,81,84,88,93,95,99,136,141,142,153,57]
    events_bool[remove_events_nums] = False
    events_removed_label = '_events_removed'
    
else:
    
    events_bool = np.full(len(MCS_dt), True) 
    events_removed_label = ''

MCS_dt = np.array(MCS_dt)[events_bool]
MCS_centroid_elevation = np.array(MCS_centroid_elevation)[events_bool]
MCS_center_lons_initiation_filtered = np.array(MCS_center_lons_initiation_filtered)[events_bool]
MCS_center_lats_initiation_filtered = np.array(MCS_center_lats_initiation_filtered)[events_bool]
MCS_ccs_area_growth_filtered = np.array(MCS_ccs_area_growth_filtered)[events_bool]
MCS_growth_stage_time_length_filtered = np.array(MCS_growth_stage_time_length_filtered)[events_bool]

MCS_ccs_area_growth_rate_filtered = MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered
median_MCS_ccs_area_growth_rate_filtered = np.nanmedian(MCS_ccs_area_growth_rate_filtered)

### plot histogram by first storms time, seperate by growth rate, if chosen ###

fig, ax = plt.subplots()

if seperation_plot == True:
    
    ############ Change variable for seperation ############

    var_seperation = MCS_ccs_area_growth_rate_filtered
    var_seperation_name = 'ccs_area_growth_rate'
    seperation_threshold = median_MCS_ccs_area_growth_rate_filtered

    ############### Seperate by variable ################
    
    var_seperation_label = '_%s_by_%s%d' %(growth_type, var_seperation_name, seperation_threshold)

    if growth_type == 'rapid':
        mask = var_seperation >= seperation_threshold
    elif growth_type == 'slow':
        mask = var_seperation < seperation_threshold
    else:
        print('Please enter valid growth type')

    MCS_dt_mask = MCS_dt[mask]
    MCS_centroid_elevation_mask = MCS_centroid_elevation[mask]
    MCS_center_lons_initiation_filtered_mask = MCS_center_lons_initiation_filtered[mask]
    MCS_center_lats_initiation_filtered_mask = MCS_center_lats_initiation_filtered[mask]
    
else: # seperation_plot == False

    MCS_dt_mask = MCS_dt[:]
    MCS_centroid_elevation_mask = MCS_centroid_elevation[:]
    MCS_center_lons_initiation_filtered_mask = MCS_center_lons_initiation_filtered[:]
    MCS_center_lats_initiation_filtered_mask = MCS_center_lats_initiation_filtered[:]
    
    var_seperation_label = ''
    
fig, axs = plt.subplots()

hours = [dt.hour + dt.minute/60 for dt in MCS_dt_mask]

print('hours', hours)

print('sorted hours', sorted(hours))

plt.hist(hours, bins=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24], alpha=0.3)

print('# of MCS', len(MCS_center_lons_initiation_filtered_mask))

# ---------- plot number of MCS  ----------

plt.text(.01, .99, 'n = %d' %(len(MCS_center_lons_initiation_filtered_mask)), ha='left', va='top', transform=axs.transAxes, fontsize=30)

general_outpath = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/firstStorms_InitialGrowthRate_shear_direction_MCS_centroids_plus0.5deg'

specific_outpath = '/data/'

plt.savefig(general_outpath + '/5%s%s%s_diurnal_MCSfirstStorms_%s%s.png' %(MCS_file_label, offset_label, var_seperation_label, MCS_init_area, events_removed_label), dpi=600)

print('saved')

plt.close()
    