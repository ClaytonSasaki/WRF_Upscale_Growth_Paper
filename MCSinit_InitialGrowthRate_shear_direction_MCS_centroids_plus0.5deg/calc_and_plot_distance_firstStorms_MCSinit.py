#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Dec 5 2024

@author: crs326

Calculates the distances between points given arrays of starting and ending latitudes and longitudes. It uses the Haversine formula, which is commonly used to calculate the great-circle distance between two points on the Earth's surface.

"""

# Use non-interactive backend for environments without display (e.g., remote servers)
import matplotlib
matplotlib.use('Agg') 

import matplotlib.pyplot as plt
import numpy as np
import pickle

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points on the Earth using the Haversine formula.
    Parameters are in degrees.
    """
    # Convert latitude and longitude from degrees to radians
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])

    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371  # Radius of the Earth in kilometers
    return c * r

def calculate_distances(start_lats, start_lons, end_lats, end_lons):
    """
    Calculate distances for arrays of coordinates.
    """
    if len(start_lats) != len(start_lons) or len(start_lons) != len(end_lats) or len(end_lats) != len(end_lons):
        raise ValueError("All input arrays must have the same length")

    distances = []
    for lat1, lon1, lat2, lon2 in zip(start_lats, start_lons, end_lats, end_lons):
        distance = haversine(lat1, lon1, lat2, lon2)
        distances.append(distance)
    
    return distances

# ----------------- USER-DEFINED SETTINGS ------------------

filter_label = '_filtered_init_loc_and_start_type'
MCS_type = 'all_MCS'               # Options: 'all_MCS', 'robustMCS'
MCS_init_area = 'large_area1'      # Choose region for MCS initiation
MCS_init_location_filter = True
MCS_start_type_filter = True
offset_MCS_and_conditions = False
hours_offset = -1
events_removed = True

# ----------------------------------------------------------

######### get MCS initation times and centroid locations #########

if MCS_type == 'all_MCS':
    MCS_file_label = 'MCS'
elif MCS_type == 'robustMCS':
    MCS_file_label = 'robustMCS'
else:
    print('MUST choose valid MCS_type')
    
if offset_MCS_and_conditions == True:
    offset_label = '_%dhrPrior' %(abs(hours_offset))
else: #offset_MCS_and_conditions == False
    offset_label = ''

# Input paths
MCSinit_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/MCSinit_InitialGrowthRate_shear_direction_MCS_centroids_plus0.5deg'

MCSinit_specific_inpath = '/data/'

# Load MCS initiation data
MCSinit_center_lons_initiation_filtered = pickle.load(open(MCSinit_path + MCSinit_specific_inpath + "%s%s_MCS_center_lons%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCSinit_center_lats_initiation_filtered = pickle.load(open(MCSinit_path + MCSinit_specific_inpath + "%s%s_MCS_center_lats%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))

# Input paths
general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/firstStorms_InitialGrowthRate_shear_direction_MCS_centroids_plus0.5deg'

specific_inpath = '/data/'

# Load first storm data
MCS_centroid_elevation = pickle.load(open(general_path + specific_inpath + "%s%s_MCS_centroid_elevation%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCSfirstStorms_center_lons_initiation_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_MCS_center_lons%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCSfirstStorms_center_lats_initiation_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_MCS_center_lats%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))

# Filter events
if events_removed == True:

    events_bool = np.full(len(MCS_centroid_elevation), True)
    remove_events_nums = [1,7,10,154,72,73,81,84,88,93,95,99,136,141,142,153,57]
    events_bool[remove_events_nums] = False   
    events_removed_label = '_events_removed'
    
else:
    
    events_bool = np.full(len(MCS_centroid_elevation), True)   
    events_removed_label = ''

MCS_centroid_elevation = np.array(MCS_centroid_elevation)[events_bool]
MCSinit_center_lons_initiation_filtered = np.array(MCSinit_center_lons_initiation_filtered)[events_bool]
MCSinit_center_lats_initiation_filtered = np.array(MCSinit_center_lats_initiation_filtered)[events_bool]
MCSfirstStorms_center_lons_initiation_filtered = np.array(MCSfirstStorms_center_lons_initiation_filtered)[events_bool]
MCSfirstStorms_center_lats_initiation_filtered = np.array(MCSfirstStorms_center_lats_initiation_filtered)[events_bool]

# remame centroid lat / lons for readability
start_lats = MCSfirstStorms_center_lats_initiation_filtered
start_lons = MCSfirstStorms_center_lons_initiation_filtered
end_lats = MCSinit_center_lats_initiation_filtered
end_lons = MCSinit_center_lons_initiation_filtered

distances = calculate_distances(start_lats, start_lons, end_lats, end_lons)
print("Distances (in km):", distances)

print("Sorted Distances (in km):", sorted(distances))

print('Median Distance: ', np.nanmedian(distances))
print('Mean Distance: ', np.nanmean(distances))

### subset by locations: those over the SDC vs those away ###

# define area near the SDC
lon_min = -66.0
lon_max = -64.0
lat_min = -34.0
lat_max = -29.75

# subset by firstStorms locations in chosen area
condition_lat_min = lat_min <= MCSfirstStorms_center_lats_initiation_filtered
condition_lat_max = MCSfirstStorms_center_lats_initiation_filtered <= lat_max
condition_lon_min = lon_min <= MCSfirstStorms_center_lons_initiation_filtered
condition_lon_max = MCSfirstStorms_center_lons_initiation_filtered <= lon_max

mask = condition_lat_min & condition_lat_max & condition_lon_min & condition_lon_max

distances_SDC = np.array(distances)[mask]

print("Sorted Distances SDC (in km):", sorted(distances_SDC))

print('Median Distance SDC: ', np.nanmedian(distances_SDC))
print('Mean Distance SDC: ', np.nanmean(distances_SDC))

### subset by locations: those over 500m vs those below 500m ###

print('MCS_centroid_elevation', MCS_centroid_elevation)

condition_elevation = MCS_centroid_elevation >= 450

mask = condition_elevation

print(mask)

distances_grt500m = np.array(distances)[mask]

distances_grt500m = distances_grt500m[distances_grt500m < 110 ]

print("Sorted Distances >=500m (in km):", sorted(distances_grt500m))

print('Median Distance >=500m: ', np.nanmedian(distances_grt500m))
print('Mean Distance >=500m: ', np.nanmean(distances_grt500m))

### subset by locations: those over 500m vs those below 500m ###

condition_elevation = MCS_centroid_elevation < 450

mask = condition_elevation

distances_less500m = np.array(distances)[mask]

distances_less500m = distances_less500m[distances_less500m > 15]

print("Sorted Distances <500m (in km):", sorted(distances_less500m))

print('Median Distance <500m: ', np.nanmedian(distances_less500m))
print('Mean Distance <500m: ', np.nanmean(distances_less500m))

############## plotting (type 1) #################

#fig, ax = plt.subplots()
#
#data = [distances_less500m,distances_grt500m]
#
#n, bins, patches = plt.hist(x=data, bins=np.arange(0,261,20), color=['#0504aa', 'red'], alpha=0.7, stacked=True, label=['<500m AGL', '>=500m AGL'])
#plt.legend()
#plt.xlabel('Distance between first storms and MCSinit (km)')
#plt.ylabel('Counts')
#
#outpath = MCSinit_path
#
#plt.savefig(outpath + '/5%s%s_hist_distance_MCSfirstStorms_MCSinit_%s%s.png' %(MCS_file_label, offset_label, MCS_init_area, events_removed_label), dpi=600)
#
#print('saved')
#
#plt.close()

############## plotting (type 2) #################

# Create a figure and a set of subplots
fig, axs = plt.subplots(3, 1, figsize=(4, 4))

distances = np.array(distances) # convert list to array
distances = distances[distances !=0]

data = [distances,distances_less500m,distances_grt500m]

label_list = ['All', '<500m AGL', '>=500m AGL']

color_list = ['gray', '#0504aa', 'red']

for i in range(len(data)):
    axs[i].hist(x=data[i], bins=np.arange(0,261,20), color=color_list[i], alpha=0.7, label=label_list[i])
    
    median_val = np.nanmedian(data[i])
    axs[i].axvline(median_val, color=color_list[i], linestyle='dashed', linewidth=2, label=f'Median: {median_val:.0f}')
    
    axs[i].legend()
    
    axs[i].set_ylim([0,30])
    axs[i].set_yticks([0,10,20,30])
    axs[i].set_yticklabels([0,10,20,30])
    
axs[-1].set_xlabel('Distance between first storms and MCSinit (km)')
axs[-1].set_ylabel('Counts')
    
outpath = MCSinit_path

plt.tight_layout()

plt.savefig(outpath + '/%s%s_hist_distance_MCSfirstStorms_MCSinit_%s%s.png' %(MCS_file_label, offset_label, MCS_init_area, events_removed_label), dpi=600)

print('saved')

plt.close()

