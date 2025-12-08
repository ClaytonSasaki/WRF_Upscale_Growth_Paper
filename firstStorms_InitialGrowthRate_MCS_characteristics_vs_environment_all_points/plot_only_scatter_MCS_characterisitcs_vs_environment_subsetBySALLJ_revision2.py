#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on fri Sept 30 11:59:01 2022

@author: crs326

Plots distributions of mean/median environmental characteristics at first storms time for all MCS and those seperated by growth rate

"""

# Use non-interactive backend for environments without display (e.g., remote servers)
import matplotlib
matplotlib.use('Agg') 

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pickle
from scipy import stats

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

plot_type = 'violin'

var_plot_str = 'q_flux_850_v' # Possible options: q_850, q_flux_850_v, srh_0_1km, srh_0_3km, srh_0_1km_noLatCorrection, srh_0_3km_noLatCorrection, MUCAPE, MUCIN, LCL, LFC, bulk_shear_0_3km_0, bulk_shear_0_6km_0, bulk_shear_2_6km, MLCAPE_100mb, MLCIN_100mb, rh_700, lapse_rate_3_8km, theta_e_diff_500_3000m

var_label = '850mb v q flux (g $kg^{-1}$ m $s^{-1}$)' # Possible options: '850mb q (g/kg)', srh_0_1km (m2/ms2), srh_0_3km (m2/ms2), 'bulk 1-6 km shear (kts)', 'MUCAPE (J/kg)', '3-8 km shear (kts)', r'$\Delta\theta_{e}$ 500-3000m (K)', '850mb v q flux (g $kg^{-1}$ m $s^{-1}$)'

# Note: to do LFC-LCL, 1) var_plot_str = 'LFC', 2) var_label = 'LFC - LCL (m)', 3) add -LCL to filetype_label under violin plot, and 4) uncomment lines labeled 'For LFC-LCL'

stat_type = 'mean' # Options: mean, meadian

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
median_SALLJ_max_wind_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_median_SALLJ_max_wind_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#median_SALLJ_height_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_median_SALLJ_height_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
var_plot_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_%s_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, var_plot_str, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#LCL_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_LCL_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb")) # For LFC-LCL
#MCS_refl_coverage_grt0_byEvent = pickle.load( open(general_path + specific_inpath + "%s_refl_coverage_grt0_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_refl_coverage_grt30_byEvent = pickle.load( open(general_path + specific_inpath + "%s_refl_coverage_grt30_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))

# Filter events
if events_removed == True:

    events_bool = np.full(len(var_plot_all_points_byEvent), True)
    remove_events_nums = [1,7,10,154,72,73,81,84,88,93,95,99,136,141,142,153,57]
    events_bool[remove_events_nums] = False
    events_removed_label = '_events_removed'
    
else:
    
    events_bool = np.full(len(var_plot_all_points_byEvent), True)
    events_removed_label = ''
    
MCS_prop_area_SALLJ_all_points_byEvent = np.array(MCS_prop_area_SALLJ_all_points_byEvent, dtype="object")[events_bool]
#median_SALLJ_max_wind_all_points_byEvent = np.array(median_SALLJ_max_wind_all_points_byEvent)[events_bool]
MCS_ccs_area_growth_filtered_all_points_byEvent = np.array(MCS_ccs_area_growth_filtered_all_points_byEvent, dtype="object")[events_bool]
MCS_growth_stage_time_length_filtered_all_points_byEvent = np.array(MCS_growth_stage_time_length_filtered_all_points_byEvent, dtype="object")[events_bool]
var_plot_all_points_byEvent = np.array(var_plot_all_points_byEvent, dtype="object")[events_bool]
#LCL_all_points_byEvent = np.array(LCL_all_points_byEvent, dtype="object")[events_bool] # For LFC-LCL

#var_plot_all_points_byEvent = var_plot_all_points_byEvent - LCL_all_points_byEvent # For LFC-LCL

# ---------- Collapse spatial data and calculate stats ----------

# calculate mean/median in variables and save them into 1D array (for MCS and SALLJ characteristics, first value is taken as does not vary spatially)

MCS_prop_area_SALLJ = np.array([arr[0,0] for arr in MCS_prop_area_SALLJ_all_points_byEvent])
#median_SALLJ_max_wind = np.array([arr[0,0] for arr in median_SALLJ_max_wind_all_points_byEvent])
MCS_ccs_area_growth_filtered = np.array([arr[0,0] for arr in MCS_ccs_area_growth_filtered_all_points_byEvent])
MCS_growth_stage_time_length_filtered = np.array([arr[0,0] for arr in MCS_growth_stage_time_length_filtered_all_points_byEvent])
if stat_type == 'mean':
    var_plot = np.array([np.nanmean(arr) for arr in var_plot_all_points_byEvent])
elif stat_type == 'median':
    var_plot = np.array([np.nanmedian(arr) for arr in var_plot_all_points_byEvent])
else:
    print('Please choose a valid stat_type')

# calculate growth rate

MCS_ccs_area_growth_rate_filtered = MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered
print(MCS_ccs_area_growth_rate_filtered)
#print(mean_q_850)

# remove events where growth rate is nan (because MCS was found at time 0)

nan_indices = np.isnan(MCS_ccs_area_growth_rate_filtered)

MCS_prop_area_SALLJ = np.delete(MCS_prop_area_SALLJ, nan_indices)
#median_SALLJ_max_wind = np.delete(median_SALLJ_max_wind, nan_indices)
MCS_ccs_area_growth_filtered = np.delete(MCS_ccs_area_growth_filtered, nan_indices)
MCS_growth_stage_time_length_filtered = np.delete(MCS_growth_stage_time_length_filtered, nan_indices)
MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, nan_indices)
var_plot = np.delete(var_plot, nan_indices)

# remove events where chosen var is nan

nan_indices = np.isnan(var_plot)

MCS_prop_area_SALLJ = np.delete(MCS_prop_area_SALLJ, nan_indices)
#median_SALLJ_max_wind = np.delete(median_SALLJ_max_wind, nan_indices)
MCS_ccs_area_growth_filtered = np.delete(MCS_ccs_area_growth_filtered, nan_indices)
MCS_growth_stage_time_length_filtered = np.delete(MCS_growth_stage_time_length_filtered, nan_indices)
MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, nan_indices)
var_plot = np.delete(var_plot, nan_indices)


# remove incorrect values in srh_0_1km_noLatCorrection
if var_plot_str == 'srh_0_1km_noLatCorrection':

    incorrect_indices = np.where(var_plot < -300)
    var_plot = np.delete(var_plot, incorrect_indices)

    MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, incorrect_indices)
    MCS_prop_area_SALLJ = np.delete(MCS_prop_area_SALLJ, incorrect_indices)
    
elif var_plot_str == 'srh_0_1km':

    incorrect_indices = np.where(var_plot < -400)
    var_plot = np.delete(var_plot, incorrect_indices)

    MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, incorrect_indices)
    MCS_prop_area_SALLJ = np.delete(MCS_prop_area_SALLJ, incorrect_indices)
    
elif var_plot_str == 'LFC':

    incorrect_indices = np.where(var_plot> 10000)
    var_plot = np.delete(var_plot, incorrect_indices)

    MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, incorrect_indices)
    MCS_prop_area_SALLJ = np.delete(MCS_prop_area_SALLJ, incorrect_indices)
    
elif var_plot_str == 'theta_e_diff_500_3000m':

    incorrect_indices = np.where(var_plot> 5)
    var_plot = np.delete(var_plot, incorrect_indices)

    MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, incorrect_indices)
    MCS_prop_area_SALLJ = np.delete(MCS_prop_area_SALLJ, incorrect_indices)
    
    incorrect_indices = np.where((var_plot == 4.8719735) | (var_plot == 3.8457193) | (var_plot == 3.532486) | (var_plot == 4.9246745))
    var_plot = np.delete(var_plot, incorrect_indices)

    MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, incorrect_indices)
    MCS_prop_area_SALLJ = np.delete(MCS_prop_area_SALLJ, incorrect_indices)
    
elif var_plot_str == 'q_flux_850_v':
    
    var_plot = var_plot * 1000 # multiple by 1000 to get in units of g/kg from kg/kg
    
else:
    pass


############ Change variable for seperation ############

var_seperation = MCS_ccs_area_growth_rate_filtered # MCS_ccs_area_growth_filtered, MCS_growth_stage_time_length_filtered
var_seperation_name = 'growth_rate' # 'growth_rate', 'growth_stage_time_length'
seperation_threshold = np.nanmedian(var_seperation)

# create masked array based on MCS css area growth rate and SALLJ presence

rapid_growth_condition = var_seperation >= seperation_threshold
slow_growth_condition = var_seperation < seperation_threshold
SALLJ_condition = MCS_prop_area_SALLJ >= .3
noSALLJ_condition = MCS_prop_area_SALLJ < .3

mask_rapid_SALLJ = rapid_growth_condition & SALLJ_condition
mask_rapid_noSALLJ = rapid_growth_condition & noSALLJ_condition
mask_slow_SALLJ = slow_growth_condition & SALLJ_condition
mask_slow_noSALLJ = slow_growth_condition & noSALLJ_condition

var_plot_rapid_growth_and_SALLJ = var_plot[mask_rapid_SALLJ]
var_plot_rapid_growth_and_noSALLJ = var_plot[mask_rapid_noSALLJ]
var_plot_slow_growth_and_SALLJ = var_plot[mask_slow_SALLJ]
var_plot_slow_growth_and_noSALLJ = var_plot[mask_slow_noSALLJ]

print('mask_rapid_SALLJ', mask_rapid_SALLJ)
print('mask_rapid_noSALLJ', mask_rapid_noSALLJ)
print('mask_slow_SALLJ', mask_slow_SALLJ)
print('mask_slow_noSALLJ', mask_slow_noSALLJ)

# calculate numbers in each growth rate group

num_MCS_all_growth = len(MCS_prop_area_SALLJ)
num_MCS_rapid_growth_SALLJ = len(var_plot_rapid_growth_and_SALLJ)
num_MCS_rapid_growth_noSALLJ = len(var_plot_rapid_growth_and_noSALLJ)
num_MCS_slow_growth_SALLJ = len(var_plot_slow_growth_and_SALLJ)
num_MCS_slow_growth_noSALLJ = len(var_plot_slow_growth_and_noSALLJ)

#print('MCS of ALL growth rate, n = %d' %(num_MCS_all_growth))
#print('MCS growth rate >= %.0f km^2/hr, n = %d' %(seperation_threshold, num_MCS_rapid_growth))
#print('MCS growth rate < %.0f km^2/hr, n = %d' %(seperation_threshold, num_MCS_slow_growth))

# -------------------- Plotting --------------------

fig, ax = plt.subplots()


if plot_type == 'scatter':

    #------ MCS Length Growth vs 850hPa q -------

    #x_var = mean_q_850
    #y_var = MCS_majoraxislength_growth_filtered
    #
    #filename_x = 'mean_q_850'
    #filename_y = 'MCS_majoraxislength_growth_filtered'
    #
    #xlabel = 'mean q at 850hPa'
    #ylabel = 'MCS_majoraxislength_growth (km/h)'

    #------ MCS Area Growth vs 850hPa q -------

    #x_var = mean_q_850
    #y_var = MCS_ccs_area_growth_filtered
    #
    #filename_x = 'mean_q_850'
    #filename_y = 'MCS_ccs_area_growth_filtered'
    #
    #xlabel = 'mean q at 850hPa'
    #ylabel = 'MCS_ccs_area_growth_filtered (km2/h)'

    #------ MCS Area Growth Rate vs 850hPa q -------

#    x_var = mean_q_850
#    y_var = MCS_ccs_area_growth_rate_filtered
#
#    filename_x = 'mean_q_850'
#    filename_y = 'MCS_ccs_area_growth_rate_filtered'
#
#    xlabel = 'mean q at 850hPa'
#    ylabel = 'ccs growth rate (km2/hr)'

    #------ MCS Area Growth Rate vs 850hPa q flux v -------

    x_var = mean_q_flux_850_v
    y_var = MCS_ccs_area_growth_rate_filtered

    filename_x = 'mean_q_flux_850_v'
    filename_y = 'MCS_ccs_area_growth_rate_filtered'

    xlabel = 'mean q v flux at 850hPa'
    ylabel = 'ccs growth rate (km2/hr)'

    #------ MCS Area Growth Rate vs MUCAPE -------

#    x_var = median_MUCAPE
#    y_var = MCS_ccs_area_growth_rate_filtered
#    
#    filename_x = 'median_MUCAPE'
#    filename_y = 'MCS_ccs_area_growth_rate_filtered'
#    
#    xlabel = 'median MUCAPE'
#    ylabel = 'ccs growth rate (km2/hr)'

    #------ MCS Length Growth vs SALLJ Coverage -------

    #x_var = MCS_prop_area_SALLJ
    #y_var = MCS_majoraxislength_growth_filtered
    #
    #filename_x = 'MCS_prop_area_SALLJ'
    #filename_y = 'MCS_majoraxislength_growth_filtered'
    #
    #xlabel = 'Coverage proportion w/ SALLJ'
    #ylabel = 'MCS_majoraxislength_growth (km/h)'

    #------ MCS Area Growth 2hr vs SALLJ Coverage -------

    #x_var = MCS_prop_area_SALLJ
    #y_var = MCS_ccs_area_growth_filtered_2hr
    #
    #filename_x = 'MCS_prop_area_SALLJ'
    #filename_y = 'MCS_ccs_area_growth_filtered_2hr'
    #
    #xlabel = 'Coverage proportion w/ SALLJ'
    #ylabel = 'MCS_ccs_area_growth_filtered_2hr (km2/2h)'

    #------SALLJ Coverage vs Shear Strength -------

    #x_var = mean_bulk_shear_2_6km
    #y_var = MCS_prop_area_SALLJ
    #
    #filename_x = 'mean_bulk_shear_2_6km'
    #filename_y = 'MCS_prop_area_SALLJ'
    #
    #xlabel = 'Mean 2-6 km shear (kts)'
    #ylabel = 'Coverage proportion w/ SALLJ'

    #------ MCS Length Growth vs Shear Strength -------

    #x_var = mean_bulk_shear_2_6km
    #y_var = MCS_majoraxislength_growth_filtered
    #
    #filename_x = 'mean_bulk_shear_2_6km'
    #filename_y = 'MCS_majoraxislength_growth_filtered'
    #
    #xlabel = 'Mean 2-6 km shear (kts)'
    #ylabel = 'MCS_majoraxislength_growth (km/h)'

    #------ MCS Duration vs Shear Strength -------

    #x_var = mean_bulk_shear_0_3km
    #y_var = MCS_duration_filtered
    #
    #filename_x = 'mean_bulk_shear_0_3km'
    #filename_y = 'MCS_duration_filtered'
    #
    #xlabel = 'Mean 0-3 km shear (kts)'
    #ylabel = 'MCS duration (h)'

    #------ MCS Area Growth Rate vs Shear Strength -------

#    x_var = median_bulk_shear_2_6km
#    y_var = MCS_ccs_area_growth_rate_filtered
#    
#    filename_x = 'median_bulk_shear_2_6km'
#    filename_y = 'MCS_ccs_area_growth_rate_filtered'
#    
#    xlabel = 'median 2-6 km shear (kts)'
#    ylabel = 'ccs growth rate (km2/hr)'

    #------ Shear Strength vs SALLJ Stength (for only SALLJ times) ------- # NOTE: must turn off second scatter where c=median_SALLJ_max_wind below and change first scatter to c=MCS_prop_area_SALLJ_onlySALLJ

    #indices_SALLJ = np.where(~np.isnan(median_SALLJ_max_wind))
    #
    #mean_bulk_shear_2_6km_onlySALLJ = mean_bulk_shear_2_6km[indices_SALLJ]
    #median_SALLJ_max_wind_onlySALLJ = median_SALLJ_max_wind[indices_SALLJ]
    #MCS_prop_area_SALLJ_onlySALLJ = MCS_prop_area_SALLJ[indices_SALLJ]
    #
    #x_var = median_SALLJ_max_wind_onlySALLJ
    #y_var = mean_bulk_shear_2_6km_onlySALLJ
    #
    #filename_x = 'median_SALLJ_max_wind_onlySALLJ'
    #filename_y = 'mean_bulk_shear_2_6km_onlySALLJ'
    #
    #xlabel = 'median_SALLJ_max_wind (kts)'
    #ylabel = 'Mean 2-6 km shear (kts)'
    #
    #ax.set_aspect(0.5)

    # ----------------- ADD SALLJ --------------------

    if add_SALLJ == True:

        ############################## Color by SALLJ coverage ###############################

        plotted_fig = ax.scatter(x_var,y_var, c=MCS_prop_area_SALLJ, cmap='Reds', zorder=2) # MCS_prop_area_SALLJ_onlySALLJ

        cbar = fig.colorbar(plotted_fig)
        cbar.ax.set_ylabel('Proportion SALLJ')


        ############################## Outline color by SALLJ strength ###############################
        #cc = []
        #
        #for SALLJ_max in median_SALLJ_max_wind:
        #    if SALLJ_max >= 20 and SALLJ_max < 30:
        #        cc.append('green')
        #    elif SALLJ_max >= 30 and SALLJ_max < 40:
        #        cc.append('orange')
        #    elif SALLJ_max >= 40:
        #        cc.append('red')
        #    else:
        #        cc.append('black')

        ax.scatter(x_var, y_var, marker='o', s=90, c=median_SALLJ_max_wind, cmap='Blues', zorder=1)

        SALLJ_added_label = '_SALLJ_added'

    else:

        ax.scatter(x_var, y_var, marker='o', s=90)

        SALLJ_added_label = ''

    # ----------------- ADD STATISTICS --------------------

    if add_statistics == True:

        #######################  polynomial fit ###################### 
        ## Fit a polynomial of degree 3
        #degree = 1
        #coeffs = np.polyfit(x_var, y_var, degree)
        #
        ## Create a polynomial function from the coefficients
        #poly_fit = np.poly1d(coeffs)
        #
        ## Generate values for the fitted curve
        #x_fit = np.linspace(min(x_var), max(x_var), 100)
        #y_fit = poly_fit(x_fit)

        ####################### linear regression ###################### 

        # Perform linear regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_var, y_var)

        # R-squared value
        r_squared = r_value**2

        # Generate the regression line
        x_regress = np.linspace(min(x_var), max(x_var), 100)
        y_regress = slope * x_regress + intercept

        # Plot fitted line
        plt.plot(x_regress, y_regress, color='red', linewidth=2)

        print('r_squared', r_squared)

        plt.text(0.1, 0.9, f'$R^2 = {r_squared:.4f}$')

        stats_added_label = '_stats_added'

    else:

        stats_added_label = ''

    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    
    filetype_label = '%s_%s_vs_%s_scatter' %(stat_type, filename_y, filename_x)
    
    flags = '%s%s%s%s' %(filter_label, events_removed_label, SALLJ_added_label, stats_added_label)

    
elif plot_type == 'violin':
    
#    print('sorted rapid', np.sort(var_plot_rapid_growth))
#    print('sorted slow', np.sort(var_plot_slow_growth))
#    
#    # run Mann–Whitney nonparametric test and plot p value on plot
#    _, p = stats.mannwhitneyu(var_plot_rapid_growth, var_plot_slow_growth)
#    
#    if p >= .001:
#        ax.text(0.82, 0.96, 'p = %.3f' %(p), transform=ax.transAxes, fontsize=12)
#    else:
#        ax.text(0.82, 0.96, 'p = %.e' %(p), transform=ax.transAxes, fontsize=12)
    
    data_list = [var_plot, var_plot_rapid_growth_and_SALLJ, var_plot_rapid_growth_and_noSALLJ, var_plot_slow_growth_and_SALLJ, var_plot_slow_growth_and_noSALLJ]
    color_list = ['gray', 'darkblue', 'blue', 'darkred', 'red']
    
    #legend_labels = [(mpatches.Patch(color='gray', alpha=0.2), 'all MCS'), (mpatches.Patch(color='blue', alpha=0.2), 'rapid growth MCS'), (mpatches.Patch(color='red', alpha=0.2), 'slow growth MCS')]
    
    x_labels = ['all MCS', 'rapid \nSALLJ', 'rapid \nnoSALLJ', 'slow \nSALLJ', 'slow \nnoSALLJ']
    
    alpha = 0.2
    
    def adjacent_values(vals, q1, q3):
        upper_adjacent_value = q3 + (q3 - q1) * 1.5
        upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

        lower_adjacent_value = q1 - (q3 - q1) * 1.5
        lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
        return lower_adjacent_value, upper_adjacent_value
    
    def set_axis_style(ax, labels):
        ax.set_xticklabels(x_labels)
        ax.set_xticks(np.arange(0,len(x_labels)))
        
    ax.axhline(y=0, color='lightgray', linestyle='--', alpha=0.5) # Add a dashed horizontal line at y=0
    
    for i, data in enumerate(data_list):
    
        violin = plt.violinplot(data, positions=[i], showmeans=False, showmedians=False, showextrema=False)
        
        for pc in violin['bodies']:
            pc.set_facecolor(color_list[i])
            #pc.set_edgecolor('black')
            pc.set_alpha(alpha)
            
        #violin['cmedians'].set_color(color_list[i])
        #violin['cmedians'].set_linestyle('--') # Set the line style to dashed
        
        # shows quartiles and median, comment out lines above (violin['cmedians']) and set 'showmedians' above to False to use 
        quartile1, median, quartile3 = np.percentile(data, [25, 50, 75])
        #mean = np.mean(data)
        whiskers = np.array([adjacent_values(data, quartile1, quartile3)])
        whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
        
        #ax.hline(median, i-.3, i+.3, color=color_list[i], s=50, zorder=3)
        ax.scatter(i, median, marker='o', color=color_list[i], s=50, zorder=3)
        #ax.scatter(i, mean, marker='*', color='black', s=50, zorder=3)
        ax.vlines(i, quartile1, quartile3, color=color_list[i], linestyle='-', lw=5, alpha=0.8)
        #ax.vlines(1, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)
        
        ax.text(i+.05, quartile1, '%.1f' %(quartile1), horizontalalignment='left', verticalalignment='center', color='black', fontsize=12)
        ax.text(i+.05, median, '%.1f' %(median), horizontalalignment='left', verticalalignment='center', color='black', fontsize=12)
        ax.text(i+.05, quartile3, '%.1f' %(quartile3), horizontalalignment='left', verticalalignment='center', color='black', fontsize=12)

        # adding labels
        plt.ylabel('%s %s' %(stat_type, var_label), fontsize=16)
        set_axis_style(ax, x_labels)
        
        # adding legend
        #plt.legend(*zip(*legend_labels), loc=2)
        
        filetype_label = '%s_%s_violin_by_%s%d' %(stat_type, var_plot_str, var_seperation_name, seperation_threshold)
        
        flags = '%s%s' %(filter_label, events_removed_label)


plt.tight_layout()

print('saving')

#specific_outpath = '%sarea_%s%s%s/plots/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

plt.savefig(general_path + '%s%s_wValueLabels_wSALLJ_revision2.png' %(filetype_label, flags), dpi=200)

print('saved')
    