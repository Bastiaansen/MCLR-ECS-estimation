# This file contain python codes to compute the best estimate for the ultimate 
# warming of the model. This is based on information from the last 15% of warming
# following a similar approach to [Rugenstein et al (2019)]. It also includes a
# function to perform a resampling of this data to obtain a range of the possible
# equilibrium warming.



##############################
### IMPORT Python packages ###
##############################

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import os
import tikzplotlib # for saving figures to tikz files -- remove at will

########################
###### RESAMPLING ######
########################

def performResampling(DT,DR, name, model_name):
    # Function that performs a resampling of the datapoints and computes some statistics from 
    # the resampling to get a bit of a spread in the predicted values for DT_* as there can be
    # quite a large uncertainty in the true value -- even with very long runs
    # Also produces a figure of the statistics

    # First perform analysis on the complete data (of the last 15% of warming)
    # Linear 'Gregory' regression
    mu_best, F_best = np.polyfit(DT,DR,1)
    DT_best = - F_best / mu_best
    
    
    N = 10000 # How many times do we resample?
    M = 500 # How many samples do we use each time
    
    DTs = np.zeros(N) # prelocate space for resampling estimates
    DTs[:] = np.nan
    
    for i in range(0,N):
        sample_datapoints = np.random.randint(0,len(DT),M)
        mu, F = np.polyfit(DT[sample_datapoints],DR[sample_datapoints],1)
        DTs[i] = - F / mu
        
    # Statistics on the resampled estimates for DT_*    
    DT_low = np.percentile(DTs,5)
    DT_high = np.percentile(DTs,95)
    
    # Make a figure
    plt.figure()
    plt.hist(DTs,1000)
    plt.vlines(DT_best,0,40, 'r')
    plt.vlines(DT_low,0,40, 'b')
    plt.vlines(DT_high,0,40, 'b')
    plt.xlabel('$DT_*^\mathrm{est,best}$', fontsize = 40)
    plt.ylabel('Numer of occurrences', fontsize = 40)
    plt.title(model_name, fontsize = 40)
    
    # Save figure  
    parent_directory = 'longrunmip_figs'
    directory = parent_directory + '/' + name
    # Create directory if it does not exist yet
    if not os.path.isdir(parent_directory):
        os.mkdir(parent_directory)    
    if not os.path.isdir(directory):
        os.mkdir(directory)
    plt.savefig(directory + '/best_estimation_range.png')
    
    # Convert to tikz
    try:
        if not os.path.isdir(directory + '/tikz'):
            os.mkdir(directory + '/tikz')
        tikzplotlib.save(directory + '/tikz/best_estimation_range.tikz')
    except:
        pass
        
    plt.close() # do not show figure in Notebook
    
    # Save resampling statistics
    if not os.path.isdir(directory + '/data'):
        os.mkdir(directory + '/data')
    DTs = xr.DataArray(DTs)
    DTs.to_netcdf(directory + '/data/best_estimate_resampling.nc')
    
    return DT_best,[DT_low,DT_high], [mu_best, F_best]


###############################################
### COMBINATION FUNCTION 'compute_best_fit' ###
###############################################

def computeBestFit(DT,DR, name, model_name):
    # Computes best fitted value for equilibrium warming based on last 15% of warming
    # and determines a range of possible values for DT_*
    # Outputs best value, range of DT_* and fitted coefficients of linear regression for best estimate

    # Take 15% of the mean remaining warming of the last 50 years of simulation,
    # EXCEPT for in case of the 'IPSCLM5A' model. There, the 1000 years are not nearly enough to be in the
    # last warming, and using only the last 50 years results in a dataset that is almost horizontal, thus
    # leading to a range spanning from values going down to -10^4 and up to +10^4, so that's not useful.
    # Therefore, instead we use the complete dataset (instead of the last 50 years only) to determine the
    # last warming years. In principle, this should lead to a bound that is too low, but in practice, it is
    # still higher then any techniques predict.
    last_warming_R = 0.85 * np.mean(DR[-50:])
    if name == 'IPSLCM5A':
        last_warming_R = 0.85 * np.mean(DR)
    # Find first year that has radiative imbalance below that value
    last_warming = next(ind for ind, val in enumerate(DR) if val < last_warming_R)
    # Fit a linear line through data from that point onward and apply resampling
    DT_star_best, DT_star_bounds, coeff_best = performResampling(DT[last_warming:],DR[last_warming:], name, model_name)
    
    return DT_star_best, DT_star_bounds, coeff_best