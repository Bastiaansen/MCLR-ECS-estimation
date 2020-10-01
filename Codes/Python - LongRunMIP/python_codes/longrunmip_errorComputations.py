# This file contain python codes to compute estimations for the real model
# equilibrium warming based on data. Estimates are based as a function of
# time t, where only data up to time t has been used for such estimate.
#
# The following estimations have been programmed
# 1: Raw Temperature time series
# 2: Linear 'Gregory' regression on all data
# 3: Linear 'Gregory' regression on data from year 20 onwards
# 4: Double linear fit 'Double Gregory' on all data
# 5: System Fit [DR, DALBd] = A [DT, DALB] + F
# 6: System Fit [DR, DEMMd] = A [DT, DEMM] + F
# 7: System Fit [DR, DALBd, DEMMd] = A [DT, DALB, DEMM] + F
#
# The functions are ultimately all handled via the function 'do_estimations',
# which is the only function that should be included from this python file.
#
# This function also outputs several figures on evolution of estimations,
# regression lines etc.
#
# Output consists of an array containing for each estimation technique a
# dictionary that contains the fields 'method', containing a short-hand
# name for the method, and 'ests' containing the estimated equilibrium
# values for each time point.


##############################
### IMPORT Python packages ###
##############################

import numpy as np
from matplotlib import pyplot as plt
import os
import json
import tikzplotlib # for saving figures to tikz files -- remove at will


#################################################
######   FUNCTIONS FOR ERROR COMPUTATIONS  ######
#################################################

def comp_error_to_range(ests, val_range):
    # Computes the distance between each element in ests and the range in
    # best_est_range. Output is this distance in an array of same size as 
    # ests.
    
    val_min = np.min(val_range)
    val_max = np.max(val_range)
        
    error = np.maximum( ests - val_max, val_min - ests ) # Above or below range
    error = np.maximum( 0, error ) # If value is in the range set error to 0
    
    rel_error = error / ests # relative error
    
    return rel_error

def comp_remaining_error(error):
    # This function computes the remaining error. That is, for every element in
    # the array 'error' the maximum is taken over that value plus all remaining
    # values in the array.
    
    rem_error = np.zeros(len(error))

    for i in range(0, len(error)):
        rem_error[i] = np.max(error[i:])
        
    return rem_error

def comp_errors(estimates, best_est_range):
    # Loop over all models and compute error
    
    remaining_errors = []

    for est_technique in estimates:
        ## Computation of errors
        ests = est_technique['DT_ests']
        errs = comp_error_to_range(ests, best_est_range)
        rem_errs = comp_remaining_error(errs)
        
        ## Put in error array
        remaining_errors.append({'name': est_technique['name'], 'color': est_technique['color'], 'marker': est_technique['marker'], 'rem_err' : rem_errs})
        
    return remaining_errors


######################################
######   FUNCTIONS FOR FIGURES  ######
######################################

def make_errorPlots(t, rem_errs, name):
    # Makes a figure for the remaining error
    
    plt.figure()
    for err in rem_errs:
        plt.loglog(t, err['rem_err'], err['color'], label = err['name'])
    plt.xlabel('last year used in estimate', fontsize = 40)
    plt.ylabel('Remaining relative error', fontsize = 40)
    plt.title(name)
    plt.legend(fontsize = 20)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    
    plt.ylim([10**-3, 10**1])
    plt.xlim(10**1,t[-100])
    
    # Save fig
    # Create save directory if it does not exist yet
    directory = 'longrunmip_figs/' + name
    if not os.path.isdir(directory):
        os.mkdir(directory)
    
    plt.savefig(directory + '/rem_error.png')
    try:
        if not os.path.isdir(directory + '/tikz'):
            os.mkdir(directory + '/tikz')
        tikzplotlib.save(directory + '/tikz/rem_error.tikz')
    except:
        pass
    plt.close()
        
    return
    
############################
## Comparing model errors ##
############################

def compute_comparisonErrors(t, rem_errors, comparison_times):
    # Finds the values for the remaining errors at the given times in the array
    # 'comparison_times'. Output will be an array with a dictionary in each entry
    # that stipulates the remaining error at the given times for each different
    # estimation technique available in 'rem_errors'
    
    comp_rem_errs = []
    for rem_errors_technique in rem_errors:
        rem_errs = []
        t_errs = []
        for comp_t in comparison_times:
            ii = t.index(comp_t)
            rem_errs.append(rem_errors_technique['rem_err'][ii])
            t_errs.append(comp_t)
        comp_rem_errs.append({'name': rem_errors_technique['name'], 'color': rem_errors_technique['color'], 'marker': rem_errors_technique['marker'], 't': t_errs, 'rem_err': rem_errs})
    return comp_rem_errs



############################################
### COMBINATION FUNCTION 'computeErrors' ###
############################################

def computeErrors(t,estimates, best_est_range, name):
    # Function that combines all previous functions. Contains a loop over all
    # estimation techniques in 'estimates', computes errors w.r.t. the best_est_range
    # values for the real equilibrium, saves these and constructs plots.

    ### Compute remaining relative errors
    rem_errs = comp_errors(estimates, best_est_range)
    
    ### Make plots
    make_errorPlots(t, rem_errs, name)
    
    ### Save the data on the estimations
    parent_directory = 'longrunmip_figs'
    directory = parent_directory + '/' + name
    if not os.path.isdir(directory + '/data'):
        os.mkdir(directory + '/data')
    # Since numpy arrays cannot be dumped directly to json files...
    errs_save = rem_errs.copy()
    for err in errs_save:
        for x in err.keys():
            if isinstance(err[x], np.ndarray):
                err[x] = err[x].tolist()
    json.dump( errs_save, open(directory + '/data/relative_errors.json', 'w'))
    
    return rem_errs