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
from sklearn.linear_model import LinearRegression # to fit linear systems to data
from scipy.optimize import curve_fit # to fit a double linear to data
import tikzplotlib # for saving figures to tikz files -- remove at will


#################################################
###### FUNCTIONS FOR ESTIMATION TECHNIQUES ######
#################################################


########################
# Gregory-like methods #
########################

def estimationGregory(DT,DR):
    # Use a linear 'Gregory' regression on DT and DR to fit the linear
    # function DR = a DT + f to the data. Estimated equilibrium warming
    # is then given by DT_*^est = - f / a. Outputted is an array of 
    # estimated equilibrium warmings as function of t, where only data
    # up to time t has been used.
    
    DT_ests = np.empty(len(DT))
    DT_ests[:] = np.nan
    
    for i in range(1,len(DT)):
        mu, F = np.polyfit(DT[:i],DR[:i],1) # Linear fit
        DT_est = - F / mu
        DT_ests[i] = DT_est
        
    return DT_ests
    
    
def estimationGregory20(DT,DR):
    # Use a linear 'Gregory' regression on DT and DR to fit the linear
    # function DR = a DT + f to the data, but ignoring the first 20
    # years of data. Estimated equilibrium warming is then given by 
    # DT_*^est = - f / a. Outputted is an array of estimated equilibrium 
    # warmings as function of t, where only data up to time t has been used.
    
    DT_ests = np.empty(len(DT))
    DT_ests[:] = np.nan
    
    i_start = 20 # How many years to ignore at the start
    
    for i in range(i_start+1,len(DT)):
        mu, F = np.polyfit(DT[i_start:i],DR[i_start:i],1) # Linear fit
        DT_est = - F / mu
        DT_ests[i] = DT_est
        
    return DT_ests
    
    
def double_linear(x, l1, l2, f1, f2):
    # Help function: a double linear function.
    # This function is first linear according to
    # y1 = f1 + l1 * x (for x < x_switch)
    # Then linear acording to
    # y2 = f2 + l2 * x (for x < x_switch)
    # We want to stipulate continuity everywhere, so we need to have
    # y1(x_switch) = y2(x_switch)
    # This yields x_switch = - (f2-f1)/(l2-l1)
    
    x_switch = - (f2-f1)/(l2-l1)
    y = (f1 + l1 * x) * (x < x_switch) + (f2 + l2 * x) * (x >= x_switch)
    return y

def estimationDoubleGregory(DT,DR):
    # Uses a double linear ('Double Gregory') fit on DT and DR.
    # Uses the slope and intercept of the latter part of this 
    # double linear to estimate equilibrium warming.
    
    DT_ests = np.empty(len(DT))
    DT_ests[:] = np.nan
    
    for i in range(4,len(DT)): # does not work with too few data points
        mu_s, F_s = np.polyfit(DT[:4],DR[:4],1) # Initial guess for initial linear
        mu_e, F_e = np.polyfit(DT[:i],DR[:i],1) # Initial guess for final linear
        try:
            popt, pcov = curve_fit(double_linear, DT[:i], DR[:i], [mu_s,mu_e,F_s,F_e])
            DT_est = - popt[3] / popt[1]
        except: # Then the curve fit algorithm could not fit the parameters so we put nan
            DT_est = np.nan
        DT_ests[i] = DT_est
        
    return DT_ests
    
  
###########################
# System-Fit-like methods #
###########################


def system_fit_method(X,Y):
    # General function for all system fit methods. Using same-size arrays X and Y
    # linear regression to Y = A X + F is performed. Different kind of linear
    # regression models can be used by changing the model-variables (see documentation
    # of sklearn and specifically sklearn.linear_model). As found matrix A is quite
    # small taking the inverse is OK time-wise.
    # Input:
    # X: m x N data array of m observables (that get predicted in equilibrium)
    # Y: m x N data array of m observables (that tend to 0 in equilibrium)
    # Output:
    # estimated values for all observables
    
    model = LinearRegression(normalize = True)
    model.fit(X,Y)
    A = model.coef_
    F = model.intercept_
    X_est = - np.dot( np.linalg.inv(A),F )
    return X_est


def system_fit_2vars(X1,X2, Y1,Y2):
    # Perform system fit with 2 variables (m = 2). They get lumped together in arrays X
    # and Y and estimates are made by calling function 'system_fit_method'. Estimations
    # are made as function of time t, where only data up to time t has been used
    
    X1_ests = np.zeros(len(X1))
    X1_ests[:] = np.nan
    X2_ests = np.zeros(len(X1))
    X2_ests[:] = np.nan
    
    X = np.array([X1,X2]).transpose() # to get in right format
    Y = np.array([Y1,Y2]).transpose()
    
    for i in range(5,len(X1)): # Does not work with too few data points
        X_est = system_fit_method(X[:i],Y[:i])
        X1_ests[i] = X_est[0]
        X2_ests[i] = X_est[1]
    
    return X1_ests, X2_ests


def system_fit_3vars(X1,X2,X3, Y1,Y2,Y3):
    # Perform system fit with 3 variables (m = 3). They get lumped together in arrays X
    # and Y and estimates are made by calling function 'system_fit_method'. Estimations
    # are made as function of time t, where only data up to time t has been used
    
    X1_ests = np.zeros(len(X1))
    X1_ests[:] = np.nan
    X2_ests = np.zeros(len(X2))
    X2_ests[:] = np.nan
    X3_ests = np.zeros(len(X3))
    X3_ests[:] = np.nan
    
    X = np.array([X1,X2,X3]).transpose() # to get in right format
    Y = np.array([Y1,Y2,Y3]).transpose()
    
    for i in range(5,len(X1)): # Does not work with too few data points
        X_est = system_fit_method(X[:i],Y[:i])
        X1_ests[i] = X_est[0]
        X2_ests[i] = X_est[1]
        X3_ests[i] = X_est[2]
    
    return X1_ests, X2_ests, X3_ests
    


###################################
###### FUNCTIONS FOR FIGURES ######
###################################

def make_scatterPlot(X,Y, xlabel, ylabel, title):
    # Function that makes a scatter plot with X on horizontal axis and Y on the vertical axis.
    # xlabel ylabel & title will be applied automatically Also handles minimal make-up of the figure.
    
    plt.figure()
    plt.scatter(X,Y, marker = '.', c = 'r')
    plt.xlabel(xlabel, fontsize = 40)
    plt.ylabel(ylabel, fontsize = 40)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(title, fontsize = 40)
    
    return
    
def make_scatterPlots(DT,DALB,DEMM,DR,DALBd,DEMMd, name, coeff_best_fit):
    # Function that construct all possible variations of scatter plots for the model data DT,
    # DALB, DEMM, DR, DALBd and DEMMd.
    
    # Create save directory if it does not exist yet
    directory = 'longrunmip_figs/' + name
    if not os.path.isdir(directory):
        os.mkdir(directory)
    
    
    ########### DR 
    
    ## DT vs DR (+ best fitted line)
    make_scatterPlot(DT,DR, '$\Delta T$', '$\Delta R$', name)
    DR_best = coeff_best_fit[1] + coeff_best_fit[0] * DT
    plt.plot(DT,DR_best, 'k:', label = 'Best Fit')
    plt.legend(fontsize = 20)  
    plt.savefig(directory + '/plot_DT_DR.png')
    try:
        if not os.path.isdir(directory + '/tikz'):
            os.mkdir(directory + '/tikz')
        tikzplotlib.save(directory + '/tikz/plot_DT_DR.tikz')
    except:
        pass
    plt.close()
    
    ## DALB vs DR
    make_scatterPlot(DALB, DR, '$\Delta ALB$', '$\Delta R$', name)
    plt.savefig(directory + '/plot_DALB_DR.png')
    try:
        tikzplotlib.save(directory + '/tikz/plot_DALB_DR.tikz')
    except:
        pass
    plt.close()    
    
    ## DEMM vs DR
    make_scatterPlot(DEMM, DR, '$\Delta EMM$', '$\Delta R$', name)
    plt.savefig(directory + '/plot_DEMM_DR.png')
    try:
        tikzplotlib.save(directory + '/tikz/plot_DEMM_DR.tikz')
    except:
        pass
    plt.close()    
    
    
    ########### DALBd 
    
    ## DT vs DALBd
    make_scatterPlot(DT, DALBd, '$\Delta T$', '$d/dt \Delta ALB$', name)
    plt.savefig(directory + '/plot_DT_DALBd.png')
    try:
        tikzplotlib.save(directory + '/tikz/plot_DT_DALBd.tikz')
    except:
        pass
    plt.close()    
    
    ## DALB vs DALBd
    make_scatterPlot(DALB, DALBd, '$\Delta ALB$', '$d/dt \Delta ALB$', name)
    plt.savefig(directory + '/plot_DALB_DALBd.png')
    try:
        tikzplotlib.save(directory + '/tikz/plot_DALB_DALBd.tikz')
    except:
        pass
    plt.close()    
    
    ## DEMM vs DALBd
    make_scatterPlot(DEMM, DR, '$\Delta EMM$', '$d/dt \Delta ALB$', name)
    plt.savefig(directory + '/plot_DEMM_DALBd.png')
    try:
        tikzplotlib.save(directory + '/tikz/plot_DEMM_DALBd.tikz')
    except:
        pass
    plt.close()    
    
    
    ########### DEMMd 
    
    ## DT vs DEMMd
    make_scatterPlot(DT, DEMMd, '$\Delta T$', '$d/dt \Delta EMM$', name)
    plt.savefig(directory + '/plot_DT_DEMMd.png')
    try:
        tikzplotlib.save(directory + '/tikz/plot_DT_DEMMd.tikz')
    except:
        pass
    plt.close()
    
    ## DALB vs DEMMd
    make_scatterPlot(DALB, DEMMd, '$\Delta ALB$', '$d/dt \Delta EMM$', name)
    plt.savefig(directory + '/plot_DALB_DEMMd.png')
    try:
        tikzplotlib.save(directory + '/tikz/plot_DALB_DEMMd.tikz')
    except:
        pass
    plt.close()
    
    ## DEMM vs DEMMd
    make_scatterPlot(DEMM, DEMMd, '$\Delta EMM$', '$d/dt \Delta EMM$', name)
    plt.savefig(directory + '/plot_DEMM_DEMMd.png')
    try:
        tikzplotlib.save(directory + '/tikz/plot_DEMM_DEMMd.tikz')
    except:
        pass
    plt.close()
    
    
    return

def make_estimationPlot(t,ests, name, est_type, DT_best_range):
    # Function plots all estimated values for the estimated variable in est_type if available
    for est_technique in ests:
        try:
            plt.plot(t, est_technique[est_type + '_ests'], est_technique['color'], label = est_technique['name'])
        except:
            pass
    plt.xlabel('last year used in estimate', fontsize = 40)
    plt.ylabel('estimated change' + est_type + '_*^{est}', fontsize = 40)
    plt.title(name)
    plt.legend(fontsize = 20)
    
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    
    # Put bounds on y-axis depending on type of plot
    if est_type == 'DT':
        plt.ylim([0, np.ceil(np.max(DT_best_range))+1])
    elif est_type == 'DEMM':
        plt.ylim([-4*10**-9, -1*10**-9])
    
    # Save fig
    # Create save directory if it does not exist yet
    directory = 'longrunmip_figs/' + name
    if not os.path.isdir(directory):
        os.mkdir(directory)
    
    plt.savefig(directory + '/estimates_' + est_type + '.png')
    try:
        if not os.path.isdir(directory + '/tikz'):
            os.mkdir(directory + '/tikz')
        tikzplotlib.save(directory + '/tikz/estimates_' + est_type + '.tikz')
    except:
        pass
    
    # Also save first 500 years only plot
    plt.xlim([0, 500])
    plt.savefig(directory + '/estimates_' + est_type + '_500.png')
    try:
        if not os.path.isdir(directory + '/tikz'):
            os.mkdir(directory + '/tikz')
        tikzplotlib.save(directory + '/tikz/estimates_' + est_type + '_500.tikz')
    except:
        pass 
    plt.close()
    
    return
    
    
def make_estimationPlots(t, ests, best_est, best_est_range, name):
    # Make all the plots for the estimated values. Automatically uses the right colors,
    # labels and only plots method which have actually estimated the value that is
    # specified. Here, it is set-up to plot estimates for DT, DALB and DEMM. It plots
    # these both for all data points and only for the first 500 years.

    ### estimated DT plots
    plt.figure()
    plt.fill_between([0, t[-1]], [np.min(best_est_range), np.min(best_est_range)], [np.max(best_est_range), np.max(best_est_range)], facecolor = 'black', alpha = 0.5)
    plt.plot([0, t[-1]], [best_est, best_est], 'k:', label = 'best estimate')
    make_estimationPlot(t,ests, name, 'DT', best_est_range)
    
    ## estimated DALB plots
    plt.figure()
    make_estimationPlot(t,ests, name, 'DALB', [])
    
    ## estimated DEMM plots
    plt.figure()
    make_estimationPlot(t,ests, name, 'DEMM', [])    
    
    return

   
#############################################
### COMBINATION FUNCTION 'do_estimations' ###
#############################################

def do_estimations(t,DT, DALB, DEMM, DR, DALBd, DEMMd, name, best_est, best_est_range, coeff_best_fit):
    # function that combines all functions in this python file. It uses
    # the processed model output to compute estimations for the equilibrium
    # warming based on several extrapolation techniques.
    # These are outputted in an array containing for each estimation technique a
    # dictionary that contains the fields 'method', containing a short-hand
    # name for the method, and 'ests' containing the estimated equilibrium
    # values for each time point, and 'color' for their color in the plots

    
    ### Perform estimations

    ests = []
    
    # Method: 'raw time series'
    ests.append({'name': 'raw simulation value', 'color': 'red', 'marker': 'X', 'DT_ests': DT})
    
    # Method: 'Gregory'
    ests.append({'name': 'Gregory', 'color': 'blue', 'marker': 'o', 'DT_ests': estimationGregory(DT,DR)})
    
    # Method: 'Gregory20'
    ests.append({'name': 'Gregory (ignoring y1-20)', 'color': 'yellow', 'marker': 's', 'DT_ests': estimationGregory20(DT,DR)})
    
    # Method: 'Double Gregory'
    ests.append({'name': 'Double Gregory', 'color': 'green', 'marker': 'D', 'DT_ests' : estimationDoubleGregory(DT,DR)})
    
    # Method: 'System Fit [T,ALB]'
    DT_ests, DALB_ests = system_fit_2vars(DT,DALB,DR,DALBd)  
    ests.append({'name': 'System Fit [T,ALB]', 'color': 'cyan', 'marker': '<', 'DT_ests': DT_ests, 'DALB_ests': DALB_ests})
    
    # Method: 'System Fit [T,EMM]'
    DT_ests, DEMM_ests = system_fit_2vars(DT,DEMM,DR,DEMMd)
    ests.append({'name': 'System Fit [T,EMM]', 'color': 'black', 'marker': '>', 'DT_ests': DT_ests, 'DEMM_ests': DEMM_ests})
    
    # Method: 'System Fit [T,ALB, EMM]'
    DT_ests, DALB_ests, DEMM_ests = system_fit_3vars(DT,DALB,DEMM, DR,DALBd, DEMMd)
    ests.append({'name': 'System Fit [T,ALB, EMM]', 'color': 'magenta', 'marker': '^', 'DT_ests': DT_ests, 'DALB_ests': DALB_ests, 'DEMM_ests': DEMM_ests})
    
    
    
    ### Make some figures
    
    # Diagnostic scatter plots for DT, DALB, DEMM, DR, DALBd and DEMMd
    make_scatterPlots(DT,DALB,DEMM,DR,DALBd,DEMMd, name, coeff_best_fit)
    
    # Plot of the evolution of the estimated warming by the estimation methods
    make_estimationPlots(t,ests, best_est, best_est_range, name)
    
    
    ### Save the data on the estimations
    parent_directory = 'longrunmip_figs'
    directory = parent_directory + '/' + name
    if not os.path.isdir(directory + '/data'):
        os.mkdir(directory + '/data')
    # Since numpy arrays cannot be dumped directly to json files...
    ests_save = ests.copy()
    for est in ests_save:
        for x in est.keys():
            if isinstance(est[x], np.ndarray):
                est[x] = est[x].tolist()
    json.dump( ests_save, open(directory + '/data/estimates.json', 'w'))
    
    
    # Return the found estimates
    return ests