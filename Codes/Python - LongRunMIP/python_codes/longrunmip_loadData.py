# This file contain python codes to load in data from experiment runs and control runs.
# Moreover, it also contain functions that compute initial values (T0, R0, ALB0, EMM0)
# using the control runs, and anomalies (DT,DR,DALB,DEMM) using the experiment runs and
# these initial values. Also a function to numerically compute the derivative using
# forward differences is included. A combination of all the data processing that is 
# needed for each model is included in the function 'loadData', which should be all 
# that one needs to include in the jupyter notebook

# NOTE: code takes care of an edge-case in the control run for 'GISSE2R' (see load_control)
# NOTE: code takes care of an edge-case for the model 'CCSM3' (see compute_anomalies)
   


##############################
### IMPORT Python packages ###
##############################

import xarray as xr
import numpy as np

#####################################
### FUNCTIONS FOR LOADING IN DATA ###
#####################################


def load_experiment(direct, name, experiment):
    # Load in data for the fields 'rlut', 'rsut', 'rsdt' and 'tas' for the experiment run

    abr_rlut = xr.open_dataarray(direct + '/rlut_' + name + "_" + experiment + ".nc")
    abr_rsut = xr.open_dataarray(direct + '/rsut_' + name + "_" + experiment + ".nc")
    abr_rsdt = xr.open_dataarray(direct + '/rsdt_' + name + "_" + experiment + ".nc")
    abr_tas = xr.open_dataarray(direct + '/tas_' + name + "_" + experiment + ".nc")
    
    return abr_tas, abr_rlut, abr_rsut, abr_rsdt
 
 
def load_control(direct, name, control):
    # loads all datasets for the control run
    
    # Since 'GISSE2R' does have a discrepency, we need to take care of that here
    # The problem is that the fields 'rlut' and 'tas' have data for 5525 years, 
    # and 'rsut' and 'rsdt' only have data for 5225 years. According to the data
    # instruction the last years of the longer sets match with the data from the
    # shorter sets. Hence, for simplicity we simply discard the first 300 years
    # of the longer sets to circumvent issues (it should not matter too much, as
    # the control runs are used only to construct initial values and we then still
    # have >5000 years for that.) See also the LongRunMIP data readme.
    if name == 'GISSE2R':
        ctrl_rlut = xr.open_dataarray(direct + '/rlut_' + name + "_" + control + ".nc")
        ctrl_rsut = xr.open_dataarray(direct + '/rsut_' + name + "_control_5225.nc")
        ctrl_rsdt = xr.open_dataarray(direct + '/rsdt_' + name + "_control_5225.nc")
        ctrl_tas = xr.open_dataarray(direct + '/tas_' + name + "_" + control + ".nc")
        ctrl_rlut = ctrl_rlut[-len(ctrl_rsut):]
        ctrl_tas = ctrl_tas[-len(ctrl_rsut):]
    else: # for all other models just open the normal data sets
        ctrl_rlut = xr.open_dataarray(direct + '/rlut_' + name + "_" + control + ".nc")
        ctrl_rsut = xr.open_dataarray(direct + '/rsut_' + name + "_" + control + ".nc")
        ctrl_rsdt = xr.open_dataarray(direct + '/rsdt_' + name + "_" + control + ".nc")
        ctrl_tas = xr.open_dataarray(direct + '/tas_' + name + "_" + control + ".nc")
    
    return ctrl_tas, ctrl_rlut, ctrl_rsut, ctrl_rsdt
    
    
    
#####################################
### FUNCTIONS FOR PROCESSING DATA ###
#####################################    
    
    
def compute_anomalies(abr_tas, abr_rlut, abr_rsut, abr_rsdt, ctrl_tas, ctrl_rlut, ctrl_rsut, ctrl_rsdt, name):
    # This function computes effective albedo, emissivity and radiative imbalance for experiment
    # and control run. It computes initial values for these fields based on the control run and
    # computes anomalies as difference between experiment run and these initial values.
    
    # Albedo: ALB = rsut / rsdt_
    # Emissivity: EMM = rlut / T^4
    # Imbalance: R = rsdt - rsut - rlut

    ## Computation of T,ALB, EMM and R for control and experiment run
    ctrl_T = ctrl_tas
    ctrl_ALB = ctrl_rsut / ctrl_rsdt
    ctrl_EMM = ctrl_rlut / (ctrl_tas**4)
    ctrl_R = ctrl_rsdt - ctrl_rsut - ctrl_rlut
    
    abr_T = abr_tas
    abr_ALB = abr_rsut / abr_rsdt
    abr_EMM = abr_rlut / (abr_tas**4)
    abr_R = abr_rsdt - abr_rsut - abr_rlut

    ## Computation of initial values and anomalies
    # NOTE: The 'CCSM3' run is out of equilibrium in the first years of both the experiment and the 
    # control runs according to the longrunmip readme file. Hence, for 'CCSM3' we cannot use the mean 
    # of the  whole control run, but will use the year-to-year values for the first 20 years instead, 
    # as stipulated by the longrunmip readme file.
    if name == 'CCSM3':
        DT = np.empty(len(abr_T))
        DR = np.empty(len(abr_R))
        DALB = np.empty(len(abr_ALB))
        DEMM = np.empty(len(abr_EMM))
        # First 20 years
        DT[:20] = abr_T[:20] - ctrl_T[:20]
        DR[:20] = abr_R[:20] - ctrl_R[:20]
        DALB[:20] = abr_ALB[:20] - ctrl_ALB[:20]
        DEMM[:20] = abr_EMM[:20] - ctrl_EMM[:20]
        # Remaining years (here mean over all remaining years is taken)
        DT[20:] = abr_T[20:] - ctrl_T[20:].mean()
        DR[20:] = abr_R[20:] - ctrl_R[20:].mean()
        DALB[20:] = abr_ALB[20:] - ctrl_ALB[20:].mean()
        DEMM[20:] = abr_EMM[20:] - ctrl_EMM[20:].mean()
    else:
        # Initial values
        T0 = ctrl_T.mean()
        R0 = ctrl_R.mean()
        ALB0 = ctrl_ALB.mean()
        EMM0 = ctrl_EMM.mean()
        # Anomalies
        DT = abr_T - T0
        DR = abr_R - R0
        DALB = abr_ALB - ALB0
        DEMM = abr_EMM - EMM0
    
    ## Convert to np.array's (as not all python packages can support xarrays)
    DT = np.array(DT)
    DR = np.array(DR)
    DALB = np.array(DALB)
    DEMM =np.array(DEMM)
    
    return DT, DALB, DEMM, DR
    
    
def comp_der(y):
    # Function computes numerical derivative using forward differences and also
    # output the mean of subsequent entries to have same-size arrays
    der = np.diff(y)
    y_m = (y[:-1]+y[1:])/2
    return der,y_m
  
  
def comp_ders(DT,DALB,DEMM,DR):
    # Takes derivatives of observables and outputs same-size arrays to work with
    [DTd,DT] = comp_der(DT)
    [DRd,DR] = comp_der(DR)
    [DALBd,DALB] = comp_der(DALB)
    [DEMMd,DEMM] = comp_der(DEMM)
    return DT, DTd, DR, DRd, DALB, DALBd, DEMM, DEMMd
    
   
   
#######################################
### COMBINATION FUNCTION 'loadData' ###
####################################### 

def loadData(directory, name, experiment, control):
    # function that combines all previous functions and outputs ready-to-use
    # datasets for t (#years), DT, DALB, DEMM, DR, DALBd and DEMMd
    
    
    abr_tas, abr_rlut, abr_rsut, abr_rsdt = load_experiment(directory, name, experiment)
    ctrl_tas,ctrl_rlut,ctrl_rsut,ctrl_rsdt= load_control(directory,name, control)
    DT, DALB, DEMM, DR = compute_anomalies(abr_tas, abr_rlut, abr_rsut, abr_rsdt, ctrl_tas, ctrl_rlut, ctrl_rsut, ctrl_rsdt, name)   
    
    # Rework datasets to include derivatives and interpolation values
    DT, DTd, DR, DRd, DALB, DALBd, DEMM, DEMMd = comp_ders(DT,DALB,DEMM,DR)
    
    # Could also be catched from the xarray datasets instead; however, not all models store
    # this the same way, so that is much more of a hassle
    # range starts at 2 since we have lost one entry in the arrays because of taking the derivative
    t = range(2,len(DT)+2)
    
    return t, DT, DALB, DEMM, DR, DALBd, DEMMd
    