<h1> MCLR-ECS-estimation </h1>
Multi-Component Linear Regression method to obtain equilibrium climate sensitivity from short transient warming model simulations

This repository contains the MATLAB (2018a), Jupyter (3.0) and Jupyter Notebook files used for the research in the paper 'Improved Equilibrium Climate Sensitivity Estimations using Short Transient Warming Model Simulations' (submitted)

These codes can be used to perform and compare Equilibrium Climate Sensitivity (ECS) estimation methods. In the paper, these tests have been performed on a conceptual 0D Global Energy Budget Model (GEBM) as well as on model simulations from LongRunMIP (longrunmip.org).

Below follows a short description and instruction for the various codes.

<h2> MATLAB - 3 component model </h2>
This folder contains the MATLAB code to use a 0D GEBM that models global mean surface temperature, albedo and emissivity using a stochastic differential equations (SDE). The files 'doSimulationGEM_ensemble.m' and 'doSimulationGEM_once.m' should be run for respectively an ensemble run or a single run of the model, the parameters can be specified in that file. The file 'ThreeComponentGEM.m' contains the code for the actual SDE. Finally, the file 'perform_estimatoins.m' contains code for equilibrium estimation methods that are compared.

<h2> MATLAB - fig1 - DTDR </h2>
This folder contains MATLAB code to perform various estimation techniques on data. To use these, first the code 'makeDTDRfig1_v2.m' should be run and then subsequently the other scripts can be run as desired.

It is necessary to provide your own files for the globally averaged and yearly averaged observables 'tas' (near surface atmospheric temperature) , 'rsdt' (top-of-atmsophere downwelling shortwave radiation), 'rsut' (top-of-atmsophere upwelling shortwave radiation) and 'rlut' (top-of-atmosphere upwelling longwave radiation). Paths to these files should be put into the relevant variables in the 'make_fig1_v2.m' file (in the code section 'Paths'). These data files are available directly from the LongRunMIP database (see longrunmip.org; [Rugenstein et al, 2010, https://doi.org/10.1175/BAMS-D-19-0068.1]).

Methods provided here are:
* 3EXP: fit a sum of three exponentials to data on global mean temperature and radiative imbalance [Proistosecu and Huybers, 2017, https://doi.org/10.1126/sciadv.1602821]
* EBMeps: fit a two-component GEBM with ocean heat uptake [Geoffroy et al, 2013, https://doi.org/10.1175/JCLI-D-12-00195.1]
* Gregory: fits using Gregory-like linear regressions [Gregory et al, 2004, https://doi.org/10.1029/2003GL018747]
* sysFit: the here newly introduced Multi-Component Linear (system) Regression (MC-LR) [Bastiaansen et al, 2020, submitted]

<h2> MATLAB - rationale images </h2>
Straightforward code to produce one of the illustrative example figures for the rationale section in the Supporting Information

<h2> Python - LongRunMIP </h2>
This folder contains all codes necessary to (re)produce the tests made on the models participating in LongRunMIP. The Jupyter Notebook 'LongRunMIP_ECS_estimations.ipynb' should be run to perform the computations and produce some of the figures. This also automatically creates the folder 'longrunmip_figs'; here, subfolders are created for each model that contain png files of the result figures, as well as a folder for tikz code for the figures and another folder for the results in nc and json file formats. Also some overview figures are created along with relevant data. The notebook does not contain the code for the compuations; the actual python files that are used for the computations are collected in the subfolder 'python_codes' ('bestEstimate' contains code to find true equilibrium values, 'errorComputations' contains code to compute errors of the estimations, 'estimations' contains code to produce the estimations and 'loadData' contains code that loads-in the relevant data files; all files are commented and should be self-explanatory.)

The Jupyter notebook 'remaining_error_plot.ipynb' can be used to quickly perform additional visualisations of the overview results comparing methods and models at fixed times using the automatically created files from the 'LongRunMIP_ECS_estimations.ipynb' notebook.

Finally, 'list_abrupt4x.json' is a json file that contains information on the input data the codes should use (formatted as an array of dictionaries, each array element being a new dictionary for a different model). Specifically, the entry 'model' should use the string that should be used in the figures, 'name' should be the internal name for the model (as used in the input file formats and map structure), 'control' should be the name suffix for the control run and 'experiment' the name suffix for the experiment run that one wants to use. In this repositry this file is formatted for use with models in LongRunMIP that have a long enough 'abrupt-4xCO2' run.

<h3> Instructions -- How to use </h3>

* Create a new folder 'longrunmip_data' with another folder 'global'
* For each model one wants to consider create another folder with the internal name of the model
* In each model, put the nc files (for rlut, rsut, rsdt and tas) for both the control and experiment runs (for LongRunMIP data see longrunmip.org)
* Files should be named like this: 'variablename_modelinternalname_experimentname.nc'; for example 'rlut_CESM104_abrupt4x_5900.nc'
* For each model, create a new array entry in the 'list_abrupt4x.json' file with the required information
* Once this has been done, set-up for the input data is complete

<h4> Python Packages </h4>
The notebooks and python files more contain detailed information on their package dependecies. Specifically, the following packes are used

* numpy
* xarray
* os
* warnings
* json
* sys
* matplotlib
* sklear.linear_model
* scipy.optimize
* tikzplotlib

The package 'tikzplotlib' is probably the least-used and least-common package. This package is used to obtain tikz code (https://github.com/pgf-tikz/pgf) that can be used along with LaTeX to produce (vector) graphics. The code has been set-up in such a way that codes should still work if this package is not installed, provided that one comments out all 'import tikzplotlib' lines in all notebooks and python files (all code referencing the package is put into try except blocks).
