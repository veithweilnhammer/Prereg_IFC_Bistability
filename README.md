## README
This README-file describes the use of this repository as supplementary information for our pre-registered study “The Role of Inferior Frontal Cortex in Bistable Perception”. Specifically, we provide information on how to install and use our Matlab-scripts for behavioural pilot data, simulations, power analysis, PsychToolbox-based presentation of the experiments as well as fMRI-preprocessing.

## Table of Contents
1. General remarks

2. Analysis of Pilot Data and Simulations

3. Experimental Scripts

3.1. Pilot Experiment

3.2. Main Experiment

3.3. Control Experiment

3.4. Pre-Test Experiment

3.5. Hetereochromatic Flicker Photometry

4. fMRI-preprocessing

5. Power analysis


## 1. General Remarks
This repository contains all materials necessary to review and repeat our experimental procedures, simulations and analyses. With this, we would like to enable reviewers to retrace all steps relevant to our pre-registration “The Role of Inferior Frontal Cortex in Perceptual Bistability”. 

We created all scripts on 18.04 Ubuntu LTS with Matlab 2014b® (all available toolboxes included). Importantly, the experimental scripts make use of the Psychtoolbox (Version 3.0.8, see <http://psychtoolbox.org/> to download and install the toolbox). Furthermore, we use functions from Matlab® Statistics Toolbox to analyze pilot and simulation data. For Random-Effects Bayesian model comparison, we use the SPM12 toolbox (see <https://www.fil.ion.ucl.ac.uk/spm/software/spm12/> to download and install the toolbox). Please do not hesitate to contact us if you miss any function necessary for running our scripts.

To download the repository, you can either clone or download the repository at <https://github.com/veithweilnhammer/Prereg_IFC_Bistability>. An easy way to obtain all files at once is to click on “Clone or download” on the summary page of our repository. You can either clone the repository to your machine using Git and the URL (<https://github.com/veithweilnhammer/Prereg_IFC_Bistability.git>) or save the repository as a ZIP-file and unpack it to your preferred location on your computer.

You should then be able to open all “.m” files using an Editor such as the Matlab®-Editor. We provide any pre-existing data (Results from the pilot experiment, estimated behavioral models and results from the stimulation data described in the manuscript) in the “.mat” format. 


## 2. Analysis of Pilot Data and Simulations
In this section, we describe how you can retrace our analyses of the pilot experiment, invert the predictive coding models of bistable perception, compare models using Bayesian model selection, simulate from the model and analyze simulated data. 

To do this, please go to the repository subfolder “Analysis” in Matlab®, add this folder to your matlab path and open “analysis_RKD_pilot_sim.m”. In the section “Find results files”, you should update the variable “root_dir” according to the location of the repository on your computer. It should point to the location of the Results files in the subfolder “Results/Pilot/”. 

At the top of the script, you can choose whether to invert the predictive coding model of bistable perception yourself or to load the inverted models (“Run_Modelling = 1” vs. “Run_Modelling = 0”). Please note that inverting the models takes time. You can also choose whether to simulate data yourself or to load a set of simulated data. Once you have chosen these options, you can run “analysis_RKD_pilot_sim.m” to get all results and figures at once or go through the script step-by-step.


## 3. Experimental Scripts
In this section, we would like to walk you through the presentation of our experimental paradigms. The relevant scripts are located in the folder “Experimental/Scripts”. Pilot, main and control experiment as well as heterochromatic flicker photometry are located in the respective subfolder. 

In each subfolder, you will find an “experimental_script_X.m”-file to run the respective experiment and a subfolder “Settings” containing a “Settings_*.m”-file, which will allow you to specify a number of parameters relevant to the presentation of the experiment (i.e. timing, monitor settings etc.). Please note that the Settingsfile contains the parameters needed for presentation specific to your equipment. It may be necessary to change them according to your setup (i.e. viewing distance, keyboard, monitor width/resolution/frequency etc.) to allow for optimal viewing.

To run the control experiment, which controls for temporal imbalance in the main fMRI experiment, you need behavioural data from the fMRI-experiment. To this end, we provide a pilot participant ("observer_1_veith").

For debugging purposes (missing functions or other unexpected problems reviewers might encounter when running the psychtoolbox scripts), we limit the experimental screen to an area of 900x900 pixels at the top-left of the primary screen. This can be easily changed by editing the function "presentation_*_RDK.m" in the subfolder "Functions" in each experimental folder. 

To view the experiment full-screen, please change the line "[w, rect] = Screen('OpenWindow', screenNumber, 0,[0 0 900 900], 32, doublebuffer+1, [], 128)" to "[w, rect] = Screen('OpenWindow', screenNumber, 0,[ ], 32, doublebuffer+1, [], 128)".


To conduct an experiment, go to the respective subfolder in Matlab®, remove all previously used experimental folders from your matlab path, add the current experimental folder to our matlab path and run “experimental_script_X.m”. Here is a list of the exact location of the relevant files: 


**3.1. Pilot Experiment**

Experimental_Scripts/Experimental_Scripts_Pilot/experimental_script_RDK.m

Experimental_Scripts/Experimental_Scripts_Pilot/Settings/Settings_RDK.m 

**3.2. Main Experiment**

Experimental_Scripts/Experimental_Scripts_fMRI/experimental_script_pRDK.m

Experimental_Scripts/Experimental_Scripts_fMRI/Settings/Settings_pRDK.m 

**3.3. Control Experiment**

Experimental_Scripts/Experimental_Scripts_fMRI_Control/experimental_script_C_pRDK.m 

Experimental_Scripts/Experimental_Scripts_fMRI_Control/Settings/Settings_C_pRDK.m 

**3.4. PreTest Experiment**

Experimental_Scripts/Experimental_Scripts_fMRI_PreTest/experimental_script_pRDK.m 

Experimental_Scripts/Experimental_Scripts_fMRI_PreTest/Settings/Settings_pRDK_pretest.m 

**3.5. Hetereochromatic Flicker Photometry**

Experimental_Scripts/Heterochromatic/experimental_script_heterochromatic.m

## 4. fMRI-preprocessing
We define our routines for preprocessing of fMRI-images in “Preprocessing/preproc_fmri_vw.m”. This script is a modification of the preprocessing-routines of the SPM12-toolbox (see “preproc_fmri.m”, © 2014 Wellcome Trust Center for Neuroimaging, Ged Ridgway) and builds a matlabbatch-structure, which you can view in the Batch-Editor of SPM12 GUI.

## 5. Power Analysis
Power analysis was performed with fMRI-power within an antomical mask for inferior frontal cortex (i.e., combined neuromorphometrics mask for right-hemispherical anterior insula and inferior frontal gyrus). The mask is available in “/Power_analysis/Neuromorphometrics_mask_IFG_INSULA.nii”. The results from power analysis are stored within “Power_analysis/pow_results.mat”. They can be displayed using “Power_analysis/display_power_results.m”.
