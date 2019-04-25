%% example_script

clear all
close all

load('test_input.mat')

only_ambiguous = 1; % set to 1 if you would like to ignore unambiguous overlaps ... 
% this should point you to parts in the model that you could skip


%% run model
subject = 1;
session = 1;
    
if only_ambiguous
   y(u(:,3) ~= 0.5,:) = [];
   u(u(:,3) ~=0.5,:) = [];
end

%% explanation y/u

% y: this is what the participant perceived at every overlap of the
% Lissajous figre

% u: this is the input variable and has several columns
    % u(:,1): this is identical to y, i.e. the perveiced rotation at every
    % overlap
    % u(:,2): this is always 0.5 (initially included to account for bias,
    % you can surely leave it)
    % u(:,3): 0.5 if the stimulation is ambiguous, 0/1 if there is
    % stereo-disparity
    % u(:,4)= timing of the overlap

%% Fitting the model with the tapas toolbox

% When assessing our model, Model 1 to 4 are important. Here, we
% systematically remove parts of the model, regarding the parameter InitPi
% (the inital precision of a percept) and StereoPi (the precision of the
% stereo-disparity). Since you don't have stereodisparity, the most
% appropriate model should be Model 2. You could compare it to Model 1 (the
% null model) using Bayesian model comparison; Model 3 and 4 don't make
% sense in your case, since you do not have a stereodisparity condition. 

% Model 5 - 7 include alternative models of bistable perception to compare
% the model to. Model 8 changes an additional parameter of the model (zeta).
% 


Model{1}.subject{subject}.session{session}=tapas_fitModel(y,u,'tapas_hgf_binary_Lissajous_config1', 'tapas_categorical_config')

Model{2}.subject{subject}.session{session}=tapas_fitModel(y,u,'tapas_hgf_binary_Lissajous_config2', 'tapas_categorical_config')    

%Model{3}.subject{subject}.session{session}=tapas_fitModel(y,u,'tapas_hgf_binary_Lissajous_config3', 'tapas_categorical_config')

%Model{4}.subject{subject}.session{session}=tapas_fitModel(y,u,'tapas_hgf_binary_Lissajous_config4', 'tapas_categorical_config')


%Model{5}.subject{subject}.session{session}=tapas_fitModel(y,u,'tapas_hgf_binary_Lehky_config', 'tapas_categorical_Wilson_config')

%Model{6}.subject{subject}.session{session}=tapas_fitModel(y,u,'tapas_hgf_binary_Wilson_config_full', 'tapas_categorical_Wilson_config')

%Model{7}.subject{subject}.session{session}=tapas_fitModel(y,u,'tapas_hgf_binary_Moreno_config', 'tapas_categorical_Moreno_config')

%Model{8}.subject{subject}.session{session}=tapas_fitModel(y,u,'tapas_hgf_binary_Lissajous_config4', 'tapas_categorical_free_config')



