%% display power analysis

%% Author: Veith Weilnhammer
%% Last edited: 17.1.2019

clear all
close all

root_dir = '/home/veithweilnhammer/Public/Git/scripts/Matlab/RDK/parametric_RDK/Power_Analysis'
cd(root_dir)

load('pow_results.mat')

samples = find(pow_results.power>95)+10;
minimum_sample_size = min(samples)

figure(), set(gcf, 'Color', 'w')
plot([10:60], pow_results.power), xlabel('sample size'), ylabel('power')