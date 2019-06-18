%% Pipeline for analysis of behavioural data of the fMRI experiment.

%% Author: Veith Weilnhammer
%% Last edited: 17.1.2019

clear all
close all

%% Main experiment

%% Find results files

root_dir = '/home/veithweilnhammer/Public/Git/Prereg_IFC_Bistability/Code/Experimental_Scripts/Experimental_Scripts_fMRI/Results/'
cd(root_dir)
observers = dir(fullfile(root_dir, '*_run_1*'))


%%
%% Prepare behavioural data + Apply Exclusion Criteria
%%

Exclusion.Criteria.PhaseDuration = 118; % do not consider blocks in which participants did not perceive perceptual transitions;
Exclusion.Criteria.AveragePhase = 35; % do not consider participants which show perceptual phases longer than 35 sec on average;
Exclusion.Criteria.Pcorrect = 0.75; % do not consider participants who did not achieve an average accuracy above 75% in he

[Pcorrect rBias Frequency PhaseDur Unclear GLM Exclusion] = prepare_data_main(observers, Exclusion);

%Pcorrect(:,Exclusion.participants_to_exclude,:) = [];
%rBias(:,Exclusion.participants_to_exclude,:) = [];
%PhaseDur(:,Exclusion.participants_to_exclude,:) = [];
%Unclear(:,Exclusion.participants_to_exclude,:) = [];


%%
%% Report conventional data
%%

%% P correct
Conditions1 = {'d1'; 'd2'; 'd3'; 'd4'; 'd5'; 'd6'};

Mean_P_Correct = mean(nanmean(Pcorrect,3),2);
SEM_P_Correct = std(nanmean(Pcorrect,3),0,2)/sqrt(size(Pcorrect,2));
T1 = table(Conditions1,Mean_P_Correct);

figure(), set(gcf, 'Color', 'w')
boxplot(nanmean(Pcorrect,3)'), set(gca, 'XTick', [1:7],'XTickLabel', Conditions1), xlabel('Disambiguation'), ylabel('Percent Correct')


%
%% Perform one-way Repeated measures ANOVA: Conventional Data
%%

%% Dependent Variables:

%% Pcorrect
for idx = 1:length(Conditions1)
    Y(:,idx) = nanmean(Pcorrect(idx,:,2:4),3)';
end
t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6),'VariableNames', {'d1', 'd2', 'd3', 'd4', 'd5', 'd6'});
Disambiguation_Level = [1:6]';

rm = fitrm(t,'d1-d6 ~ 1','WithinDesign',Disambiguation_Level)
ranova(rm)

%% rBias
for idx = 1:length(Conditions1)
    Y(:,idx) = nanmean(rBias(idx,:,2:4),3)';
end
t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6),'VariableNames', {'d1', 'd2', 'd3', 'd4', 'd5', 'd6'});
Disambiguation_Level = [1:6]';

rm = fitrm(t,'d1-d6 ~ 1','WithinDesign',Disambiguation_Level)
ranova(rm)

%% Unclear
for idx = 1:length(Conditions1)
    Y(:,idx) = nanmean(Unclear(idx,:,2:4),3)';
end
t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6),'VariableNames', {'d1', 'd2', 'd3', 'd4', 'd5', 'd6'});
Disambiguation_Level = [1:6]';

rm = fitrm(t,'d1-d6 ~ 1','WithinDesign',Disambiguation_Level)
ranova(rm)

%% Pdur
for idx = 1:length(Conditions1)
    Y(:,idx) = nanmean(PhaseDur(idx,:,2:4),3)';
end
t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6),'VariableNames', {'d1', 'd2', 'd3', 'd4', 'd5', 'd6'});
Disambiguation_Level = [1:6]';

rm = fitrm(t,'d1-d6 ~ 1','WithinDesign',Disambiguation_Level)
ranova(rm)
