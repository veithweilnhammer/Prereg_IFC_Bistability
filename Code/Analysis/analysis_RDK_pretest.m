%% Pipeline for analysis of behavioural data of the pre-test experiment.

%% Author: Veith Weilnhammer
%% Last edited: 17.1.2019

clear all
close all

%% Main experiment

%% Find results files

root_dir = '/home/veithweilnhammer/Public/Git/Prereg_IFC_Bistability/Code/Experimental_Scripts/Experimental_Scripts_PreTest/Results'
cd(root_dir)

observers = dir(fullfile(root_dir, '*_run_1*'))


%%
%% Prepare behavioural data + Apply Exclusion Criteria
%%

Exclusion.Criteria.PhaseDuration = 118; % do not consider blocks in which participants did not perceive perceptual transitions;
Exclusion.Criteria.AveragePhase = 35; % do not consider participants which show perceptual phases longer than 35 sec on average;
Exclusion.Criteria.Pcorrect = 0.75; % do not consider participants who did not achieve an average accuracy above 75% in he

[Pcorrect rBias Frequency PhaseDur Unclear congruent_Rating incongruent_Rating Exclusion] = prepare_data_pretest(observers, Exclusion);


%%
%% Report conventional data
%%

%% P correct
Conditions1 = {'d1'; 'd2'; 'd3'; 'd4'; 'd5'; 'd6'};

figure(), set(gcf, 'Color', 'w')
subplot(2,1,1)
boxplot(Pcorrect(:,:,3)'), set(gca, 'XTick', [1:7],'XTickLabel', Conditions1), xlabel('Disambiguation'), ylabel('Percent Correct')

subplot(2,1,2)
boxplot(PhaseDur(:,:,3)'), set(gca, 'XTick', [1:7],'XTickLabel', Conditions1), xlabel('Disambiguation'), ylabel('PhaseDur'), ylim([0 40])

%% Perform one-way Repeated measures ANOVA: Conventional Data
%%

%% Dependent Variables:

%% Pcorrect

Y = Pcorrect(:,:,3)';
t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6),'VariableNames', {'d1', 'd2', 'd3', 'd4', 'd5', 'd6'});
Disambiguation_Level = [1:6]';

rm = fitrm(t,'d1-d6 ~ 1','WithinDesign',Disambiguation_Level)
ranova(rm)

%% rBias
    Y = rBias(:,:,3)';

t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6),'VariableNames', {'d1', 'd2', 'd3', 'd4', 'd5', 'd6'});
Disambiguation_Level = [1:6]';

rm = fitrm(t,'d1-d6 ~ 1','WithinDesign',Disambiguation_Level)
ranova(rm)

%% Unclear
Y = Unclear(:,:,3)';
t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6),'VariableNames', {'d1', 'd2', 'd3', 'd4', 'd5', 'd6'});
Disambiguation_Level = [1:6]';

rm = fitrm(t,'d1-d6 ~ 1','WithinDesign',Disambiguation_Level)
ranova(rm)

%% Pdur
Y = PhaseDur(:,:,3)';

t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6),'VariableNames', {'d1', 'd2', 'd3', 'd4', 'd5', 'd6'});
Disambiguation_Level = [1:6]';

rm = fitrm(t,'d1-d6 ~ 1','WithinDesign',Disambiguation_Level)
ranova(rm)
