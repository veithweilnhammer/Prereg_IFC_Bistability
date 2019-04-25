%% Analysis of Behavioural Pilot Data and Simulation

%% Author: Veith Weilnhammer
%% Last edited: 17.1.2019

clear all
close all

%% Options

Run_Modelling = 0; % if set to 0, this script will load to estimated Models; if set to 1, models will be inverted in this script.
Run_Simulation = 0;  % if set to 0, this script will load an exmample simulation; if set to 1, models will be done.

save_images = 0;

%% Find results files
paper_dir = '/home/veithweilnhammer/Public/Git/Prereg_IFC_BIstability/Manuscript';
root_dir = '/home/veithweilnhammer/Public/Git/Prereg_IFC_BIstability/Code/Experimental_Scripts/Experimental_Scripts_Pilot/Results'
cd(root_dir)

observers = dir(fullfile(root_dir, '*_run_1*'))

%%
%% Prepare behavioural data + Apply Exclusion Criteria
%%

Exclusion.Criteria.PhaseDuration = 58; % do not consider blocks in which participants did not perceive perceptual transitions;
Exclusion.Criteria.AveragePhase = 35; % do not consider participants which show perceptual phases longer than 35 sec on average;
Exclusion.Criteria.Pcorrect = 0.75; % do not consider participants who did not achieve an average accuracy above 75% in the fully diambiguated condition
[Pcorrect rBiasAmb PhaseDurAmb PhaseDurReplay UnclearAmb  UnclearReplay Exclusion] = prepare_data_pilot(observers, Exclusion);

%%
%% Report conventional data
%%

Conditions1 = {'d1'; 'd2'; 'd3'; 'd4'; 'd5'; 'd6'; 'd7'};

Mean_P_Correct = mean(nanmean(Pcorrect,3),2);
SEM_P_Correct = std(nanmean(Pcorrect,3),0,2)/sqrt(size(Pcorrect,2));

Mean_Bias = mean(nanmean(rBiasAmb,3),2) - 0.5;
SEM_Bias = std(nanmean(rBiasAmb,3),0,2)/sqrt(size(rBiasAmb,2));

T1 = table(Conditions1,Mean_P_Correct, SEM_P_Correct, Mean_Bias, SEM_Bias);

figure(), set(gcf, 'Color', 'w')
boxplot(nanmean(Pcorrect,3)'), set(gca, 'XTick', [1:7],'XTickLabel', Conditions1), xlabel('Disambiguation'), ylabel('Percent Correct'), ylim([0.25 1.01])
if save_images, plot2svg(fullfile(paper_dir, 'percent_correct_human.svg')); end

Conditions2 = {'a'; 'd1'; 'd2'; 'd3'; 'd4'; 'd5'; 'd6'; 'd7'};

Mean_PDur = [mean(mean(nanmean(PhaseDurAmb,3)));  nanmean(nanmean(PhaseDurReplay,3),2)];
SEM_PDur = [std(mean(nanmean(PhaseDurAmb,3)))/sqrt(size(PhaseDurAmb,2)); nanstd((nanmean(PhaseDurReplay,3)),0,2)/sqrt(size(PhaseDurReplay,2))];

Mean_Unclear = [mean(mean(nanmean(UnclearAmb,3)));  nanmean(nanmean(UnclearReplay,3),2)];
SEM_Unclear = [std(mean(nanmean(UnclearAmb,3)))/sqrt(size(UnclearAmb,2)); nanstd((nanmean(UnclearReplay,3)),0,2)/sqrt(size(UnclearReplay,2))];

T2 = table(Conditions2,Mean_PDur, SEM_PDur, Mean_Unclear , SEM_Unclear);


%
%% Perform one-way Repeated measures ANOVA: Conventional Data
%%
%% Dependent Variables:

%% Pcorrect
for idx = 1:length(Conditions1)
    Y(:,idx) = nanmean(Pcorrect(idx,:,:),3)';
end
t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6), Y(:,7),'VariableNames', {'d1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7'});
Disambiguation_Level = [1:7]';

rm = fitrm(t,'d1-d7 ~ 1','WithinDesign',Disambiguation_Level)
ranova(rm)

%% Mean Phase Duration

Y = [mean(nanmean(PhaseDurAmb,3)); nanmean(PhaseDurReplay,3)]';
t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6), Y(:,7), Y(:,8),'VariableNames', {'a','d1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7'});
Sensory_Evidence = [1:8]';

rm = fitrm(t,'a-d7 ~ 1','WithinDesign',Sensory_Evidence)
ranova(rm);

%% Frequency Unclear

Y = [mean(nanmean(UnclearAmb,3)); nanmean(UnclearReplay,3)]';
t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6), Y(:,7), Y(:,8),'VariableNames', {'a','d1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7'});
Sensory_Evidence = [1:8]';

rm = fitrm(t,'a-d7 ~ 1','WithinDesign',Sensory_Evidence)
ranova(rm);

%% Modelling

%% Run modelling
if Run_Modelling
    for idx = 1:length(observers) %% loop over participants
        
        for iidx = 1:size(Pcorrect,3) %% loop over runs
            
            %% load results file
            clear Results
            load([observers(idx).name(1:22) '_run_' num2str(iidx) '.mat'])
            
            for trial = 1:length(Results.discrete_steps) %% loop over blocks
                
                % binary perceptual response: 1 = rightward rotation, 0:
                % leftward rotation
                Results.u{trial}(:,1) = Results.discrete_steps{trial}'; Results.u{trial}(Results.u{trial}(:,1) == -1,1) = 0;
                
                % u(:,2) amount of sensory evidence (O = full ambiguity to 1 =
                % full disambiguation)
                % u(:,3) direction of disambiguation: 0.5 = ambiguous stimulation, 1 =
                % rightward stimulation; 0 = leftward stimulation
                if mod(trial,2)
                    Results.u{trial}(:,2) = zeros(length(Results.discrete_steps{trial}),1);
                    Results.u{trial}(:,3) = ones(length(Results.discrete_steps{trial}),1)-0.5;
                    
                else
                    Results.u{trial}(:,2) = repmat(Results.disambiguation(trial),length(Results.discrete_steps{trial}),1);
                    Results.u{trial}(:,3) = Results.u{trial-1}(:,1);
                    
                end
                
                % timing of the overlap
                Results.u{trial}(:,4) = (Results.TrialStartTime{trial}:(Results.TrialEndTime{trial}-Results.TrialStartTime{trial})/ length(Results.u{trial}):...
                    Results.TrialEndTime{trial}-(Results.TrialEndTime{trial}-Results.TrialStartTime{trial})/ length(Results.u{trial}))-Results.SessionStartTime;
                
                % Perceptual response
                Results.y{trial} = Results.u{trial}(:,1);
                
            end
            
            % blockwise model estimation
            for block = 1:length(Results.u)/2;
                
                clear y u
                
                u = vertcat(Results.u{[block*2-1 block*2]});
                y = u(:,1);
                
                %% Estimate Models to be compared using RFX-BMC
                [idx, iidx, block]
                
                % No Stereodisparity, no stability prior
                Model{1}.subject{idx}.session{iidx}.block{block}=tapas_fitModel(y,u,'tapas_hgf_binary_Lissajous_config1', 'tapas_categorical_config')
                
                % No Stereodisparity, stability prior
                Model{2}.subject{idx}.session{iidx}.block{block}=tapas_fitModel(y,u,'tapas_hgf_binary_Lissajous_config2', 'tapas_categorical_config')
                
                % Stereodisparity, no stability prior
                Model{3}.subject{idx}.session{iidx}.block{block}=tapas_fitModel(y,u,'tapas_hgf_binary_Lissajous_config3', 'tapas_categorical_config')
                
                % Stereodisparity, stability prior
                Model{4}.subject{idx}.session{iidx}.block{block}=tapas_fitModel(y,u,'tapas_hgf_binary_Lissajous_config4', 'tapas_categorical_config')
                
                %% Extract log model evidence for RFX-BMC
                for model = 1:4
                    lme(model, find(Results.disambiguation(block*2) == sort(Results.disambiguation(Results.disambiguation~=0))), idx, iidx) =  Model{model}.subject{idx}.session{iidx}.block{block}.optim.LME;
                    StereoPi(model, find(Results.disambiguation(block*2) == sort(Results.disambiguation(Results.disambiguation~=0))), idx, iidx) =  Model{model}.subject{idx}.session{iidx}.block{block}.p_prc.StereoPi;
                    InitPi(model, find(Results.disambiguation(block*2) == sort(Results.disambiguation(Results.disambiguation~=0))), idx, iidx) =  Model{model}.subject{idx}.session{iidx}.block{block}.p_prc.InitPi;
                end
            end
        end
    end
else
    
    % load pre-existing modelling file
    load(fullfile(root_dir, 'Modelling_Pilot.mat'), 'Model', 'lme', 'StereoPi','InitPi')
    
end


%% Model Analysis

% transfer exclusion criteria on subject level from conventional
% analysis
lme(:,:,Exclusion.participants_to_exclude,:) = [];
StereoPi(:,:,Exclusion.participants_to_exclude,:) = [];
InitPi(:,:,Exclusion.participants_to_exclude,:) = [];

%% Model level inference (RFX-BMC using SPM12-routines)

% We assume fixed effects within participants across blocks (i.e.,
% summing up log-model evidences). spm_BMS_gibb computes
% Random effects BMC across participants. you need SPM12 on your
% matlab - path to

[rfx_all.exp_r,rfx_all.xp,rfx_all.r_samp,rfx_all.g_post]=spm_BMS_gibbs(reshape(sum(sum(lme,2),4),size(lme,3), size(lme,1)))

% Identify optimal model
model_to_extract = find(rfx_all.xp == max(rfx_all.xp));

%% Parameters

%% Extract Parameters and apply block-wise exlusion from conventional analysis
IPS(:,:,:) = log(InitPi(model_to_extract,:,:,:));
DIS(:,:,:) = log(StereoPi(model_to_extract,:,:,:));

IPS(isnan(Pcorrect)) = NaN;
DIS(isnan(Pcorrect)) = NaN;

%% Overview plot of posterior model parameters
figure(), set(gcf, 'Color', 'w')
subplot(2,1,1)
boxplot(nanmean(IPS,3)'), ylabel('Init Pi')
subplot(2,1,2)
boxplot(nanmean(DIS,3)'), ylabel('Stereo Pi'), xlabel('Disambiguation Level')
if save_images, plot2svg(fullfile(paper_dir, 'posterior_IPS_DIS.svg')); end

%% rmANOVA IPS
Y = [nanmean(IPS,3)'];
t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6), Y(:,7),'VariableNames', {'d1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7'});
Sensory_Evidence = [1:7]';

rm = fitrm(t,'d1-d7 ~ 1','WithinDesign',Sensory_Evidence)
ranova(rm);

%% rmANOVA DIS
Y = [nanmean(DIS,3)'];
t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6), Y(:,7),'VariableNames', {'d1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7'});
Sensory_Evidence = [1:7]';

rm = fitrm(t,'d1-d7 ~ 1','WithinDesign',Sensory_Evidence)
ranova(rm);

%% Correlating Conventional with Modelling

%% Perceptual Phase Durations vs. Initial prior stability

% This correlation is calculated across participants
figure(), set(gcf, 'Color', 'w')
scatter(mean(nanmean(PhaseDurAmb,3),1),  nanmean(nanmean(IPS,3),1)), lsline
[rho pval] = corr(mean(nanmean(PhaseDurAmb,3),1)',  nanmean(nanmean(IPS,3),1)')
xlabel(['Average Phase Duration, pval:' num2str(pval) ', rho:' num2str(rho)]), ylabel('Posterior Init Pi')
if save_images, plot2svg(fullfile(paper_dir, 'corr_posterior_IPS.svg')); end
%% Congruent perceptual decisions vs. precision of sensory evidence

% Here, this correlation is calculated within participants; the resulting
% rhos are tested against 0 on the second-level (two-sided one-sample
% t-test)
for idx = 1:size(Pcorrect,2)
    [R(idx) pval] = corr(nanmean(Pcorrect(:,idx,:),3), nanmean((DIS(:,idx,:)),3));
end
[H,P,CI,STAT] = ttest(R)

figure(), set(gcf, 'Color', 'w')
boxplot(R), ylim([0 1])
xlabel(['Rho averaged across particiants, 2-sided one-sample t-test:' num2str(P)]), ylabel('Rho')
if save_images, plot2svg(fullfile(paper_dir, 'corr_posterior_DIS.svg')); end
%%
%% Simulation
n_sim_participants = 47;
if Run_Simulation
    %% Values for simulation
    
    % participants to simulate
    
    
    % Here, we pick parameters for simulation. In analogy to our pre-registered
    % experiment, we simulate blocks of 120 sec in with inter-overlap intervals
    % of 1.5 sec. We consider 6 blocks with increasing levels of sensory
    % evidence.
    
    experimental_duration = 120;
    sampling_rate = 1/1.5;
    
    % random inital precision of stability prior between 0.30 and 0.70 quantile.
    sim_IPS = rand(1,n_sim_participants).*(quantile(nanmean(IPS(:),3),0.70)-quantile(nanmean(IPS(:),3),0.30))+quantile(nanmean(IPS(:),3),0.30)';
    
    for idx = 1:n_sim_participants
        % random inital precision of disambiguating signal between 0.30 and 0.70 quantile.
        sim_DIS(:,idx) = rand(1,6).*(quantile(nanmean(DIS(2:end,:,:),3)',0.70)-quantile(nanmean(DIS(2:end,:,:),3)',0.30))+quantile(nanmean(DIS(2:end,:,:),3)',0.30);
    end
    
    % show parameters for simulation
    figure(), set(gcf, 'Color', 'w'), title('Prior parameters for simulation'), hold on
    
    subplot(2,1,1)
    boxplot(sim_IPS), xlabel('Initial stability precision'), ylabel('Average Prior Stability Precision'), hold on
    errorbar(1,mean(sim_IPS), std(sim_IPS))
    
    subplot(2,1,2)
    boxplot(sim_DIS')
    xlabel('Levels of Disambiguation'), ylabel('Average Prior Precision of disambiguation')
    
    %% run simulation
    
    for idx = 1:n_sim_participants
        
        % ambiguous blocks
        for iidx = 1      % loop over ambiguous runs SR1
            for iiidx = 1:6     % loop over blocks within an ambiguous run
                [SIM_correct(iiidx, iidx,idx) SIM_probability(iiidx, iidx,idx) mean_c_daq(iiidx, iidx,idx) mean_ic_daq(iiidx, iidx,idx) congruency(:,iiidx, iidx,idx) daq(:,iiidx, iidx,idx) pd(:,iiidx, iidx,idx)]  =  ...
                    simulate_data(experimental_duration, sampling_rate, exp(sim_IPS(idx)), 0);
            end
        end
        
        switch_rate = SIM_probability(:, :,idx); switch_rate(switch_rate == 0) = NaN;
        
        % disambiguated blocks
        for iidx = 2:4      % loop over parametrically disambiguated runs SR3 to SR6
            for iiidx = 1:6     % loop over blocks within an ambiguous run
                [SIM_correct(iiidx, iidx,idx) SIM_probability(iiidx, iidx,idx) mean_c_daq(iiidx, iidx,idx) mean_ic_daq(iiidx, iidx,idx) congruency(:,iiidx, iidx,idx) daq(:,iiidx, iidx,idx) pd(:,iiidx, iidx,idx)]  =  ...
                    simulate_data(experimental_duration, sampling_rate, exp(sim_IPS(idx)), exp(sim_DIS(iiidx,idx)), nanmean(nanmean(switch_rate)));
            end
        end
    end
    
else
    % load pre-existing simulation
    load(fullfile(root_dir, 'Simulation_Pilot.mat'))
end

%% Plot simulation

% Pcorrect
figure(), set(gcf, 'Color', 'w')
boxplot(reshape(mean(SIM_correct(:,2:4,:),2),6,n_sim_participants)'), ylim([0.25 1.01])
xlabel('Levels of Disambiguation'), ylabel('Simulated percent correct')
if save_images, plot2svg(fullfile(paper_dir, 'percent_correct_sim.svg')); end

%% rmANOVA
Y = [reshape(mean(SIM_correct(:,2:4,:),2),6,n_sim_participants)'];
t = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6),'VariableNames', {'d1', 'd2', 'd3', 'd4', 'd5', 'd6'});
Sensory_Evidence = [1:6]';

rm = fitrm(t,'d1-d6 ~ 1','WithinDesign',Sensory_Evidence)
ranova(rm);
%%

% Average prediction errors separately for ambiguity, incongruency and
% congruency
figure(), set(gcf, 'Color', 'w')
subplot(3,1,1)
boxplot(reshape(mean(mean(mean(abs(daq(:,:,1,:)),1),2),3),1,n_sim_participants)), xlabel('Ambiguity'), , ylim([0 0.6])
subplot(3,1,2)
boxplot(reshape(mean(mean_c_daq(:,2:4,:),2),6,n_sim_participants)'), xlabel('Congruency: Levels of Disambiguation'), ylabel('Average Prediction Error')
subplot(3,1,3)
boxplot(reshape(mean(mean_ic_daq(:,2:4,:),2),6,n_sim_participants)'), xlabel('Incongruency: Levels of Disambiguation')
if save_images, plot2svg(fullfile(paper_dir, 'model_quantities.svg')); end

%% two_way_repeated_meausures ANOVA repeated_measures ANOVA (code by Aaron Schurger 2005)

clear Y S F1 F2 % Y: dependent variable (PEs), S: subject factor, F1: congruency, F2: level of disparity

counter = 0;
for idx = 1:n_sim_participants %% loop over participants
    
    % incongruent and congruent PEs (runs R3 to R6) relative to PEs during ambiguous runs
    % R1 and R2
    C = nanmean(mean_c_daq(:,2:4,idx),2) - mean(mean(mean(abs(daq(:,:,1,idx)),1),2),3);
    I = nanmean(mean_ic_daq(:,2:4,idx),2) - mean(mean(mean(abs(daq(:,:,1,idx)),1),2),3);
    
    if ~any(isnan(C)) & ~any(isnan(I))
        for iiidx = 1:6
            counter = counter +1;
            
            Y(counter*2-1) = C(iiidx); %dependent variable average Prediction Error
            Y(counter*2) = I(iiidx);
            
            S(counter*2-1) = idx; % Subject index
            S(counter*2) = idx;
            
            F1(counter*2-1) = 1; % Factor 1: congruency
            F1(counter*2) = 2;
            
            F2(counter*2-1) = iiidx; % Factor 2: disparity
            F2(counter*2) = iiidx;
            
        end
    end
end

FACTNAMES = {'congruency', 'disparity'}
stats = rm_anova2(Y',S',F1',F2',FACTNAMES)
