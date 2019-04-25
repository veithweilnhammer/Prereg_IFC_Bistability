function [Pcorrect rBiasAmb PhaseDurAmb PhaseDurReplay UnclearAmb  UnclearReplay Exclusion] = prepare_data_pilot(observers, Exclusion)

%% prepare data from piloting experiment

%% input
% observers: participants
% Exlusion: struct for exclusion criteria

%% output
% Pcorrect: congruent perceptual transitions
% rBiasAmb: Bias to the right in ambiguous blocks
% PhaseDurAmb and PhaseDurReplay: Phae durations in ambiguity and parametric disambiguation
% UnclearAmb and UnclearReplay: Fraction of unclear perceptual decisions in ambiguity and parametric disambiguation
% Updated Exlusion-struct

%% Relevant Settings
n_runs = 4;
frames_per_rot = 1000;
rot_per_block = 5;
overlaps_per_rot = 8;
ifi = 1/85;
overlap_timing = [1:frames_per_rot/overlaps_per_rot:frames_per_rot*rot_per_block].*ifi-ifi;
n_conditions = 7;

for idx = 1:length(observers) % loop over participants

    for iidx = 1:n_runs % loop over runs
        
        clear Results
        load([observers(idx).name(1:22) '_run_' num2str(iidx) '.mat']) % load results
        
        %% Sort conditions according to blocks
        [B I] = sort(Results.disambiguation); % odd numbers in I: ambiguous blocks; even numbers in I: parametrically disambiguated blocks
        Pcorrect(:, idx,iidx) = Results.p_correct(I(mod(I,2) == 0)); % fraction of congruent perceptual decisions
        rBiasAmb(:, idx,iidx) = Results.r_bias(I(mod(I,2) == 1)); % fraction of perception of right-ward rotation in ambiguous blocks
        
        for trial = 1:length(Results.PDir) % loop over blocks
            
            % compute average phase durations
            Results.RelSwitchTime{trial} = Results.SwitchTime{trial}-Results.TrialStartTime{trial};
            
            if length(Results.RelSwitchTime{trial}) > 1 
                Results.PDur{trial} = diff(Results.RelSwitchTime{trial});
            else
                Results.PDur{trial} = Results.TrialEndTime{trial} - Results.TrialStartTime{trial};
            end
            Results.meanPhase(trial) = mean(Results.PDur{trial});
            
            % compute fraction of unclear perceptual responses
            Results.fraction_unclear_percept(trial) = length(find(Results.discrete_steps{trial} == 0)) / length(Results.discrete_steps{trial}); 
        end
        
        PhaseDurAmb(:, idx,iidx) = Results.meanPhase(I(mod(I,2) == 1)); % Phase Duration in ambiguous blocks
        PhaseDurReplay(:, idx,iidx) = Results.meanPhase(I(mod(I,2) == 0)); % Phase Duration in disambiguated blocks
        UnclearAmb(:, idx,iidx) = Results.fraction_unclear_percept(I(mod(I,2) == 1)); % Fraction unclear ambiguous blocks
        UnclearReplay(:, idx,iidx) = Results.fraction_unclear_percept(I(mod(I,2) == 0)); % Fraction unclear in disambiguated blocks
    end
    
end

        %% exlcude participants based on average phase duration or performance in fully disambiguated condition
        Exclusion.participants_to_exclude = unique([find(nanmean(nanmean(PhaseDurAmb,3),1) > Exclusion.Criteria.AveragePhase) ...
            find(mean(Pcorrect(n_conditions, :, :),3) < Exclusion.Criteria.Pcorrect)]);
        
        Pcorrect(:,Exclusion.participants_to_exclude,:) = [];
        rBiasAmb(:,Exclusion.participants_to_exclude,:) = [];
        PhaseDurAmb(:,Exclusion.participants_to_exclude,:) = [];
        PhaseDurReplay(:,Exclusion.participants_to_exclude,:) = [];
        UnclearAmb(:,Exclusion.participants_to_exclude,:) = [];
        UnclearReplay(:,Exclusion.participants_to_exclude,:) = [];
        
        % Exclusion of individual blocks because of no perceptual event
        ambiguous_blocks_to_exclude = find(PhaseDurAmb>Exclusion.Criteria.PhaseDuration);
        replay_blocks_to_exclude = find(PhaseDurReplay>Exclusion.Criteria.PhaseDuration);
       
        PhaseDurAmb(ambiguous_blocks_to_exclude) = NaN;
        UnclearAmb(ambiguous_blocks_to_exclude) = NaN;
        PhaseDurReplay(replay_blocks_to_exclude) = NaN;
        Pcorrect(replay_blocks_to_exclude) = NaN;
        rBiasAmb(replay_blocks_to_exclude) = NaN;
end
