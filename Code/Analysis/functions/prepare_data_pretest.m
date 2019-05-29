function [Pcorrect rBias Frequency PhaseDur Unclear congruent_Rating, incongruent_Rating, Exclusion] = prepare_data_pretest(observers, Exclusion)

%% prepare data from piloting experiment

%% input
% observers: participants
% Exlusion: struct for exclusion criteria

%% output
% Pcorrect: congruent perceptual transitions
% rBias: Bias to the right in ambiguous blocks
% Frequency: probability of perceptual transitions at overlaps
% PhaseDur: average phase duration
% Unclear: frequency of unclear perceptual responses (relative to number of
% overlaps)
% Updated Exlusion-struct

%% Relevant Settings
n_runs = 3;
frames_per_rot = 1020;
rot_per_block = 10;
overlaps_per_rot = 8;
ifi = 1/85;
overlap_timing = [1:frames_per_rot/overlaps_per_rot:frames_per_rot*rot_per_block].*ifi-ifi;

for idx = 1:length(observers) % loop over participants
    
    for iidx = 1:n_runs % loop over runs
        
        clear Results
        load([observers(idx).name(1:22) '_run_' num2str(iidx) '.mat']) % load results
        
        for trial = 1:length(Results.PDir)
            
            if iidx < 3
                % Sort conditions according to blocks
                Frequency(trial,idx,iidx) = length(find(diff(Results.discrete_steps{trial}) ~= 0))/length(Results.discrete_steps{trial}); % Frequency of perceptual transition
                rBias(trial,idx,iidx) = length(find(Results.discrete_steps{trial} == 1))/length(Results.discrete_steps{trial}); % rightward Bias
                PhaseDur(trial,idx,iidx) = mean(diff(Results.discrete{trial})*ifi); % Phase Duration
                Unclear(trial,idx,iidx) = length(find(Results.discrete_steps{trial} == -2))/length(Results.discrete_steps{trial}); % unclear
                
                if isnan(Results.transition_probability)
                    Pcorrect(trial,idx,iidx) = NaN;
                else
                    Pcorrect(trial,idx,iidx) = length(find(Results.template.discrete_steps{trial} == Results.discrete_steps{trial}))/length(Results.discrete_steps{trial});
                end
                
            else
                % Sort conditions according to blocks
                Frequency(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx,iidx) = length(find(diff(Results.discrete_steps{trial}) ~= 0))/length(Results.discrete_steps{trial}); % Frequency of perceptual transition
                rBias(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx,iidx) = length(find(Results.discrete_steps{trial} == 1))/length(Results.discrete_steps{trial}); % rightward Bias
                PhaseDur(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx,iidx) = mean(diff(Results.discrete{trial})*ifi); % Phase Duration
                Unclear(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx,iidx) = length(find(Results.discrete_steps{trial} == -2))/length(Results.discrete_steps{trial}); % unclear
                
                congruent_Rating(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx) = Results.Congruent_Rating(trial);
                incongruent_Rating(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx) = Results.Incongruent_Rating(trial);
                
                if isnan(Results.transition_probability)
                    Pcorrect(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx,iidx) = NaN;
                else
                    Pcorrect(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx,iidx) = length(find(Results.template.discrete_steps{trial} == Results.discrete_steps{trial}))/length(Results.discrete_steps{trial});
                end
            end
            
        end
        
    end
end

%% Exclusion of individual blocks because of no perceptual event
blocks_to_exclude = unique([find(PhaseDur>Exclusion.Criteria.PhaseDuration) find(isnan(PhaseDur) == 1)]);

PhaseDur(blocks_to_exclude) = NaN;
Unclear(blocks_to_exclude) = NaN;
Pcorrect(blocks_to_exclude) = NaN;
rBias(blocks_to_exclude) = NaN;
Frequency(blocks_to_exclude) = NaN;

% ambiguous blocks exlcuded
for idx = 1:length(observers)
    Exclusion.amb_blocks_excluded(idx) = length(find(isnan(PhaseDur(:,idx,1))==1));
end

%% exlcude participants based on average phase duration or performance in fully disambiguated condition
Exclusion.participants_to_exclude = unique([find(mean(nanmean(PhaseDur(2:end,:,1)),3) > Exclusion.Criteria.AveragePhase) find(mean(Pcorrect(6,:,2),3) < Exclusion.Criteria.Pcorrect) find(Exclusion.amb_blocks_excluded >=3)]);



