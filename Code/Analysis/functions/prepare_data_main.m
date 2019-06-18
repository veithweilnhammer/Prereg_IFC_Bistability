function [Pcorrect rBias Frequency PhaseDur Unclear GLM Exclusion] = prepare_data_main(observers, Exclusion)

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
n_runs = 4;
frames_per_rot = 720;
rot_per_block = 10;
overlaps_per_rot = 8;
ifi = 1/60;
overlap_timing = [1:frames_per_rot/overlaps_per_rot:frames_per_rot*rot_per_block].*ifi-ifi;

for idx = 1:length(observers) % loop over participants
    
    for iidx = 1:n_runs % loop over runs
        
        clear Results
        load([observers(idx).name(1:28) '_run_' num2str(iidx) '.mat']) % load results
        
        for trial = 1:length(Results.PDir)
            
            % Sort conditions according to blocks
            Frequency(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx,iidx) = length(find(diff(Results.discrete_steps{trial}) ~= 0))/length(Results.discrete_steps{trial}); % Frequency of perceptual transition
            rBias(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx,iidx) = length(find(Results.discrete_steps{trial} == 1))/length(Results.discrete_steps{trial}); % rightward Bias
            PhaseDur(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx,iidx) = mean(diff(Results.discrete{trial})*ifi); % Phase Duration
            Unclear(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx,iidx) = length(find(Results.discrete_steps{trial} == -2))/length(Results.discrete_steps{trial}); % unclear
            
            if isnan(Results.transition_probability)
                Pcorrect(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx,iidx) = NaN;
            else
                Pcorrect(find(sort(Results.disambiguation) == Results.disambiguation(trial)),idx,iidx) = length(find(Results.template.discrete_steps{trial} == Results.discrete_steps{trial}))/length(Results.discrete_steps{trial});
            end
            
            GLM{idx}.run{iidx}.block{find(sort(Results.disambiguation) == Results.disambiguation(trial))}.T(1,:) = Results.discrete{trial}.*ifi+Results.TrialStartTime{trial}-Results.SessionStartTime;
            GLM{idx}.run{iidx}.block{find(sort(Results.disambiguation) == Results.disambiguation(trial))}.T(2,:) = zeros(1,length(Results.discrete{trial}.*ifi+Results.TrialStartTime{trial}-Results.SessionStartTime));
            
            if size(GLM{idx}.run{iidx}.block{find(sort(Results.disambiguation) == Results.disambiguation(trial))}.T,2) <= 1
                GLM{idx}.run{iidx}.block{find(sort(Results.disambiguation) == Results.disambiguation(trial))}.B(1,:) = Results.discrete{trial}(1).*ifi+Results.TrialStartTime{trial}-Results.SessionStartTime;
                GLM{idx}.run{iidx}.block{find(sort(Results.disambiguation) == Results.disambiguation(trial))}.B(2,:) = Results.TrialEndTime{trial}-Results.TrialStartTime{trial};
            else
                
                if ~isnan(Results.transition_probability)
                    clear C1 C2 congruency congruency_change end_of_events
                    congruency = Results.discrete_steps{trial}-Results.template.discrete_steps{trial};
                    congruency_change = find([1 diff(congruency)]~=0);
                    end_of_events = [[overlap_timing(congruency_change) frames_per_rot*rot_per_block*ifi] + Results.TrialStartTime{trial}-Results.SessionStartTime];
                    
                    if ~isempty(congruency_change(congruency(congruency_change) == 0))
                    C1(1,:) = overlap_timing(congruency_change(congruency(congruency_change) == 0)) + Results.TrialStartTime{trial}-Results.SessionStartTime;  % onset congruency
                    
                    % duration congruency
                    for events = 1:size(C1,2)
                        C1(2, events) = min(end_of_events(end_of_events > C1(1, events))) - C1(1, events);
                    end
                    else
                        C1 = [];
                    end
                    
                    if ~isempty(congruency_change(congruency(congruency_change) ~= 0))
                    C2(1,:) = overlap_timing(congruency_change(congruency(congruency_change) ~= 0)) + Results.TrialStartTime{trial}-Results.SessionStartTime;  % onset incongruency
                    
                    % duration incongruency
                    for events = 1:size(C2,2)
                        C2(2, events) = min(end_of_events(end_of_events > C2(1, events))) - C2(1, events);
                    end
                    
                    else
                    C2 = [];    
                    end
                    
                    GLM{idx}.run{iidx}.block{find(sort(Results.disambiguation) == Results.disambiguation(trial))}.C1 = C1;
                    GLM{idx}.run{iidx}.block{find(sort(Results.disambiguation) == Results.disambiguation(trial))}.C2 = C2;
                    
                end
                
                
            end
        end
        
    end
end




%% exlcude participants based on average phase duration or performance in fully disambiguated condition

% Exclusion of individual blocks because of no perceptual event
blocks_to_exclude = find(PhaseDur>Exclusion.Criteria.PhaseDuration);

PhaseDur(blocks_to_exclude) = NaN;
Unclear(blocks_to_exclude) = NaN;
Pcorrect(blocks_to_exclude) = NaN;
rBias(blocks_to_exclude) = NaN;


Exclusion.participants_to_exclude = unique([find(nanmean(nanmean(PhaseDur(:,:,:)),3) > Exclusion.Criteria.AveragePhase) find(nanmean(Pcorrect(6,:,2:4),3) < Exclusion.Criteria.Pcorrect)]);

