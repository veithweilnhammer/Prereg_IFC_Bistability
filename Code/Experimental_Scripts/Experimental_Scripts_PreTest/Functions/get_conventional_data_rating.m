
function [frequency correct rating] = get_conventional_data_rating(Results)
for trial = 1:length(Results.PDir)

    % transition probability per overlap
    frequency(trial,1) = length(find(diff(Results.discrete_steps{trial}) ~= 0))/length(Results.discrete_steps{trial});

    
    
    
if isnan(Results.transition_probability)
 correct(trial, 1) = NaN;
else

  
    % fraction of perceptual decisions congruent to current sensory
    % evidence; 1st row: lowest level of disambigation; last row: highest
    % level of disambiguation
    correct(find(sort(Results.disambiguation) == Results.disambiguation(trial)), 1) = length(find(Results.template.discrete_steps{trial} == Results.discrete_steps{trial}))/length(Results.discrete_steps{trial});
    rating(find(sort(Results.disambiguation) == Results.disambiguation(trial)), 1) = Results.Congruent_Rating(trial);
    rating(find(sort(Results.disambiguation) == Results.disambiguation(trial)), 2) = Results.Incongruent_Rating(trial);
end
end
