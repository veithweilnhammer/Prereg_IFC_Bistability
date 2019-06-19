
function [frequency correct] = get_conventional_control_data(Results)
for trial = 1:length(Results.PDir)
    frequency(trial,1) = length(find(diff(Results.discrete_steps{trial}) ~= 0))/length(Results.discrete_steps{trial});

if isnan(Results.transition_probability)
 correct(trial, 1) = NaN;
else
    correct(trial, 1) = length(find(Results.template.discrete_steps{trial} == Results.discrete_steps{trial}))/length(Results.discrete_steps{trial});
end
end
