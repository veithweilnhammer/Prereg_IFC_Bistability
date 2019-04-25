function [Pcorrect Tprobability mean_c_daq mean_ic_daq congruency daq pd] = simulate_data(experimental_duration, sampling_rate, IPS, DIS,Transition_Prob)

%% inputs
% experimental_duration: duration of block to sample 
% sampling_rate: rate of overlaps
% IPS: Initial precision of stability prior 
% DIS: precision of sensory evidence
% Transition_Prob: probabiltiy of change in sensory evidence at overlap

%% outputs:
% Pcorrect: congruent perceptual transitions
% Tprobability: probability of transition at overlap
% mean_c_daq: Mean Prediction errors during congruency between sensory evidence and perceptual decision
% mean_ic_daq: Mean Prediction errors during incongruency between sensory evidence and perceptual decision
% congruency: 1: incongruent perceptual decision; 0: congruent perceptual decision
% daq: PE timecourse
% pd: perceptual decision timecourse

percept2=1;  percept1=0;
outcomes = [percept2 percept1];

% modelling procedures equivalent to inversion

%% Initialize Inputs
u(:,1) = NaN(length(1:1/sampling_rate:experimental_duration),1); %% perceptual decisions
u(:,2) = repmat(DIS, length(1:1/sampling_rate:experimental_duration),1); %% precision of disambiguation

if DIS == 0
    u(:,3) = ones(length(1:1/sampling_rate:experimental_duration),1)-0.5; %% direction of disambiguation
else
    u(2,3) = randp([0.5 0.5],1,1) - 1;
    for k = 3:length(1:1/sampling_rate:experimental_duration)
        if randp([Transition_Prob 1-Transition_Prob],1,1) == 1; % introduce changes in the direction of disambiguation according to Transition_Prob
            u(k,3) = abs(u(k-1,3)-1);
        else
            u(k,3) = u(k-1,3);
        end
    end
    
end

for k = 1: length(1:1/sampling_rate:experimental_duration)
    
    %% Perceptual Stability
    
    if k==1
        Prediction(k,1)=0;
        Precision(k,1)=IPS; % initial precision of stability prior
        PerceptualAge(k,1) = 0;
    else
        Prediction(k,1)=u(k-1,1);
        
        if IPS~=0
            if PerceptualAge(k-1,1)==0
                Precision(k,1)=IPS; % reset precision of stability prior immediately after change in perceptual decision
            else
                Precision(k,1)=Precision(k-1,1)*exp(-abs(daq(k-1))); % reduce precision of the stability prior if no new perceptual deicision has been adopted
            end
        else
            Precision(k,1)=0;
        end
    end
    
    %% Disambiguation
    Precision(k,2)=DIS; % Precision of disambiguation
    Prediction(k,2)=u(k,3); % Direction of disambiguation
    
    %% Combine stability prior and disambiguation
    inform=find(Prediction(k,:)~=0.5);
    prior_precision(k,1)=sum(Precision(k, inform)); prior_sa(k,1)=1./prior_precision(k,1);
    prior_mean(k,1)=sum(Precision(k,inform).*Prediction(k,inform))/prior_precision(k,1);
    
    if prior_sa(k,1)==Inf
        ratio(k,1)=1;
    else
        ratio(k,1)=exp(((percept2-prior_mean(k,1)).^2-(percept1-prior_mean(k,1)).^2)...
            ./(2.*(prior_sa(k,1)).^2));
    end
    
    q(k,1)= 1/(ratio(k,1)+1); % joint posterior probability of right-ward rotatino
    
    %% Simulate response
    ind = randp([q(k,1) 1-q(k,1)],1,1);
    y(k,1) = outcomes(ind);
    u(k,1) = y(k,1);
    
    %% Compute prediction error
    daq(k,1) = u(k,1) - q(k,1);
    if k ~= 2
        w_daq(k,1) = daq(k)*(1/Precision(k,1));
    else
        w_daq(k,1) = 0;
    end
    
    %% Determine current duration of perceptual phase
    if k==1, NewPerceptualDecision(k,1)=1;
    else
        if u(k,1)==u(k-1,1), NewPerceptualDecision(k,1)=0;
            PerceptualAge(k,1)=k-max(find(NewPerceptualDecision(1:k-1,1)==1));
        else NewPerceptualDecision(k,1)=1;
            PerceptualAge(k,1) = 0;
        end
    end
    
    if k == 1
        PerceptualAge(k,1) == 0;
    else
        if u(k,1) == u(k-1,1)
            PerceptualAge(k,1)=k-max(find(NewPerceptualDecision(1:k-1,1)==1));
        else
            PerceptualAge(k,1) = 0;
        end
    end
end

%% Get measures for further analysis
% perceptual decisions
pd = u(:,1);

% Probability of congruent perceptual decisions
if u(:,2) == 0
Pcorrect = 0;
else    
Pcorrect =  length((find(u(:,1) == u(:,3))))/length(u);
end

% Trial by trial congruency
congruency = abs(u(:,1) - u(:,3)); % 0: congruent; 1: incongruent.

% mean_precision_error in (in-)congruency
mean_c_daq = NaN; mean_c_daq = nanmean(abs(daq(congruency == 0))); % average prediction error at congruent overlap
mean_ic_daq = NaN; mean_ic_daq = nanmean(abs(daq(congruency == 1))); %average prediction error at incongruent overlap

% Transisition probability
Tprobability = length(find(diff(u(:,1)~= 0)))/length(u);
if isempty(Tprobability)
    Tprobability = 0;
end
