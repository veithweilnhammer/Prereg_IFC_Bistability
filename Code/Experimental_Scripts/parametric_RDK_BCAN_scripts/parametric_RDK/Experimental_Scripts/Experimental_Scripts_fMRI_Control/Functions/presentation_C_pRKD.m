
function [Results] = presentation_C_pRKD(root_dir, SettingsName, ObserverName, BlueValue, session, balance, transition_probability)


%% get oberserver

Results.ObserverName = ObserverName;
Results.session=session;

%% Get balance between left/right
balance = balance(randperm(length(balance)));
Results.balance = balance; % temporal balance across conditions in randomized order
trial_n = length(Results.balance); Results.trial_n = trial_n;

Results.transition_probability = transition_probability; % probability of a change of planar motion direction at timepoints corresponding to overlaps in the main exerperiment

%% run Settings
Settingsfile= [root_dir 'Settings/' 'Settings_' SettingsName '.m'];
run(Settingsfile);

%% Create perceptual timecourse
percepts = [-1 1];
for trial = 1:trial_n
    
    a = repmat(percepts(randi(2)), 1, round(length(overlaps)*balance(trial)));
    
    if length(a) == length(overlaps) % full imbalance
        template.discrete_steps{trial} = a;
    else
        
        b = repmat(a(1)*(-1), 1, length(overlaps)-round(length(overlaps)*balance(trial)));
        
        %% find optimal combination of imbalance and transition probability
        
        for counter = 1:1000 
           
            xx(1, counter) = percepts(randi(2));
            
            for idx = 2:length(overlaps)
             if randp([transition_probability 1-transition_probability],1,1) == 1;
                xx(idx, counter) = (-1)*xx(idx-1, counter);
            else
                xx(idx, counter) = xx(idx-1, counter);
            end
            end
        error(counter) = abs(length(find(xx(:,counter) == mode(xx(:,counter))))-length(a));    
        end
        
        point_to=find(error==min(error));  template.discrete_steps{trial} = xx(:,point_to(1))';
        %%
    end
    
    template.discrete{trial} = overlaps(find([1 diff(template.discrete_steps{trial})] ~= 0));
    template.PDir{trial} = template.discrete_steps{trial}(find([1 diff(template.discrete_steps{trial})] ~= 0));
    template.balance(trial) = length(find(template.discrete_steps{trial} == mode(template.discrete_steps{trial})))/length(template.discrete_steps{trial});
end
accuracy = Results.balance - template.balance;
Results.template = template; 
Results.accuracy = accuracy;


%% Initialize Psychtoolbox & Keyboard

AssertOpenGL;
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE'); 

%% KbNames PC
 left=KbName('1!'); 
 down=KbName('2@');
 right=KbName('3#'); 


%% KbNames Scannner
% left=KbName('1!'); right=KbName('2@');
% down=KbName('3#'); 

%%
trigger=KbName('5%');
%%
RestrictKeysForKbCheck([escapeKey left right down trigger]);

%% open screen & get monitor parameters

doublebuffer=1;
screens=Screen('Screens');
screenNumber=max(screens);

%[w, rect] = Screen('OpenWindow', screenNumber, 0,[1920 0 2*1920 1080], 32, doublebuffer+1, [], 128); % adjust antialiasing by last number
[w, rect] = Screen('OpenWindow', screenNumber, 0,[], 32, doublebuffer+1, [], 128); % adjust antialiasing by last number
Results.Monitor.rect = rect;


% Enable alpha blending with proper blend-function
% for drawing of smoothed points:
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

[center(1), center(2)] = RectCenter(rect);

fps=Screen('FrameRate',w);      % frames per second
ifi=Screen('GetFlipInterval', w); Results.Monitor.ifi=ifi;
if fps==0
    fps=1/ifi;
end;

%% Colors

black = BlackIndex(w);
white = WhiteIndex(w);
red = [str2num(RedValue) 0 0];
blue = [0 0 str2num(BlueValue)];
green = [0 100 0];

%% ramp_colors: ramp in/out and beginning/end of block

ramp_red = [[0:red(1)/n_frames_ramp:red(1)]' zeros(n_frames_ramp+1,2); ...
    repmat(red,frames_per_cycle*rot_per_trial - 2 * (n_frames_ramp + 1),1); ...
    [red(1):-red(1)/n_frames_ramp:0]' zeros(n_frames_ramp+1,2)];

ramp_blue = [zeros(n_frames_ramp+1,2) [0:blue(3)/n_frames_ramp:blue(3)]'; ...
    repmat(blue,frames_per_cycle*rot_per_trial - 2 * (n_frames_ramp + 1),1); ...
    zeros(n_frames_ramp+1,2) [blue(3):-blue(3)/n_frames_ramp:0]'];

ramp_green = [zeros(n_frames_ramp+1,1) [0:green(2)/n_frames_ramp:green(2)]' zeros(n_frames_ramp+1,1); ...
    repmat(green,frames_per_cycle*rot_per_trial - 2 * (n_frames_ramp + 1),1); ...
    zeros(n_frames_ramp+1,1) [green(2):-green(2)/n_frames_ramp:0]' zeros(n_frames_ramp+1,1)];

%% Start presentation
HideCursor;	% Hide the mouse cursor
Priority(MaxPriority(w));

% Do initial flip...
vbl=Screen('Flip', w);

%%  prepare fixation and frame

ppd = pi * (rect(3)-rect(1)) / atan(mon_width/v_dist/2) / 360; % pixels per degree
Results.Monitor.ppd = ppd;

s = dot_w* ppd;   % dot size (pixels)
s_fix_r = fix_r*ppd; % fixation dot size (pixels)


%% wait for scanner trigger
while 1
    [ keyIsDown, seconds, keyCode ] = KbCheck;
    if keyIsDown
        if keyCode(trigger)
            break;
        end
    end
end

%% start session
Results.SessionStartTime=GetSecs;

for trial = 1:trial_n
    n_switch = 0;
    
    % Initialize Results structures
    Results.PDir{trial} = [];
    Results.PosSwitch{trial} = [];
    Results.SwitchTime{trial} = [];
    
    % draw fixation dot
    Screen('DrawDots', w, [0; 0], s_fix_r, ramp_blue(100,:)+ramp_red(100,:), center,dot_type);
    
    Screen('Flip', w);
    WaitSecs(fixation_interval);
    Results.TrialStartTime{trial}=GetSecs;
    
    
    for frames = 1:frames_per_cycle*rot_per_trial
        
        % draw fixation dot
        Screen('DrawDots', w, [0; 0], s_fix_r, ramp_blue(frames,:)+ramp_red(frames,:), center,dot_type);
        
        % move dots to the left or to the right
        x = x + shift_per_frame*template.PDir{trial}(max(find(template.discrete{trial} <= frames)));
        
        % make sure that dots are within bounds
        outside = find(((x).^2 + (z).^2) >= 1);
        
        %relocate dots
        x(outside) = -x(outside) + shift_per_frame*template.PDir{trial}(max(find(template.discrete{trial} <= frames)));
        
        % draw Stimulus
        
        Screen('DrawDots', w, [x; z].*size_factor, s, ramp_blue(frames,:)+ramp_red(frames,:), center,dot_type);
        Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
        
        %% response collection
        
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        if keyIsDown
            
            if length(find(keyCode==1)) == 1 % against multiple events
                if keyCode(escapeKey)
                    break; close all
                    
                elseif keyCode(right) | keyCode(left) | keyCode(down)
                    n_switch=n_switch+1;
                    Results.PDir{trial}(n_switch)= find(keyCode == 1);
                    Results.PosSwitch{trial}(n_switch)=frames;
                    Results.SwitchTime{trial}(n_switch)=seconds;
                end
            end
            
        end
        
        %% Flip
        if (doublebuffer==1)
            vbl=Screen('Flip', w, vbl + (waitframes-0.5)*ifi);
        end;
    end
    
    Results.TrialEndTime{trial}=seconds;
    
    %% if no response was given
    if isempty(Results.PDir{trial})
        Results.PDir{trial} = percepts(randi(2));
    end
    
    %% clear response entries
    
    RepDir=diff(Results.PDir{trial});
    Results.PDir{trial}(find(RepDir==0)+1)=[];
    Results.PosSwitch{trial}(find(RepDir==0)+1)=[];
    Results.SwitchTime{trial}(find(RepDir==0)+1)=[];
    
    %% rename directions
    
    Results.PDir{trial}(Results.PDir{trial} == right) = +1; %right
    Results.PDir{trial}(Results.PDir{trial} == left) = -1; %left
    Results.PDir{trial}(Results.PDir{trial} == down) = -2; %mixed
    
    %% round to overlaps
    
    Results.discrete{trial}(1) = 1;
    for idx = 2:length(Results.PosSwitch{trial})
        Results.discrete{trial}(idx) = max(overlaps(overlaps <= Results.PosSwitch{trial}(idx)-minimal_response_frames));
    end
    Results.RT{trial} =  (Results.PosSwitch{trial} - Results.discrete{trial})*ifi;
    
    Results.discrete_steps{trial} = zeros(1,length(overlaps));
    
    for idx = 1:length(Results.discrete_steps{trial})
        Results.discrete_steps{trial}(idx) = Results.PDir{trial}(max(find((Results.discrete{trial}<=overlaps(idx)))));
    end
    
    %% break out of experiment
    
    [ keyIsDown, seconds, keyCode ] = KbCheck;
    if keyIsDown
        if keyCode(escapeKey)
            break; close all
        end 
    end  
end

%% Final fixation

Screen('DrawDots', w, [0; 0], s_fix_r, ramp_blue(100,:)+ramp_red(100,:), center,dot_type);
Screen('Flip', w);
WaitSecs(fixation_interval);
Results.SessionEndTime=GetSecs;

%% Save results
save([root_dir  'Results/Results_' Results.ObserverName '_' num2str(session) '.mat'], 'Results')

%% Back to Matlab Screen
Priority(0);
ShowCursor
Screen('CloseAll');

end

