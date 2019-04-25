
% Copyright (C) 2019 Veith Weilnhammer, CCM

clear all
close all
root_dir = '/home/veithweilnhammer/Public/Git/Prereg_IFC_BIstability/Code/Experimental_Scripts/Experimental_Scripts_Pilot/';

%% get oberserver
ObserverName= input('Please enter observer: ', 's');
Results.ObserverName = ObserverName;
session=input('Please enter session: ', 's');
Results.session=session;
if exist ([root_dir Results.ObserverName '\Results_' Results.ObserverName '_run_' num2str(session) '.mat'], 'file')
    error='Results file already exists. Please change ObserverName or Session'
    return;
end

%% run Settings
SettingsName= 'RDK'; %input('Please enter name of Settings_File: ', 's');
Settingsfile= [root_dir 'Settings/' 'Settings_' SettingsName '.m'];
run(Settingsfile);

%% save Settings
save([root_dir  'Settings/Settings_' Results.ObserverName '_run_' num2str(session) '.mat'])

%% get blue channel
BlueValue = input('Put the blue glass over the right eye. Please enter value of blue channel (red/green: 100): ', 's');

%% Initialize Psychtoolbox & Keyboard
AssertOpenGL;
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE'); left=KbName('LeftArrow'); right=KbName('RightArrow'); mixed=KbName('DownArrow');

RestrictKeysForKbCheck([escapeKey left right mixed]);

%% open screen & get monitor parameters
doublebuffer=1;
screens=Screen('Screens');
screenNumber=max(screens);

[w, rect] = Screen('OpenWindow', screenNumber, 0,[0 0 900 900], 32, doublebuffer+1, [], 128); % adjust antialiasing by last number

% Enable alpha blending with proper blend-function
% for drawing of smoothed points:
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

[center(1), center(2)] = RectCenter(rect);

fps=Screen('FrameRate',w);      % frames per second
ifi=Screen('GetFlipInterval', w); Results.ifi=ifi;
if fps==0
    fps=1/ifi;
end;

%% Colors
black = BlackIndex(w);
white = WhiteIndex(w);
red = [100 0 0];
blue = [0 0 str2num(BlueValue)];
green = [0 100 0];

%% ramp_colors
ramp_red = [[0:red(1)/n_frames_ramp:red(1)]' zeros(n_frames_ramp+1,2); ...
    repmat(red,frames_per_rot*rot_per_trial - 2 * (n_frames_ramp + 1),1); ...
    [red(1):-red(1)/n_frames_ramp:0]' zeros(n_frames_ramp+1,2)];

ramp_blue = [zeros(n_frames_ramp+1,2) [0:blue(3)/n_frames_ramp:blue(3)]'; ...
    repmat(blue,frames_per_rot*rot_per_trial - 2 * (n_frames_ramp + 1),1); ...
    zeros(n_frames_ramp+1,2) [blue(3):-blue(3)/n_frames_ramp:0]'];

ramp_green = [zeros(n_frames_ramp+1,1) [0:green(2)/n_frames_ramp:green(2)]' zeros(n_frames_ramp+1,1); ...
    repmat(green,frames_per_rot*rot_per_trial - 2 * (n_frames_ramp + 1),1); ...
    zeros(n_frames_ramp+1,1) [green(2):-green(2)/n_frames_ramp:0]' zeros(n_frames_ramp+1,1)];

%% Start presentation
HideCursor;	% Hide the mouse cursor
Priority(MaxPriority(w));

% Do initial flip...
vbl=Screen('Flip', w);

%%  prepare fixation and frame
ppd = pi * (rect(3)-rect(1)) / atan(mon_width/v_dist/2) / 360;    % pixels per degree
s = dot_w* ppd;   % dot size (pixels)

fix_cord = [center-fix_r*ppd center+fix_r*ppd];
FixCross = [center(1)-1,center(2)-fix_r*ppd,center(1)+1,center(2)+fix_r*ppd;center(1)-fix_r*ppd,center(2)-1,center(1)+fix_r*ppd,center(2)+1];

%% start session
Results.SessionStartTime=GetSecs;

for trial = 1:trial_n
    n_switch = 0;
    
    amb_points_to_plot = ones(1,length(r_x)); % determine which ambiguous dots to plot
    if ~mod(trial,2)
        permuted = randperm(length(r_x));
        amb_points_to_plot(permuted(1:round(disambiguation(trial)*length(permuted)))) = 0;
    end

    Results.PDir{trial} = [];
    Results.PosSwitch{trial} = [];
    Results.SwitchTime{trial} = [];
    
    % draw fixation dot
    Screen('FillRect', w, uint8(white), FixCross');
    
    Screen('Flip', w);
    WaitSecs(fixation_interval);
    Results.TrialStartTime{trial}=GetSecs;
    
    
    idx = 0; % index for angle of rotation
    
    for frames = 1:frames_per_rot*rot_per_trial
        
        % draw fixation dot
        Screen('FillRect', w, uint8(white), FixCross');
        
        if mod(trial,2)
            idx = idx+1;
            d_idx = idx; % index for disambiguated dots
        else
            idx = idx+Results.PDir{trial-1}(max(find(Results.discrete{trial-1} <= frames)));
            d_idx = idx + stereo_d; % index for disambiguated dots
        end
        
        % keep indices within bounds
        if idx > frames_per_rot
            idx = idx-frames_per_rot;
        end
        
        if d_idx > frames_per_rot
            d_idx = d_idx-frames_per_rot;
        end
        
        if idx <= 0
            idx = idx+frames_per_rot;
        end
        
        if d_idx <= 0
            d_idx = d_idx+frames_per_rot;
        end
        
        % draw Stimulus
        if mod(trial,2) % ambiguity
            Screen('DrawDots', w, [r_x(idx,:); r_z(idx,:)].*size_factor, s, ramp_blue(frames,:)+ramp_red(frames,:), center,dot_type);
        else 
            % parametric disambiguatiib
            if ~isempty(find(amb_points_to_plot == 1))
                Screen('DrawDots', w, [r_x(idx,amb_points_to_plot == 1); r_z(idx,amb_points_to_plot == 1)].*size_factor, s, ramp_blue(frames,:)+ramp_red(frames,:), center,dot_type);
            end
            Screen('DrawDots', w, [r_x(idx,amb_points_to_plot == 0); r_z(idx,amb_points_to_plot == 0)].*size_factor, s, ramp_blue(frames,:), center,dot_type);
            Screen('DrawDots', w, [r_x(d_idx,amb_points_to_plot == 0); r_z(d_idx,amb_points_to_plot == 0)].*size_factor, s, ramp_red(frames,:), center,dot_type);
        end
       
        Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
        
        %% response collection
        
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        if keyIsDown
            if keyCode(escapeKey)
                break; 
                close all
                
            elseif keyCode(right) | keyCode(left) | keyCode(mixed)   % right arrow
                n_switch=n_switch+1;
                Results.PDir{trial}(n_switch)= find(keyCode == 1); % +1: clockwise, -1: counterclockwise, 0: unclear
                Results.PosSwitch{trial}(n_switch)=frames;
                Results.SwitchTime{trial}(n_switch)=seconds;
            end
        end
        
        %% Flip
        if (doublebuffer==1)
            vbl=Screen('Flip', w, vbl + (waitframes-0.5)*ifi);
        end;
    end
    
    Results.TrialEndTime{trial}=seconds;
    
    %% if no response was given
    percepts = [-1 1];
    if isempty(Results.PDir{trial})
        Results.PDir{trial} = percepts(randi(2))
    end
    
    %% clear response entries
    
    RepDir=diff(Results.PDir{trial});
    Results.PDir{trial}(find(RepDir==0)+1)=[];
    Results.PosSwitch{trial}(find(RepDir==0)+1)=[];
    Results.SwitchTime{trial}(find(RepDir==0)+1)=[];
    
    %% rename directions
    
    Results.PDir{trial}(Results.PDir{trial} == right) = +1; %right
    Results.PDir{trial}(Results.PDir{trial} == left) = -1; %left
    Results.PDir{trial}(Results.PDir{trial} == mixed) = 0; %mixed
    
    %% round to overlaps
    Results.discrete{trial}(1) = 1;
    for idx = 2:length(Results.PosSwitch{trial})
        Results.discrete{trial}(idx) = max(overlaps(overlaps <= Results.PosSwitch{trial}(idx)));
    end
    Results.RT{trial} =  (Results.PosSwitch{trial} - Results.discrete{trial}).*ifi;
    
    Results.discrete_steps{trial} = zeros(1,length(overlaps));
    
    for idx = 1:length(Results.discrete_steps{trial})
        Results.discrete_steps{trial}(idx) = Results.PDir{trial}(max(find((Results.discrete{trial}<=overlaps(idx)))));
    end
    
    %% get parameters
    if mod(trial,2)
        Results.p_correct(trial) = 0.5;
        Results.r_bias(trial) = length(find(Results.discrete_steps{trial} == 1))/length(overlaps);
    else
        Results.p_correct(trial) = length(find(Results.discrete_steps{trial} == Results.discrete_steps{trial-1}))/length(overlaps);
        Results.r_bias(trial) = NaN;
    end
    
    %% print Results
    Results
    
    [ keyIsDown, seconds, keyCode ] = KbCheck;
        if keyIsDown
            if keyCode(escapeKey)
                break; 
                close all
            end
        end
end

    %% transfer settings
    Results.disambiguation = disambiguation;

%% Final fixation
Screen('FillRect', w, uint8(white), FixCross');
Screen('Flip', w);
WaitSecs(fixation_interval);
Results.SessionEndTime=GetSecs;

%% Save results
save([root_dir  'Results/Results_' Results.ObserverName '_run_' num2str(session) '.mat'])

%% Back to Matlab Screen
Priority(0);
ShowCursor
Screen('CloseAll');
