
function [Results] = presentation_pRKD_certainty(root_dir, SettingsName, ObserverName, BlueValue, session, transition_probability)


%% get oberserver

Results.ObserverName = ObserverName;
Results.session=session;

%% run Settings

Settingsfile= [root_dir 'Settings/' 'Settings_' SettingsName '.m'];
run(Settingsfile);

%% Get direction of disambiguation

Results.transition_probability = transition_probability;

percepts = [-1 1];
for trial = 1:trial_n
    template.discrete_steps{trial} = zeros(1,length(overlaps));
    
    if isnan(transition_probability)
        
        template.discrete_steps{trial} = repmat(0.5,1,length(overlaps)); % transition probability not defined for ambiguous presentation
        
    else
        template.discrete_steps{trial}(1) = percepts(randi(2)); % choose initial disambiguation at random
        
        for idx = 2 : length(overlaps) % switch direction of disambiguation according to transition probability
            if randp([transition_probability 1-transition_probability],1,1) == 1;
                template.discrete_steps{trial}(idx) = (-1)*template.discrete_steps{trial}(idx-1);
            else
                template.discrete_steps{trial}(idx) = template.discrete_steps{trial}(idx-1);
            end
        end
    end
    
    template.discrete{trial} = overlaps(find([1 diff(template.discrete_steps{trial})] ~= 0)); % overlaps for transitions
    template.PDir{trial} = template.discrete_steps{trial}(find([1 diff(template.discrete_steps{trial})] ~= 0)); % directions
end

Results.template = template;

%% when to check perception
frames_to_check = overlaps(3:length(overlaps)-1) + round(mean(diff(overlaps))*0.8);
%% Initialize Psychtoolbox & Keyboard

AssertOpenGL;
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');

%% KbNames PC
left=KbName('LeftArrow'); right=KbName('RightArrow');
down=KbName('DownArrow');
one = KbName('1!'); two = KbName('2@');
three = KbName('3#'); four = KbName('4$');

%% KbNames Scannner
% left = KbName('1!'); right = KbName('2@');
% down = KbName('3#'); four = KbName('4$');

%%
trigger=KbName('5%');
%%
RestrictKeysForKbCheck([escapeKey left right down trigger one two three four]);

%% open screen & get monitor parameters

doublebuffer=1;
screens=Screen('Screens');
screenNumber=max(screens);

[w, rect] = Screen('OpenWindow', screenNumber, 0,[0 0 900 900], 32, doublebuffer+1, [], 128); % adjust antialiasing by last number
%[w, rect] = Screen('OpenWindow', screenNumber, 0,[], 32, doublebuffer+1, [], 128); % adjust antialiasing by last number
Results.Monitor.rect = rect;


% Enable alpha blending with proper blend-function
% for drawing of smoothed points:
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

[center(1), center(2)] = RectCenter(rect);

fps=Screen('FrameRate',w);      % frames per second
ifi=Screen('GetFlipInterval', w); Results.Monitor.ifi=ifi; % inter-frame interval
if fps==0
    fps=1/ifi;
end;

%% Colors

black = BlackIndex(w);
white = WhiteIndex(w);
red = [100 0 0];
blue = [0 0 str2num(BlueValue)];
green = [0 100 0];

%% ramp_colors: colors of dots across frames (ramp stimulus in/out at beginning/end of block)

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

ppd = pi * (rect(3)-rect(1)) / atan(mon_width/v_dist/2) / 360; % pixels per degree
Results.Monitor.ppd = ppd;

s = dot_w* ppd;   % dot size (pixels)
s_fix_r = fix_r*ppd; % fixation dot size (pixels)


%% wait for scanner trigger
%while 1
%   [ keyIsDown, seconds, keyCode ] = KbCheck;
%  if keyIsDown
%     if keyCode(trigger)
%        break;
%   end
%end
%end

%% start session
Results.SessionStartTime=GetSecs;


Results.Congruent_Rating = NaN;
Results.Incongruent_Rating = NaN;

for trial = 1:trial_n
    
    
    Results.Congruent_Rating(trial) = NaN;
    Results.Incongruent_Rating(trial) = NaN;
    n_switch = 0;
    
    % Determine which points to disambiguate
    amb_points_to_plot = ones(1,length(r_x));
    if ~isnan(transition_probability)
        permuted = randperm(length(r_x));
        amb_points_to_plot(permuted(1:round(disambiguation(trial)*length(permuted)))) = 0;
    end
    
    % Initialize Results structures
    Results.PDir{trial} = [];
    Results.PosSwitch{trial} = [];
    Results.SwitchTime{trial} = [];
    
    % draw fixation dot
    Screen('DrawDots', w, [0; 0], s_fix_r, ramp_blue(100,:)+ramp_red(100,:), center,dot_type);
    
    Screen('Flip', w);
    WaitSecs(fixation_interval);
    Results.TrialStartTime{trial}=GetSecs;
    
    idx = 0; % index for angle of rotation
    
    for frames = 1:frames_per_rot*rot_per_trial
        
        % draw fixation dot
        Screen('DrawDots', w, [0; 0], s_fix_r, ramp_blue(frames,:)+ramp_red(frames,:), center,dot_type);
        
        if isnan(transition_probability) % no disambiguation for ambiguous stimulation (indeces are identical)
            idx = idx+1;
            d_idx = idx;
        else
            idx = idx+template.PDir{trial}(max(find(template.discrete{trial} <= frames)));
            d_idx = idx + stereo_d; % index for angle of rotation for disambiguated dots
        end
        
        % Make sure that indeces are whithin bounds
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
        if isnan(transition_probability)
            Screen('DrawDots', w, [r_x(idx,:); r_z(idx,:)].*size_factor, s, ramp_blue(frames,:)+ramp_red(frames,:), center,dot_type);
        else
            
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
        
        if any(frames == frames_to_check)
            
            if ~isempty(Results.PDir{trial})
                
               
                    Results.PDir{trial}(Results.PDir{trial} == right) = +1; %right
                    Results.PDir{trial}(Results.PDir{trial} == left) = -1; %left
                    Results.PDir{trial}(Results.PDir{trial} == down) = -2; %mixed
                 if abs(Results.PDir{trial}(length(Results.PDir{trial}))) ~= 2   
                    if Results.template.discrete_steps{trial}(max(find(overlaps < frames))) == Results.PDir{trial}(length(Results.PDir{trial})) && isnan(Results.Congruent_Rating(trial));
                        
                        Screen('DrawDots', w, [0; 0], s_fix_r, ramp_blue(frames,:)+ramp_red(frames,:), center,dot_type);
                        [nx, ny, textbounds] = DrawFormattedText(w, ['Perceptual Rating for Trial ' num2str(trial) '\n\n \n\n \n\nClear: Press [1] \n\nRather clear: Press [2] \n\nRather unclear: Press [3] \n\nUnclear: Press [4]'], ...
                            center(1)-rect(3)/8, center(2)-rect(4)/6 , ramp_blue(frames,:)+ramp_red(frames,:));
                        Screen('Flip',w);
                        
                        while 1
                            [ keyIsDown, seconds, keyCode ] = KbCheck;
                            if keyIsDown
                                if length(find(keyCode==1)) == 1 % against multiple events
                                    if keyCode(one) | keyCode(two) | keyCode(three) | keyCode(four)
                                        Results.Congruent_Rating(trial) = find(keyCode == 1) - 10;
                                        break;
                                    end
                                end
                            end
                        end
                        
                    elseif Results.template.discrete_steps{trial}(max(find(overlaps < frames))) ~= Results.PDir{trial}(length(Results.PDir{trial})) && isnan(Results.Incongruent_Rating(trial));
                        
                        Screen('DrawDots', w, [0; 0], s_fix_r, ramp_blue(frames,:)+ramp_red(frames,:), center,dot_type);
                        [nx, ny, textbounds] = DrawFormattedText(w, ['Perceptual Rating for Trial ' num2str(trial) '\n\n \n\n \n\nClear: Press [1] \n\nRather clear: Press [2] \n\nRather unclear: Press [3] \n\nUnclear: Press [4]'], ...
                            center(1)-rect(3)/8, center(2)-rect(4)/6 , ramp_blue(frames,:)+ramp_red(frames,:));
                        Screen('Flip',w);
                        
                        while 1
                            [ keyIsDown, seconds, keyCode ] = KbCheck;
                            if keyIsDown
                                if length(find(keyCode==1)) == 1 % against multiple events
                                    if keyCode(one) | keyCode(two) | keyCode(three) | keyCode(four)
                                        Results.Incongruent_Rating(trial) = find(keyCode == 1) - 10;
                                        break;
                                    end
                                end
                            end
                        end
                        
                    end
                    
                end
            end
        end
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

