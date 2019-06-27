clear all
close all
% instruction: adjust color by minimizing flickering (using arrow up and down keys), press arrow right key to go to next trial

%  def
screens=Screen('Screens');
scrnNum=max(screens);
present_window=[0 0 900 900];
n_trials=10;
blue_value=200; %pp_managemat(2,5);
red_channel=1; %pp_managemat(2,6);
cond_per_trial=round(rand(1,n_trials))+1;
blue_val_per_trial=repmat(blue_value,1,n_trials);
hz_per_trial=repmat(15,1,n_trials);
imSize=244; %round(bin_imSize*0.8); % 244
axis_diam=8;
frames_per_color=NaN;

%% stereoscope
stereo = 0;

%% key def

KbName('UnifyKeyNames');

% keyboard
escapeKey = KbName('ESCAPE');
left=KbName('LeftArrow');
%right=KbName('RightArrow');
%up = KbName('UpArrow');
%down=KbName('DownArrow');

% scanner
 down=KbName('1!'); 
 up=KbName('2@');
 right=KbName('3#'); 

keyswanted=[right,escapeKey,up,down];
management_buttons=[right,escapeKey];
Screen('Preference', 'SkipSyncTests', 1);

% engine
[windowPtr, windowRect]=Screen('OpenWindow', scrnNum, BlackIndex(scrnNum), present_window);


if stereo == 1
    stim_center{1}=[diff(windowRect([1,3]))*0.25 diff(windowRect([2,4]))*0.5];
    stim_center{2}=[diff(windowRect([1,3]))*0.75 diff(windowRect([2,4]))*0.5];
else
    stim_center{1} = [diff(windowRect([1,3]))*0.5 diff(windowRect([2,4]))*0.5];
    stim_center{2}=[diff(windowRect([1,3]))*0.5 diff(windowRect([2,4]))*0.5];
end

refresh_rate=Screen('NominalFrameRate', scrnNum);

[color_m square_m1 square_m2] = makeframes(round(imSize*1.1),stim_center,windowRect(3),windowRect(4));

white = WhiteIndex(windowPtr); % pixel value for white
black = BlackIndex(windowPtr); % pixel value for black
gray = (white+black)/2;
inc = white-black;
final_adjusted_color=nan(n_trials,3);
regist_trials=[];
fprintf('\nColor-Matching Flickering Task (n=%d trials; PP should take about 10-20 seconds per trial):\n',n_trials);
tic
for i=1:n_trials
    fprintf(' Trial #%d (%d remaining)',i,n_trials-i);
    if i>1
        fprintf(' [Duration of last trial was %0.1f seconds, adjusted red value was %0.2f]\n',toc,final_adjusted_color(i-1,red_channel));
        tic
    else
        fprintf('\n');
    end
    ref_color=[0 0 blue_val_per_trial(i)];
    frames_per_color=round(refresh_rate/hz_per_trial(i)/2);
    t0=Screen('Flip', windowPtr);
    adjusted_color=zeros(1,3); % init
    %random_start_value=floor(rand(1)*256); % start with a random value between 0 and 255
    random_start_value=min([255 floor(rand(1)*1*blue_val_per_trial(i))]); % start with a random value between 0 and twice the blue value (as long as it is not greater than 255)
    adjusted_color(red_channel)=random_start_value;
    while KbCheck, KbCheck; end; pause(0.1);
    regist_frames=[];
    break_trial=0;
    while break_trial==0
        
        for f=1:frames_per_color
            DrawFormattedText(windowPtr, sprintf('%d',n_trials-i+1), stim_center{2}(1)+imSize/2, stim_center{2}(2)+imSize/2, [255 255 255]);
            Screen('FillOval', windowPtr, ref_color ,[stim_center{1}-imSize/2 stim_center{1}+imSize/2])%Screen('FillOval', windowPtr, ref_color ,[windowRect(3:4)./2-imSize/2 windowRect(3:4)./2+imSize/2])
            Screen('FrameRect', windowPtr, color_m , square_m1, 2);
            Screen('FillOval', windowPtr, ref_color ,[stim_center{2}-imSize/2 stim_center{2}+imSize/2])%Screen('FillOval', windowPtr, ref_color ,[windowRect(3:4)./2-imSize/2 windowRect(3:4)./2+imSize/2])
            Screen('FrameRect', windowPtr, color_m , square_m2, 2);
            t1=Screen('Flip', windowPtr);
            regist_frames(end+1,:)=[i t1-t0 hz_per_trial(i) 1 ref_color];
        end
        
        for f=1:frames_per_color
            DrawFormattedText(windowPtr, sprintf('%d',n_trials-i+1), stim_center{2}(1)+imSize/2, stim_center{2}(2)+imSize/2, [255 255 255]);
            Screen('FillOval', windowPtr, adjusted_color ,[stim_center{1}-imSize/2 stim_center{1}+imSize/2])%[windowRect(3:4)./2-imSize/2 windowRect(3:4)./2+imSize/2])
            Screen('FrameRect', windowPtr, color_m , square_m1, 2);
            Screen('FillOval', windowPtr, adjusted_color ,[stim_center{2}-imSize/2 stim_center{2}+imSize/2])%[windowRect(3:4)./2-imSize/2 windowRect(3:4)./2+imSize/2])
            Screen('FrameRect', windowPtr, color_m , square_m2, 2);
            t1=Screen('Flip', windowPtr);
            regist_frames(end+1,:)=[i t1-t0 hz_per_trial(i) 2 adjusted_color];
            [pressed, secs, keyCode] = KbCheck;
            if pressed==1
                if keyCode(keyswanted(3)) &&  adjusted_color(red_channel)<255 % up
                    adjusted_color(red_channel)=adjusted_color(red_channel)+1;
                elseif keyCode(keyswanted(4)) && adjusted_color(red_channel)>0 % down
                    adjusted_color(red_channel)=adjusted_color(red_channel)-1;
                elseif keyCode(management_buttons(1)) % designated continue button
                    break_trial=1;
                    break
                elseif keyCode(management_buttons(2)) % ESC
                    keep_running=0;
                    break_trial=1;
                    break
                end
            end
        end
    end
    regist_trials=[regist_trials; regist_frames];
    final_adjusted_color(i,:)=adjusted_color;
end
%% check results
Screen('FrameRect', windowPtr, [0 0 0] , square_m2, 2);
t1=Screen('Flip', windowPtr);
fprintf('\nReview of results of the color matching task:\n')
fprintf('values of the red channel: M=%0.2f SD=%0.2f Median=%0.2f Min=%0.2f Max=%0.2f\n',mean(final_adjusted_color(:,red_channel)),std(final_adjusted_color(:,red_channel)),median(final_adjusted_color(:,red_channel)),min(final_adjusted_color(:,red_channel)),max(final_adjusted_color(:,red_channel)));

%% terminate
Screen('CloseAll')