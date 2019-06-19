function [colour_m square_m1 square_m2] = makeframes(SettingsStimuliSize,SettingsStimuliCenterXY,x,y)

%% square matrix
    
    %function rect = recter(siz, coord, Rect, shift)
    Settings.squares.size = SettingsStimuliSize/23; %SettingsStimuliSize/13;
    estimate_length = round(-((SettingsStimuliSize+2*Settings.squares.size)/2):1.5*Settings.squares.size:((SettingsStimuliSize+2*Settings.squares.size)/2));
    le=length(estimate_length);
    square_coord = zeros(2,4*le);
    square_coord(1,1:le) = round(-((SettingsStimuliSize+2*Settings.squares.size)/2):1.5*Settings.squares.size:((SettingsStimuliSize+2*Settings.squares.size)/2));
    square_coord(2,1:le) = round(repmat(-((SettingsStimuliSize+2*Settings.squares.size)/2), 1, le));
    square_coord(1,le+1:2*le)= round(-((SettingsStimuliSize+2*Settings.squares.size)/2):1.5*Settings.squares.size:((SettingsStimuliSize+2*Settings.squares.size)/2));
    square_coord(2,le+1:2*le) = round(repmat(((SettingsStimuliSize+2*Settings.squares.size)/2), 1, le));
    square_coord(1,2*le+1:3*le) = round(repmat(-((SettingsStimuliSize+2*Settings.squares.size)/2), 1, le));
    square_coord(2,2*le+1:3*le) = round(-((SettingsStimuliSize+2*Settings.squares.size)/2):1.5*Settings.squares.size:((SettingsStimuliSize+2*Settings.squares.size)/2));
    square_coord(1,3*le+1:4*le) = round(repmat(((SettingsStimuliSize+2*Settings.squares.size)/2), 1, le));
    square_coord(2,3*le+1:4*le) = round(-((SettingsStimuliSize+2*Settings.squares.size)/2):1.5*Settings.squares.size:((SettingsStimuliSize+2*Settings.squares.size)/2));
    
    square_m1=zeros(4,length(square_coord(1,:)));
    square_m2=zeros(4,length(square_coord(1,:)));
    
    for i= 1:length(square_coord(1,:))
        square = recter(Settings.squares.size, SettingsStimuliCenterXY{1}, [x,y], [square_coord(1,i),square_coord(2,i)]);
        square = square';
        square_m1(:,i)=square;
    end
    
    for i= 1:length(square_coord(1,:))
        square = recter(Settings.squares.size, SettingsStimuliCenterXY{2}, [x,y], [square_coord(1,i),square_coord(2,i)]);
        square = square';
        square_m2(:,i)=square;
    end
    
    leave_out=0;
    choose=randperm(length(square_coord(1,:)));
    choose=choose(1:leave_out);
    square_m1(:,choose)=[];
    square_m2(:,choose)=[];
    
    %% colour matrix
    
    v=zeros(1,length(square_coord(1,:))-leave_out);
    for i = 1:length(square_coord(1,:))-leave_out
        v(i)= (0.3+round(rand)*0.7)*256;
    end
    
    colour_m = [v;v;v];