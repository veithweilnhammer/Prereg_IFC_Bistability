%%
%% Conditions

%disambiguation = [0.15 0.30 0.45 0.60 0.75 1]; % levels of disambiguation
disambiguation = [0.95 0.96 0.97 0.98 0.99 1]; % levels of disambiguation

disambiguation = disambiguation(randperm(length(disambiguation))); % randomize order of conditions.

Results.disambiguation = disambiguation;

%% Experimental Duration

fixation_interval = 10; Results.fixation_interval = fixation_interval;

trial_n = length(Results.disambiguation); Results.trial_n = trial_n; % each level of disambiguation is presented once per run.

rot_per_trial = 10; Results.rot_per_trial = rot_per_trial; % number of rotations of the RDK per block.

desired_time_per_rotation = 12; % desired time per full rotation of the RDK

estimated_ifi = 1/60; % 0.0167 % expected inter-frame interval.

frames_per_rot = round(desired_time_per_rotation/estimated_ifi); Results.Stimulus.frames_per_rot = frames_per_rot; % number of frames per rotation

%% compute predicted duration of each run

estimated_experimental_duration = (estimated_ifi*frames_per_rot*rot_per_trial*trial_n + fixation_interval*(trial_n + 1))/60;

%% initialize sphere

max_number_of_dots = 5000; % number of dots in a full sphere
Results.Stimulus.max_numer_of_dots = max_number_of_dots;
step_size_dots = 2*pi/max_number_of_dots; 
r = 1; % radius of sphere
X_c = 0; Y_c = 0; Z_c = 0; % center of sphere

theta = [0:step_size_dots:2*pi-step_size_dots];
phi = pi.*rand(1, length(theta)) - 0.5*pi;

x = X_c + r.*cos(phi).*cos(theta); % initial dot configuration on sphere (x axis)
y = Y_c + r.*cos(phi).*sin(theta); % (y axis)
z = Z_c + r.*sin(phi); % (z axis)

%%  cut out slices
fraction_of_surface = 0.3; Results.Stimulus.fraction_of_surface = fraction_of_surface; % fraction of dots to keep.

cutout = find(abs(x) > fraction_of_surface & abs(y) > fraction_of_surface); %% kings apple
z(cutout) = [];
y(cutout) = [];
x(cutout) = [];
phi(cutout) = [];
theta(cutout) = [];

%% overlapping configurations

n_overlaps = 8; Results.n_overlaps = n_overlaps; % 8 overlaps on full rotation of the "kings apple"
lag_between_overlaps = frames_per_rot / n_overlaps;
overlaps = [1 : lag_between_overlaps : frames_per_rot*rot_per_trial]; Results.overlaps = overlaps; % frames of overlaps

minimal_response_time = 0.2; Results.minimal_response_time = minimal_response_time; % assume minimal time of respose after transition at an overlap
minimal_response_frames = round(minimal_response_time/estimated_ifi);

%% rotate sphere: computing dot positions around a full rotation of the sphere.

shift_per_frame = 2*pi/frames_per_rot;

r_theta(1,:) = theta;
r_phi(1,:) = phi;
r_x(1,:) =  x;
r_y(1,:) =  y;
r_z(1,:) =  z;

for idx = 2:frames_per_rot
    r_theta(idx,:) = r_theta(idx-1,:)+shift_per_frame;
    r_theta(idx,r_theta(idx,:)>=2*pi) = r_theta(idx,r_theta(idx,:)>=2*pi)-2*pi;
    r_phi(idx,:) = r_phi(idx-1,:);
    
    r_x(idx,:) = X_c + r.*cos(r_phi(idx,:)).*cos(r_theta(idx,:));
    r_y(idx,:) = Y_c + r.*cos(r_phi(idx,:)).*sin(r_theta(idx,:));
    r_z(idx,:) = Z_c + r.*sin(r_phi(idx,:));    
end

%% color ramp: number of frames for ramping stimuli in/out ant the beginning/end of each trial

n_frames_ramp = 40;

%% size of dots

dot_w = 0.12; Results.dot_w = dot_w; % dot-size in degree visual angle
dot_type = 1; %0: no anti-aliasing; 1: favors performance; 2: high quality anti-aliasing

% radius of fixation point (deg)
fix_r       = 0.12; Results.fix_r = fix_r; % size of fixation dot in degree visual angle

%% size of display

size_factor = 300; Results.Stimulus.size_factor = size_factor; % determines size of sphere in pixels

%% stereo

stereo_d = floor(frames_per_rot/100); Results.stereo_d = stereo_d; % shift of disambiguated dots vs. each other. 

%% Monitor

mon_width   = 39;  Results.Monitor.mon_width = mon_width; % horizontal dimension of viewable screen (cm)

v_dist      = 158;  Results.Monitor.v_dist = v_dist; % viewing distance (cm)

waitframes = 1;

%% delete unnecessary variables

clear cutout r_phi r_theta
