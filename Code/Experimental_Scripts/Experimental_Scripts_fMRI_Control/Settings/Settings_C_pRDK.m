
%% Experimental Duration

fixation_interval = 10; Results.fixation_interval = fixation_interval;

rot_per_trial = 10; Results.rot_per_trial = rot_per_trial; % number of rotations of the RDK per block.

desired_time_per_rotation = 12; % desired time per cycle.

estimated_ifi = 1/60; % 0.0167 % expected inter-frame interval.

frames_per_cycle = round(desired_time_per_rotation/estimated_ifi); Results.Stimulus.frames_per_cycle = frames_per_cycle; % number of frames per rotation

estimated_experimental_duration = (estimated_ifi*frames_per_cycle*rot_per_trial*trial_n + fixation_interval*(trial_n + 1))/60;

%% initialize sphere

max_number_of_dots = 5000 * 0.3 * 0.5; % number of dots to plot (equivalent to number of undirectionally moving dots from main experiment
Results.Stimulus.max_numer_of_dots = max_number_of_dots;
step_size_dots = 2*pi/max_number_of_dots;
r = 1; % radius of circle
X_c = 0; Y_c = 0; Z_c = 0; % center of circle

% initial positioning of random dots
theta = [0:step_size_dots:2*pi-step_size_dots];
phi = pi.*rand(1, length(theta)) - 0.5*pi;

x = X_c + r.*cos(phi).*cos(theta);
y = Y_c + r.*cos(phi).*sin(theta);
z = Z_c + r.*sin(phi);

%% overlapping configurations

% discretize perceptual timecourse in analogy to the main experiment
n_overlaps = 8; Results.n_overlaps = n_overlaps;
lag_between_overlaps = frames_per_cycle / n_overlaps;
overlaps = [1 : lag_between_overlaps : frames_per_cycle*rot_per_trial]; Results.overlaps = overlaps;

minimal_response_time = 0.2; Results.minimal_response_time = minimal_response_time;
minimal_response_frames = round(minimal_response_time/estimated_ifi);

%% move dots
shift_per_frame = 2*pi/frames_per_cycle; % translational dot-speed in analogy to the main expermient

%% color ramp: ramping in/out at beginning/end of experients
n_frames_ramp = 40;

%% size of dots
dot_w = 0.12; Results.dot_w = dot_w; % dot-size in degree visual agnle
dot_type = 1; %0: no anti-aliasing; 1: favors performance; 2: high quality anti-aliasing

% radius of fixation point (deg)
fix_r       = 0.12; Results.fix_r = fix_r;

%% size of display

size_factor = 300; Results.Stimulus.size_factor = size_factor; % determines size of circle

%% Monitor

mon_width   = 36.5;  Results.Monitor.mon_width = mon_width; % horizontal dimension of viewable screen (cm)

v_dist      = 60;  Results.Monitor.v_dist = v_dist; % viewing distance (cm)

waitframes = 1;

%% delete unnecessary variables

clear cutout r_phi r_theta
