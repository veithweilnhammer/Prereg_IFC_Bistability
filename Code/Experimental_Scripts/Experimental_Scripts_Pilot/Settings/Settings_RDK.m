%% Experimental Duration
fixation_interval = 5;
fraction_of_disambiguated_dots = [1 0.75 0.5 0.25 0.125 0.0625 0.03125];
trial_n = length(fraction_of_disambiguated_dots) * 2;

disambiguation(2:2:trial_n) = fraction_of_disambiguated_dots(randperm(length(fraction_of_disambiguated_dots)));

rot_per_trial = 5;
frames_per_rot = 1000; % number of frames per rotation

estimated_ifi = 0.0118; % 0.0167

estimated_experimental_duration = (estimated_ifi*frames_per_rot*rot_per_trial*trial_n + fixation_interval*(trial_n + 1))/60;

%% initialize sphere
max_number_of_dots = 5000;
step_size_dots = 2*pi/max_number_of_dots;
r = 1;
X_c = 0; Y_c = 0; Z_c = 0;

theta = [0:step_size_dots:2*pi-step_size_dots];
phi = pi.*rand(1, length(theta)) - 0.5*pi;

x = X_c + r.*cos(phi).*cos(theta);
y = Y_c + r.*cos(phi).*sin(theta);
z = Z_c + r.*sin(phi);

%%  cut out slices
fraction_of_surface = 0.3;
cutout = find(abs(x) > fraction_of_surface & abs(y) > fraction_of_surface); 
z(cutout) = [];
y(cutout) = [];
x(cutout) = [];
phi(cutout) = [];
theta(cutout) = [];

%% overlapping configurations
n_overlaps = 8;
lag_between_overlaps = frames_per_rot / n_overlaps;
overlaps = [1 : lag_between_overlaps : frames_per_rot*rot_per_trial];

%% rotate sphere
shift_per_frame = 2*pi/frames_per_rot;

r_theta(1,:) = theta;
r_phi(1,:) = phi;
r_x(1,:) =  x;
r_y(1,:) =  y;
r_z(1,:) =  z;

 min_color = 0.2;
 max_color = 0.8; 
 
  color_coding(1,:) = ((y+1)./2).*(max_color-min_color)+min_color;

for idx = 2:frames_per_rot
    r_theta(idx,:) = r_theta(idx-1,:)+shift_per_frame;
    r_theta(idx,r_theta(idx,:)>=2*pi) = r_theta(idx,r_theta(idx,:)>=2*pi)-2*pi;
    r_phi(idx,:) = r_phi(idx-1,:);
    
    r_x(idx,:) = X_c + r.*cos(r_phi(idx,:)).*cos(r_theta(idx,:));
    r_y(idx,:) = Y_c + r.*cos(r_phi(idx,:)).*sin(r_theta(idx,:));
    r_z(idx,:) = Z_c + r.*sin(r_phi(idx,:));
    
    color_coding(idx,:) = ((r_y(idx,:)+1)./2).*(max_color-min_color)+min_color;  
end

%% color ramp
n_frames_ramp = 40;

%% size of dots
dot_w = 0.15;
dot_type = 1; %0: no anti-aliasing; 1: favors performance; 2: high quality anti-aliasing

% radius of fixation point (deg)
fix_r       = 0.10;

%% size of display
size_factor = 300;

%% stereo
stereo_d = floor(frames_per_rot/200);

%% Monitor
mon_width   = 35;   % horizontal dimension of viewable screen (cm)

screen_res_width=1024; screen_res_height=768;
screen_width=365; screen_height=275; obs_dist=595;

v_dist      = 59.5;   % viewing distance (cm) war 60cm

waitframes = 1;

