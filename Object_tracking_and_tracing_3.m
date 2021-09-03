clear all;
close all;
clc;
set(0,'DefaultFigureWindowStyle','docked')
base_dir = '/MATLAB Drive/Test-ball-2';
cd(base_dir);

%% get listing of frames
f_list = dir('*png');

%% load tracking data
load('CM_idx_no.mat'); 

%% define main variables
dt = 0.1; %sampling rate (10fps)
S_frame = 1; %starting frame
u = 15; % define acceleration magnitude
Q= [CM_idx(S_frame,1); CM_idx(S_frame,2); 0; 0]; %initial state-it has four components: [positionX; positionY; velocityX; velocityY] of the ball
Q_estimate = Q; %estimate of initial location estimation of where the ball is 
HexAccel_noise_mag = 6; %process noise: the variability in how fast the ball is accelerating (stdv of acceleration: meters/sec^2)
tkn_x = 0.01; %measurement noise in the horizontal direction (x axis).
tkn_y = 0.01; %measurement noise in the horizontal direction (y axis).
Ez = [tkn_x 0; 0 tkn_y];
Ex = [dt^4/4 0 dt^3/2 0; ...
0 dt^4/4 0 dt^3/2; ...
dt^3/2 0 dt^2 0; ...
0 dt^3/2 0 dt^2].*HexAccel_noise_mag^2; % Ex convert the process noise (stdv) into covariance matrix
P = Ex; % estimate of initial ball position variance (covariance matrix)

%% Define update equations in 2-D (Coefficent matrices)
A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]; %state update matrix
B = [(dt^2/2); (dt^2/2); dt; dt];
C = [1 0 0 0; 0 1 0 0]; %this is measurement function C, Apply to the state estimate Q to get next/new measurement


%% initialize result variables
% Initialize for speed
Q_loc = []; % ACTUAL ball motion path
vel = []; % ACTUAL ball velocity
Q_loc_meas = []; % the ball path extracted by the camera


%% initize estimation variables
Q_loc_estimate = []; % position estimate
vel_estimate = []; % velocity estimate
P_estimate = P;
predic_state = [];
predic_var = [];
r = 5; % r is the radius of the plotting circle
j=0:.01:2*pi; %to make the plotting circle

for t = S_frame:length(f_list)
% load the image
img_tmp = double(imread(f_list(t).name));
img = img_tmp(:,:,1);
% load the given tracking
Q_loc_meas(:,t) = [ CM_idx(t,1); CM_idx(t,2)];
%% kalman filtering 
% Predict next state of the ball with the last state and predicted motion.
Q_estimate = A * Q_estimate + B * u;
predic_state = [predic_state; Q_estimate(1)] ;
%predict next covariance
P = A * P * A' + Ex;
predic_var = [predic_var; P] ;
% predicted measurement covariance
% Kalman Gain
K = P*C'*inv(C*P*C'+Ez);
% Update the state estimate.
if ~isnan(Q_loc_meas(:,t))
Q_estimate = Q_estimate + K * (Q_loc_meas(:,t) - C * Q_estimate);
end
% update covariance estimation.
P = (eye(4)-K*C)*P;
%% Store data
Q_loc_estimate = [Q_loc_estimate; Q_estimate(1:2)];
vel_estimate = [vel_estimate; Q_estimate(3:4)];
%% plot the images with the tracking
imagesc(img);
axis off
colormap(gray);
hold on;
plot(r*sin(j)+Q_loc_meas(2,t),r*cos(j)+Q_loc_meas(1,t),'.g'); % actual tracking
plot(r*sin(j)+Q_estimate(2),r*cos(j)+Q_estimate(1),'.r'); % kalman filtered tracking
hold off
pause(0.1); 
end
