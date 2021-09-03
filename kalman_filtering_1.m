%1D KALMAN FILTERING 
clc
clear all
close all
%% define variables (i.e. how long and often we will sample)
duration = 10; %how long the Bird flies
dt = .1; %Sampling interval

%% Define update equations (Coefficent matrices): A physics based model for where we expect the Bird to be
% [state transition (state + velocity)] + [input control (acceleration)]
A = [1 dt; 0 1] ; % state transition matrix: expected flight of the Bird (state prediction)
B = [dt^2/2; dt]; %input control matrix: expected effect of the input accceleration on the state.
C = [1 0]; % measurement matrix: the expected measurement given the predicted state (likelihood)
%since we are only measuring position we set the velocity variable to zero in C matrix.

%% define main variables
u = 1.5; % define acceleration magnitude
Q= [0; 0]; %initial state-has two components: [position; velocity] of the Bird
Q_estimate = Q; %x_estimate of initial location estimation of the Bird 
BirdAccel_noise_mag = 0.05; %process noise: the variation in the acceleration of the bird Bird (stdv of acceleration: meters/sec^2)
Sensor_noise_mag = 10; %measurement noise: How much error introduced by the sensor ex:RADAR (stdv of location, in meters)
Ez = Sensor_noise_mag^2;% Ez convert the measurement noise (stdv) into covariance matrix
Ex = BirdAccel_noise_mag^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]; % Ex convert the process noise (stdv) into covariance matrix
P = Ex; % estimate of initial Bird position variance (covariance matrix)

%% initialize result variables
% Initialize for speed
Q_loc = []; % ACTUAL Bird flight path
vel = []; % ACTUAL Bird velocity
Q_loc_meas = []; % Bird path that the is estimated 



%% simulate what the Sensor sees over time
figure(2);clf
figure(1);clf
for t = 0 : dt: duration

% Generate the Bird flight
BirdAccel_noise = BirdAccel_noise_mag * [(dt^2/2)*randn; dt*randn];
Q= A * Q+ B * u + BirdAccel_noise;
% Generate what the Sensor sees
SensorVision_noise = Sensor_noise_mag * randn;
y = C * Q+ SensorVision_noise;
Q_loc = [Q_loc; Q(1)];
Q_loc_meas = [Q_loc_meas; y];
vel = [vel; Q(2)];
%iteratively plot what the Sensor sees
plot(0:dt:t, Q_loc, '-r.')
plot(0:dt:t, Q_loc_meas, '-k.')
axis([0 10 -30 80])
hold on
% pause
end

%plot theoretical path of Sensor that doesn't use kalman using averaging between sampled points
plot(0:dt:t, smooth(Q_loc_meas), '-g.')

%plot velocity, just to show that it's constantly increasing, due to constant acceleration
%figure(2);
%plot(0:dt:t, vel, '-b.')



%% Do kalman filtering
%initize estimation variables
Q_loc_estimate = []; % Bird position estimate
vel_estimate = []; % Bird velocity estimate
Q= [0; 0]; % re-initized state
P_estimate = P;
P_mag_estimate = [];
predic_state = [];
predic_var = [];
for t = 1:length(Q_loc)
% Predict next state of the Bird with the last state and predicted motion.
Q_estimate = A * Q_estimate + B * u;
predic_state = [predic_state; Q_estimate(1)] ;
%predict next covariance
P = A * P * A' + Ex;
predic_var = [predic_var; P] ;
% predicted Sensor measurement covariance
% Kalman Gain
K = P*C'*inv(C*P*C'+Ez);
% Update the state estimate.
Q_estimate = Q_estimate + K * (Q_loc_meas(t) - C * Q_estimate);
% update covariance estimation.
P = (eye(2)-K*C)*P;
%Store for plotting
Q_loc_estimate = [Q_loc_estimate; Q_estimate(1)];
vel_estimate = [vel_estimate; Q_estimate(2)];
P_mag_estimate = [P_mag_estimate; P(1)];
end

% Plot the results
figure(2);
tt = 0 : dt : duration;
plot(tt,Q_loc,'-r.',tt,Q_loc_meas,'-k.', tt,Q_loc_estimate,'-g.');
axis([0 10 -30 80])

