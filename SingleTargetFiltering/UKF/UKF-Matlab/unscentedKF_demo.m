% This script runs unscentedKF example
%
% Version 1.0, January, 2019
%
% This function was written by Nathanael L. Baisa
% email: nathanaellmss@gmail.com
%
% State vector is position x,y and velcity vx, vy i.e. [x, vx, y, vy, w_turn] and 
% observation vector is position (x, y) converted into (theta, range) = (arctan(x/y); sqrt(x^2 + y^2)) to (x,y) form.
%

clc;
clear;
close all;

% basic parameters
model.x_dim= 5;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector
model.v_dim= 3;   %dimension of process noise
model.w_dim= 2;   %dimension of observation noise

% dynamical model parameters (Constant turn (CT) model)
% state transformation given by gen_newstate_fn, transition matrix is N/A in non-linear case
model.T= 1;                         %sampling period
model.sigma_vel= 5;
model.sigma_turn= (pi/180);   %std. of turn rate variation (rad/s)
model.bt= model.sigma_vel*[ (model.T^2)/2; model.T ];
model.B2= [ model.bt zeros(2,2); zeros(2,1) model.bt zeros(2,1); zeros(1,2) model.T*model.sigma_turn ];
model.B= eye(model.v_dim);
model.Q= model.B*model.B';  %process noise covariance

% observation model parameters (noisy r/theta only)
% measurement transformation given by gen_observation_fn, observation matrix is N/A in non-linear case
model.D= diag([ 2*(pi/180); 10 ]);      %std for angle and range noise
model.R= model.D*model.D';              %covariance for observation noise

% parameters for UKF
ukf_alpha= 1;                %scale parameter for UKF - choose alpha=1 ensuring lambda=beta and offset of first cov weight is beta for numerical stability
ukf_beta= 2;                 %scale parameter for UKF
ukf_kappa= 2;                %scale parameter for UKF (alpha=1 preferred for stability, giving lambda=1, offset of beta for first cov weight)

% Birth parameters 
m_birth = [ -100; 0; 250; 0; 0 ];   %mean of Gaussian birth term 
B_birth = diag([ 50; 50; 50; 50; 6*(pi/180) ]); %std of Gaussian birth term
P_birth = B_birth*B_birth';  %cov of Gaussian birth term

wturn = 2*pi/180;
xstart  = [ -150; 15; 250; -10; wturn/8 ];  
targetstate = xstart;

N_duration = 100; %length of data/number of scans
truth_X = [];  %ground truth for state of target
meas_Z = [];  %generated measurement of a target

% Initialization
m_update = m_birth;
P_update = P_birth;

% Storing estimated states
estimated_X = [];

figure(1), hold on;

for k = 1:N_duration
    
    % Generate ground truth for states of a target (track)
    targetstate = gen_newstate_fn(model,targetstate,'noise');
    truth_X = [truth_X targetstate];
    
    % Generate measurement of a target
    meas_z = gen_observation_fn(model,targetstate,'noise'); 
    meas_Z = [meas_Z meas_z];
    
    % Call uncented Kalman filter function
    [m_update, P_update] = unscentedKalmanFilter(meas_z, model, m_update, P_update, ukf_alpha, ukf_kappa, ukf_beta);
   
    estimated_X = [estimated_X m_update];
    
    % Plot the ground truth
    h1 = plot(truth_X(1,:), truth_X(3,:), '-r', 'LineWidth', 2);
    % Plot the measurement after converting from (theta, range) = (arctan(x/y) / sqrt(x^2 + y^2)) to (x,y)
    meas_Z_X = meas_Z(2,:).*sin(meas_Z(1,:));
    meas_Z_Y = meas_Z(2,:).*cos(meas_Z(1,:));
    h2 = plot(meas_Z_X, meas_Z_Y, '+b', 'LineWidth', 2);
    % Plot the estimated state
    h3 = plot(estimated_X(1,:), estimated_X(3,:), '.g', 'LineWidth', 2);
    xlabel('X coordinate');
    ylabel('Y coordinate');
    legend([h1, h2, h3], 'ground truth', 'measurement', 'estimated state')
    axis square; 
    drawnow; 
    
end






