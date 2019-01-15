% This script runs KF example
%
% Version 1.0, January, 2019
%
% This function was written by Nathanael L. Baisa
% email: nathanaellmss@gmail.com
%
%
% State vector is position x,y and velcity vx, vy i.e. [x, vx, y,vy] and 
% observation vector is position x, y i.e. [x, y]

clc;
clear;
close all;

% dynamical model parameters (CV model)
model.T= 1;                                           %sampling period
model.A0= [ 1 model.T; 0 1 ];                                            
model.F= [ model.A0 zeros(2,2); zeros(2,2) model.A0 ]; %transition matrix  
model.B0= [ (model.T^2)/2; model.T ];
model.B= [ model.B0 zeros(2,1); zeros(2,1) model.B0 ]; %process noise standard deviation
model.sigma_v = 2.5;
model.Q= (model.sigma_v)^2* model.B*model.B';   %process noise covariance

% observation model parameters (noisy x/y only)
model.H= [ 1 0 0 0 ; 0 0 1 0 ];    %observation matrix
model.D= diag([ 10; 10 ]);         %observation noise standard deviation
model.R= model.D*model.D';         %observation noise covariance

% Birth parameters 
m_birth = [ 10; 0; -10; 0 ];   %mean of Gaussian birth term 
B_birth = diag([ 10; 10; 10; 10 ]); %std of Gaussian birth term
P_birth = B_birth*B_birth';  %cov of Gaussian birth term

xstart  = m_birth; %[ 0; 0; 0; -10 ]; 
targetstate = xstart;

N_duration = 100; %length of data/number of scans
truth_X= cell(N_duration,1);  %ground truth for state of target
meas_Z= cell(N_duration,1);   %generated measurement of a target

% Initialization
m_update = m_birth;
P_update = P_birth;

% Storing estimated states
estimated_X = cell(N_duration,1);

figure(1), hold on;

for k = 1:N_duration
    %targetstate = model.F * targetstate + sqrt(model.Q)*randn(4,1); % generate ground truth for states of a target
    W = model.sigma_v*model.B*randn(size(model.B,2),size(targetstate,2)); % equivalent to 'sqrt(model.Q)*randn(4,1)'
    targetstate = model.F * targetstate + W;  % generate ground truth for states of a target
    truth_X{k}= targetstate;
    
    %meas_Z{k}= model.H*targetstate + sqrt(model.R)*randn(size(model.R,2),1); % generate measurement of a target
    V = model.D*randn(size(model.D,2),size(targetstate,2)); % equivalent to 'sqrt(model.R)*randn(size(model.R,2),1)'
    meas_Z{k}= model.H*targetstate + V; % generate measurement of a target
    
    [m_update, P_update] = KalmanFilter(meas_Z{k}, model, m_update, P_update);
    
    estimated_X{k} = m_update;
    
    %Plot the ground truth
    h1 = plot(truth_X{k}(1), truth_X{k}(3), '.r', 'LineWidth', 2);
    %Plot the measurement
    h2 = plot(meas_Z{k}(1), meas_Z{k}(2), '+b', 'LineWidth', 2);
    %Plot the estimated state
    h3 = plot(estimated_X{k}(1), estimated_X{k}(3), '.g', 'LineWidth', 2);
    xlabel('X coordinate');
    ylabel('Y coordinate');
    legend([h1, h2, h3], 'ground truth', 'measurement', 'estimated state')
    axis square;  
end













