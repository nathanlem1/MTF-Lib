% UKF DEMONSTRATION
% -----------------
%
% This is a simple nonlinear state estimation problem to demonstrate the
% use of the Unscented Kalman Filter (UKF) vs. the Extented Kalman Filter
% (EKF) on a standard nonlinear state estimation problem.
%
% AUTHOR : Rudolph van der Merwe (rvdmerwe@ieee.org)
% DATE   : 24 May 2001
%

clear all;
clc;
echo off;
path('../ukf',path);

fprintf('UKF vs. EKF  Demo 2 : State estimation experiment\n');
fprintf('-------------------------------------------------\n\n');
fprintf('In this demo the EKF and UKF are used to estimate a simple\n');
fprintf('nonlinear scalar time series.\n\n');

number_of_runs = input('Number of runs : ');

mean_RMSE_ekf = zeros(1,number_of_runs);
mean_RMSE_ukf = zeros(1,number_of_runs);

for j=1:number_of_runs,

rand('state',sum(100*clock));   % Shuffle the pack!
randn('state',sum(100*clock));   % Shuffle the pack!

u = 0;      % control inputs

N = 100;                % number of time steps

x0  = 1;  % initial state 

P0  = 1;                  % initial state covariance 

L = size(x0,1);           % state dimension

Q    = 0.1;              % process noise variance  
R    = 0.01;              % measurement noise variance

xh  = zeros(L,N+1);    % state estimate buffer
P   = zeros(L,L,1,N+1);  % state covariance buffer

xh(:,1)    = x0;       % initialize buffers
Px(:,:,1)  = P0;

xh_ukf = xh;           % create UKF buffers from template    
P_ukf  = Px;
xh_ekf = xh;           % create EKF buffers from template
P_ekf  = Px;

alpha = 1;
beta  = 0;
kappa = 2;           % 3 - state dimension

%%-------------------------------------------------------------------
%%---------------------- GENERATE DATASET ---------------------------

fprintf('\nGenerating data...\n');

x      = zeros(L,N+1);
y      = zeros(1,N+1);
v      = sqrt(Q)*randn(1,N+1);    % process noise 
n      = sqrt(R)*randn(1,N+1);    % measurement noise

x(:,1) = x0;                    % initial state condition

y(:,1) = feval('hfun2',x(:,1),u,n(:,1),1); % initial onbservation of state

for k=2:(N+1),
  x(:,k) = feval('ffun2',x(:,k-1),u,v(:,k),k);
  y(:,k) = feval('hfun2',x(:,k),u,n(:,k),k);
end

%%-------------------------------------------------------------------
%%------------------- ESTIMATE STATE USING UKF ----------------------
fprintf('\nEstimating trajectory...\n');

for k=2:(N+1),

  % Generate EKF estimate
  xPred_ekf = feval('ffun2',xh_ekf(:,k-1),u,0,k);  % EKF predicted mean
  Jx = feval('jacobian_ffun2',xh_ekf(:,k-1),u);             % Jacobian for ffun1
  PPred_ekf = Q + Jx*P_ekf(:,:,k-1)*Jx';           % EKF predicted state covariance
  yPred = feval('hfun2',xPred_ekf,u,0,k);
  Jy = feval('jacobian_hfun2',xPred_ekf,u);                % Jacobian for hfun1
  S  = R + Jy*PPred_ekf*Jy';
  Si = inv(S);
  K = PPred_ekf*Jy'*Si;
  xh_ekf(:,k)  = xPred_ekf + K*(y(:,k)-yPred);
  P_ekf(:,:,k) = PPred_ekf - K*Jy*PPred_ekf;
 
  % Generate UKF estimate
  [xh_ukf(:,k),P_ukf(:,:,k)] = ukf(xh_ukf(:,k-1),P_ukf(:,:,k-1),u,Q,'ffun2',y(:,k),R,'hfun2',k,alpha,beta,kappa);
  
end

%%-------------------------------------------------------------------
%%------------------- CALCULATE ERRORS ------- ----------------------

error_ekf = (x(:,2:end)-xh_ekf(:,2:end)).^2;
RMSE_ekf  = error_ekf.^0.5;
mean_RMSE_ekf(j) = mean(RMSE_ekf);

error_ukf = (x(:,2:end)-xh_ukf(:,2:end)).^2;
RMSE_ukf  = error_ukf.^0.5;
mean_RMSE_ukf(j) = mean(RMSE_ukf);

fprintf('\n\nEKF estimate normalized RMSE : %2.4f\n',mean_RMSE_ekf(j)/var(y));
fprintf('UKF estimate normalized RMSE : %2.4f\n\n',mean_RMSE_ukf(j)/var(y));

%%-------------------------------------------------------------------
%%------------------- DISPLAY RESULTS -------------------------------

figure(1); clf;
p1 = plot(x(2:end),'b-o','linewidth',1.5); hold on
p2 = plot(y(2:end),'k+');
p3 = plot(xh_ekf(2:end),'r');
p4 = plot(xh_ukf(2:end),'g'); hold off
legend([p1 p2 p3 p4],'true state','noisy obs','EKF estimate','UKF estimate');
xlabel('k');
title('Nonlinear time series estimation');
drawnow

end

fprintf('\n\n');
fprintf('---------------------------------------------------------\n');
fprintf('Mean & Variance of normalized RMSE over %d runs\n\n',number_of_runs);
fprintf('EKF : %2.4f (%2.4f)\n',mean(mean_RMSE_ekf/var(y)),var(mean_RMSE_ekf/var(y)));
fprintf('UKF : %2.4f (%2.4f)\n\n',mean(mean_RMSE_ukf/var(y)),var(mean_RMSE_ukf/var(y)));
fprintf('---------------------------------------------------------\n');


