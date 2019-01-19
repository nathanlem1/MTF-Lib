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

fprintf('UKF vs. EKF  Demo 1 : State estimation experiment\n');
fprintf('-------------------------------------------------\n\n');
fprintf('In this demo the EKF and UKF are used to track a vehicle which\n');
fprintf('is moving along a circular arc with a fixed radius. The speed\n');
fprintf('of the vehicle along the arc is perturbed by white Gaussian noise\n');
fprintf('source. The only observation of the vehicle position is the bearing\n');
fprintf('angle.\n\n');

number_of_runs = input('Number of runs : ');

mean_RMSE_ekf = zeros(1,number_of_runs);
mean_RMSE_ukf = zeros(1,number_of_runs);

for j=1:number_of_runs,

rand('state',sum(100*clock));   % Shuffle the pack!
randn('state',sum(100*clock));   % Shuffle the pack!

radius = 10;
speed  = 1;
dt     = 1;

u = [radius; dt];      % control inputs

N = 100;                % number of time steps

x0  = [speed; radius; 0];  % initial state 

P0  = 1*[1 0 0;
         0 1 0;
         0 0 1];                  % initial state covariance 

L = size(x0,1);                % state dimension

Q    = (.1*speed)^2;            % process noise variance  
R    = 0.1;                    % measurement noise variance

xh  = zeros(L,N+1);    % state estimate buffer
P   = zeros(L,L,1,N+1);  % state covariance buffer

xh(:,1)    = x0;       % initialize buffers
Px(:,:,1)  = P0;

xh_ukf = xh;           % create UKF buffers from template    
P_ukf  = Px;
xh_ekf = xh;           % create EKF buffers from template
P_ekf  = Px;

alpha = 1;
beta  = 2;
kappa = 0;           % 3 - state dimension

%%-------------------------------------------------------------------
%%---------------------- GENERATE DATASET ---------------------------

fprintf('\nGenerating data...\n');

x      = zeros(L,N+1);
y      = zeros(1,N+1);
v      = sqrt(Q)*randn(1,N+1);    % process noise 
n      = sqrt(R)*randn(1,N+1);    % measurement noise

x(:,1) = x0;                    % initial state condition

y(:,1) = feval('hfun1',x(:,1),u,n(:,1),1); % initial onbservation of state

for k=2:(N+1),
  x(:,k) = feval('ffun1',x(:,k-1),u,v(:,k),k);
  y(:,k) = feval('hfun1',x(:,k),u,n(:,k),k);
end


%%-------------------------------------------------------------------
%%------------------- ESTIMATE STATE USING UKF ----------------------
fprintf('\nEstimating trajectory...\n');

for k=2:(N+1),

  % Generate EKF estimate
  xPred_ekf = feval('ffun1',xh_ekf(:,k-1),u,0,k);  % EKF predicted mean
  Jx = jacobian_ffun1(xh_ekf(:,k-1),u);             % Jacobian for ffun1
  PPred_ekf = diag([Q 0 0]) + Jx*P_ekf(:,:,k-1)*Jx';           % EKF predicted state covariance
  yPred = feval('hfun1',xPred_ekf,u,0,k);
  Jy = jacobian_hfun1(xPred_ekf,u);                % Jacobian for hfun1
  S  = R + Jy*PPred_ekf*Jy';
  Si = inv(S);
  K = PPred_ekf*Jy'*Si;
  xh_ekf(:,k)  = xPred_ekf + K*(y(:,k)-yPred);
  P_ekf(:,:,k) = PPred_ekf - K*Jy*PPred_ekf;
 
  % Generate UKF estimate
  [xh_ukf(:,k),P_ukf(:,:,k)] = ukf(xh_ukf(:,k-1),P_ukf(:,:,k-1),u,Q,'ffun1',y(:,k),R,'hfun1',k,alpha,beta,kappa);
  
end

%%-------------------------------------------------------------------
%%------------------- CALCULATE ERRORS ------- ----------------------

error_ekf = (x(:,2:end)-xh_ekf(:,2:end)).^2;
RMSE_ekf  = (sum(error_ekf).^0.5);
mean_RMSE_ekf(j) = mean(RMSE_ekf);

error_ukf = (x(:,2:end)-xh_ukf(:,2:end)).^2;
RMSE_ukf  = (sum(error_ukf).^0.5);
mean_RMSE_ukf(j) = mean(RMSE_ukf);

fprintf('\n\nEKF estimate normalized RMSE : %2.4f\n',mean_RMSE_ekf(j)/radius);
fprintf('UKF estimate normalized RMSE : %2.4f\n\n',mean_RMSE_ukf(j)/radius);

%%-------------------------------------------------------------------
%%------------------- DISPLAY RESULTS -------------------------------

figure(1); clf;
subplot(311);
p1 = plot(x(2,:),x(3,:),'bo','linewidth',1.5); hold on
p2 = plot(radius*cos(y),radius*sin(y),'k+');
p3 = plot(xh_ekf(2,:),xh_ekf(3,:),'r^');
p4 = plot(xh_ukf(2,:),xh_ukf(3,:),'gd'); hold off
legend([p1 p3 p4],'true state','EKF estimate','UKF estimate');
title('Circular motion with WGN perturbed speed','fontsize',14);
axis(2*[-radius radius -radius radius]);
subplot(312);
p1=plot(RMSE_ekf,'r'); hold on;
p2=plot(RMSE_ukf,'g'); hold off
legend([p1 p2],'EKF','UKF');
title('RMS Tracking Error','fontsize',14);
xlabel('k');
ylabel('RMSE');
subplot(313);
plot(x(1,:),'b','linewidth',2);
axis tight
xlabel('k');
title('Speed');
ylabel('arc length per second');
drawnow

end

fprintf('\n\n');
fprintf('---------------------------------------------------------\n');
fprintf('Mean & Variance of normalized RMSE over %d runs\n\n',number_of_runs);
fprintf('EKF : %2.4f (%2.4f)\n',mean(mean_RMSE_ekf/radius),var(mean_RMSE_ekf/radius));
fprintf('UKF : %2.4f (%2.4f)\n\n',mean(mean_RMSE_ukf/radius),var(mean_RMSE_ukf/radius));
fprintf('---------------------------------------------------------\n');