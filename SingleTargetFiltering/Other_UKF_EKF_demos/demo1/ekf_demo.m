%%-------------------------------------------------------------------
% EKF DEMO
%
% Simple nonlinear state estimation problem to demonstrate the
% Extented Kalman Filter
%
% AUTHOR : 
% DATE   : 
%%-------------------------------------------------------------------
clear all;
clc;
echo off;

%%-------------------------------------------------------------------
  fprintf('EKF  Demo : State estimation experiment\n');
  fprintf('---------------------------------------\n\n');
  fprintf('EKF of vehichle moving along a circular arc with a fixed radius. The speed\n');
  fprintf('of the vehicle along the arc is perturbed by white Gaussian noise\n');
  fprintf('source. The only observation of the vehicle position is the bearing\n');
  fprintf('angle.\n\n');
%%-------------------------------------------------------------------

number_of_runs = input('Number of runs : ');         %

mean_RMSE_ekf = zeros(1,number_of_runs);             %

for j=1:number_of_runs,                              %

rand('state',sum(100*clock));                        % Shuffle the pack!
randn('state',sum(100*clock));                       % Shuffle the pack!

radius = 10;                                         %
speed  = 3;                                          %
dt     = 1;                                          %

u = [radius; dt];                                    % control inputs

N = 30;                                              % number of time steps

x0  = [speed; radius; 0];                            % initial state 

P0  = 1*[1 0 0;                                      %
         0 1 0;                                      %
         0 0 1];                                     % initial state covariance 

L = size(x0,1);                                      % state dimension

Q    = (0.01*speed)^2;                               % process noise variance  
%R    = 0.1;                                          % measurement noise variance
R    = 0.01;                                          % measurement noise variance

xh  = zeros(L,N+1);                                  % state estimate buffer
P   = zeros(L,L,1,N+1);                              % state covariance buffer

xh(:,1)    = x0;                                     % initialize buffers
Px(:,:,1)  = P0;

xh_ekf = xh;                                         % create EKF buffers from template
P_ekf  = Px;                                         %

alpha = 1;                                           %
beta  = 2;                                           %
kappa = 0;                                           % 3 - state dimension

%%-------------------------------------------------------------------
%% GENERATE DATASET 
%%-------------------------------------------------------------------
fprintf('\nGenerating data...\n');                   %

x      = zeros(L,N+1);                               %
y      = zeros(1,N+1);                               %
v      = sqrt(Q)*randn(1,N+1);                       % process noise 
n      = sqrt(R)*randn(1,N+1);                       % measurement noise
                                 
x(:,1) = x0;                                         % initial state condition

y(:,1) = feval('hfun1',x(:,1),u,n(:,1),1);           % initial onbservation of state
tr(:,1) = feval('hfun1',x(:,1),u,0,1);               % initial onbservation of state


for k=2:(N+1),                                       %
  x(:,k) = feval('ffun1',x(:,k-1),u,v(:,k),k);       %
  y(:,k) = feval('hfun1',x(:,k),u,n(:,k),k);         %
  tr(:,k) = feval('hfun1',x(:,k),u,0,1);             % initial onbservation of state
end                                                  % 


%%-------------------------------------------------------------------
%% ESTIMATE STATE USING UKF 
%%-------------------------------------------------------------------
fprintf('\nEstimating trajectory...\n');             %
for k=2:(N+1),                                       % Generate EKF estimate
  xPred_ekf = feval('ffun1',xh_ekf(:,k-1),u,0,k);    % Ax - predicted next step
  %xPred_ekf = feval('ffun1',xh_ekf(:,k-1),u,v(:,k),k); % Ax - predicted next step
  Jx = jacobian_ffun1(xh_ekf(:,k-1),u);              % A - Jacobian - transition matrix
  PPred_ekf = diag([Q 0 0]) + Jx*P_ekf(:,:,k-1)*Jx'; % APA+Q - propagated covariance
  %yPred = feval('hfun1',xPred_ekf,u,v(:,k),k);      % Hx - measured output
  yPred = feval('hfun1',xPred_ekf,u,0,k);            % Hx - measured output
  Jy = jacobian_hfun1(xPred_ekf,u);                  % H - Jacobian - measurement matrix
  S  = R + Jy*PPred_ekf*Jy';                         % (HPH'+R)
  Si = inv(S);                                       % inv(HPH'+R)
  K = PPred_ekf*Jy'*Si;                              % K - calculated Kalman Gain
  xh_ekf(:,k)  = xPred_ekf + K*(y(:,k)-yPred);       % 
  P_ekf(:,:,k) = PPred_ekf - K*Jy*PPred_ekf;         %  
end                                                  % 

%%-------------------------------------------------------------------
%% CALCULATE ERRORS 
%%-------------------------------------------------------------------
error_ekf = (x(:,2:end)-xh_ekf(:,2:end)).^2;         %
RMSE_ekf  = (sum(error_ekf).^0.5);                   %
mean_RMSE_ekf(j) = mean(RMSE_ekf);                   %
fprintf('\n\nEKF estimate normalized RMSE : %2.4f\n',mean_RMSE_ekf(j)/radius); %

%%-------------------------------------------------------------------
%% DISPLAY RESULTS 
%%-------------------------------------------------------------------

of = 1;                                              %
figure(1); clf;                                      %
axis([-10,10,-10,10]);                               %

p1 = plot(x(2,:),x(3,:),'bo','linewidth',0.5); hold on %

p2 = plot((radius+of)*cos(y),(radius+of)*sin(y),'k+');  %
p3 = plot((radius+of)*cos(y),(radius+of)*sin(y),'k');   %

p4 = plot(xh_ekf(2,:),xh_ekf(3,:),'r^');             %
p5 = plot(xh_ekf(2,:),xh_ekf(3,:),'r');              %

p6 = plot(0,0,'r+');                                 %

legend([p1 p2 p3 p4 p5],'true','y','y','ekf','ekf'); %
title('Circular motion with WGN perturbed speed','fontsize',14); %
%axis(2*[-radius radius -radius radius]);             %

figure(2);                                           %
p1=plot(RMSE_ekf,'r'); hold on;                      %
legend([p1],'EKF');                                  %
title('RMS Tracking Error','fontsize',14);           %
xlabel('k');                                         %
ylabel('RMSE');                                      %

drawnow                                              %

end                                                  %

