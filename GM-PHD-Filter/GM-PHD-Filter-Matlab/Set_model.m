% This function sets important model parameters
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: 
% Modified: March, 2014
% -------------------------------------------------------------------------

function model = Set_model()

% Dynamic model
I2 = eye(2); % 2x2 identify matrix, used to construct matrices
Z2 = zeros(2); % 2x2 zero matrix, used to construct matrices
dt = 1; % One-second sampling period

model.sigma_v = 1.5; %5 % Standard deviation of process noise
model.F = [ [I2, dt*I2]; [Z2 I2] ]; % State transition matrix (motion model)
model.Q = model.sigma_v^2 * [ [1/4*dt^4*I2, 1/2*dt^3*I2]; [1/2*dt^3* I2, dt^2*I2] ]; % Process noise covariance

% Observation model
model.sigma_r = 2.0; %2; 10 % Standard deviation of measurement noise

model.H = [I2, Z2]; % Observation matrix 
model.R = model.sigma_r^2 * I2; % Sensor noise covariance.

% Initial covariance for all targets.
covariance_birth = diag([100, 100, 25, 25]');
w_birthsum = 0.0001;  % The total weight of birthed targets. It is chosen depending on handling false positives.

model.P = covariance_birth;
model.w_birthsum = w_birthsum;
%% Important parameters
prob_detection = 0.98; % Probability of target detection
prob_survival = 0.99; % Probability of target survival (prob_death = 1 - prob_survival)

model.p_S = prob_survival;
model.p_D = prob_detection;

%% Generate intensity of clutter RFS K_k at time k 
nClutter = poissrnd(10); % Number of clutter points, clutter is Poisson distributed  
xrange = [-1000 1000]; % X range of measurements
yrange = [-1000 1000]; % Y range of measurements
V = ( (xrange(2)-xrange(1)) * (yrange(2)-yrange(1))); %4 * 10^6;  % Volume of surveillance region
lambda_c = nClutter/V; % Average clutter returns per unit volume (nClutter returns over the region/frame)
u_z = 1.0/( (xrange(2)-xrange(1)) * (yrange(2)-yrange(1)));
clutter_intensity = @(z_cartesian) lambda_c * V * u_z; % Generate clutter function.
model.clutterIntensity = clutter_intensity;
model.lambda_t = nClutter;
model.xrange = xrange;
model.yrange = yrange;

model.T = 10^-5; % Pruning weight threshold
model.U = 4; % Merge distance threshold
model.w_thresh = 0.5; % State extraction weight threshold

%% Calculate performance metric
% The order_p and cutoff_c control filter performance.
model.order_p = 1; % The order 
model.cutoff_c = 100; % Cutoff determines the maximum error for a point.




end

