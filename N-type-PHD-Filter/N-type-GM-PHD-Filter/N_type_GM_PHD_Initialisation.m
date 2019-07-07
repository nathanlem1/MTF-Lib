%% ------------------------------------------------------------------------
% N_type_GM_PHD_Initialisation.m
% -------------------------------------------------------------------------
%
% This file initialises most of the variables that we use for the N-type GM-PHD filter.
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

%% Control Variables
N_type_GM_PHD = 1; % Set to 1 if you want to use N-type GM-PHD filter (where N = 4) or set to 0 if you want to use N independent GM-PHD filters. This is used in update step
PLOT_ALL_MEASUREMENTS = 0; % Set to 1 to maintain the plot of the full measurement history, including clutter and error ellipses. Set to 0 to just plot the most recent measurements.
CALCULATE_OSPA_METRIC = 1; % Set to 1 to calculate and plot the OSPA performance metric for each time step. 
DRAW_VELOCITY = 0; % Set to 1 if you want to draw velocity arrows of targets.
simulateTargetStateAndMeasurement = 0; % set it to 1 if you simulate measurements or 0 if you use the saved measurement

%% Dynamic and observation models, initializing using structure
% We use a constant velocity (CV) motion model
% Play with the covariances R (for the measurement) and Q (for prediction). 

% Dynamic model:
I2 = eye(2); % 2x2 identify matrix, used to construct matrices
Z2 = zeros(2); % 2x2 zero matrix, used to construct matrices
dt = 1.0; % One-second sampling period

model.sigma_v1 = 5; % Standard deviation of process noise is 5 m/(s^2)
model.sigma_v2 = 5; 
model.sigma_v3 = 5; 
model.sigma_v4 = 5; 

% Target type 1
model.F1 = [ [I2, dt*I2]; [Z2 I2] ]; % State transition matrix (motion model)
model.Q1 = model.sigma_v1^2 * [ [1/4*dt^4*I2, 1/2*dt^3*I2]; [1/2*dt^3* I2, dt^2*I2] ]; % Process noise covariance

% Target type 2
model.F2 = [ [I2, dt*I2]; [Z2 I2] ]; 
model.Q2 = model.sigma_v2^2 * [ [1/4*dt^4*I2, 1/2*dt^3*I2]; [1/2*dt^3* I2, dt^2*I2] ]; 

% Target type 3
model.F3 = [ [I2, dt*I2]; [Z2 I2] ]; 
model.Q3 = model.sigma_v3^2 * [ [1/4*dt^4*I2, 1/2*dt^3*I2]; [1/2*dt^3* I2, dt^2*I2] ]; 

% Target type 4
model.F4 = [ [I2, dt*I2]; [Z2 I2] ];
model.Q4 = model.sigma_v4^2 * [ [1/4*dt^4*I2, 1/2*dt^3*I2]; [1/2*dt^3* I2, dt^2*I2] ]; 

% Observation model:
model.sigma_r11 = 6; % Standard deviation of measurement noise is 6m. Used in creating R matrix (below), for target type 1 by detector 1.
model.sigma_r12 = 6; % for target type 2 by detector 1.
model.sigma_r13 = 6; % for target type 3 by detector 1.
model.sigma_r14 = 6; % for target type 4 by detector 1.
model.sigma_r21 = 6; % for target type 1 by detector 2.
model.sigma_r22 = 6; % for target type 2 by detector 2.
model.sigma_r23 = 6; % for target type 3 by detector 2.
model.sigma_r24 = 6; % for target type 4 by detector 2.
model.sigma_r31 = 6; % for target type 1 by detector 3.
model.sigma_r32 = 6; % for target type 2 by detector 3.
model.sigma_r33 = 6; % for target type 3 by detector 3.
model.sigma_r34 = 6; % for target type 4 by detector 3.
model.sigma_r41 = 6; % for target type 1 by detector 4.
model.sigma_r42 = 6; % for target type 2 by detector 4.
model.sigma_r43 = 6; % for target type 3 by detector 4.
model.sigma_r44 = 6; % for target type 4 by detector 4.

% Target type 1 by detector 1
model.H11 = [I2, Z2]; % Observation matrix 
model.R11 = model.sigma_r11^2 * I2; % observation noise covariance.
% Target type 2 by detector 2
model.H22 = [I2, Z2]; 
model.R22 = model.sigma_r22^2 * I2;
% Target type 3 by detector 3
model.H33 = [I2, Z2]; 
model.R33 = model.sigma_r33^2 * I2;
% Target type 4 by detector 4
model.H44 = [I2, Z2]; 
model.R44 = model.sigma_r44^2 * I2;
% Target type 1 by detector 2, 3 & 4
model.H21 = [I2, Z2];
model.R21 = model.sigma_r21^2 * I2;
model.H31 = [I2, Z2]; 
model.R31 = model.sigma_r31^2 * I2;
model.H41 = [I2, Z2];
model.R41 = model.sigma_r41^2 * I2;
% Target type 2 by detector 1, 3 & 4
model.H12 = [I2, Z2];  
model.R12 = model.sigma_r12^2 * I2;
model.H32 = [I2, Z2]; 
model.R32 = model.sigma_r32^2 * I2;
model.H42 = [I2, Z2]; 
model.R42 = model.sigma_r42^2 * I2;
% Target type 3 by detector 1, 2 & 4
model.H13 = [I2, Z2];
model.R13 = model.sigma_r13^2 * I2;
model.H23 = [I2, Z2];
model.R23 = model.sigma_r23^2 * I2;
model.H43 = [I2, Z2];
model.R43 = model.sigma_r43^2 * I2;
% Target type 4 by detector 1, 2 & 3
model.H14 = [I2, Z2];
model.R14 = model.sigma_r14^2 * I2;
model.H24 = [I2, Z2];
model.R24 = model.sigma_r24^2 * I2;
model.H34 = [I2, Z2];
model.R34 = model.sigma_r34^2 * I2;

% Only for measurement simulation, especially for getting the same
% measurement values for the correlated targets.
model.F = [ [I2, dt*I2]; [Z2 I2] ];% State transition matrix (motion model)
model.I = eye(4);
model.sigma_v = 2;%10;
model.Q = model.sigma_v^2 * [ [1/4*dt^4*I2, 1/2*dt^3*I2]; [1/2*dt^3* I2, dt^2*I2] ]; % Process noise covariance, given in Vo&Ma.

%% These mean values are used for simulating measurements, not for birthing of targets which is done from current measurements adaptively.
birth_mean1 = [-100, 700, 0, 0]';
birth_mean2 = [-750, -100, 0, 0]';
birth_mean3 = [-200, 400, 0, 0]';
birth_mean4 = [-700, -400, 0, 0]';
birth_mean5 = [-400, 600, 0, 0]';
birth_mean6 = [-800, -600, 0, 0]';
birth_mean7 = [-500, -200, 0, 0]';
birth_mean8 = [700, 600, 0, 0]';
birth_mean9 = [-900, 100, 0, 0]';

birth_mean10 = [-800, 500, 0, 0]';
birth_mean11 = [-900, -200, 0, 0]';
birth_mean12 = [400, -600, 0, 0]';
birth_mean13 = [800, -600, 0, 0]';
birth_mean14 = [500, -700, 0, 0]';
birth_mean15 = [-700, -600, 0, 0]';
birth_mean16 = [900, -100, 0, 0]';
% Initial covariance for all targets.
covariance_birth = diag([200, 200, 100, 100]'); 
% covariance_birth = diag([100, 100, 50, 50]'); % seems better but all are almost the same!
% covariance_birth = diag([100, 100, 25, 25]'); % Original
% covariance_birth = diag([200, 200, 75, 75]');
w_birthsum = 0.00001;  % The total weight of birthed targets. It is chosen depending on handling false positives.

%% Important parameters
% For target type 1
prob_detection11 = 0.90; % Probability of detection of target 1 by detector 1
                        % Used in recalculating weights in GM_PHD_Update
prob_survival1 = 0.99; % Probability of target survival (prob_death1 = 1 - prob_survival1). Used in GM_PHD_Predict_Existing for weight calculation
% For target type 2
prob_detection22 = 0.92; % Probability of detection of target 2 by detector 2
prob_survival2 = 0.99; % Probability of target survival (prob_death2 = 1 - prob_survival2)

% For target type 3
prob_detection33 = 0.92; % Probability of detection of target 3 by detector 3
prob_survival3 = 0.99; % Probability of target survival (prob_death3 = 1 - prob_survival3)

% For target type 4
prob_detection44 = 0.91; % Probability of detection of target 4 by detector 4
prob_survival4 = 0.99; % Probability of target survival (prob_death4 = 1 - prob_survival4)

prob_detection21 = 0.6; % Probability of detection of target 1 by detector 2.
prob_detection31 = 0.6; % Probability of detection of target 1 by detector 3.
prob_detection41 = 0.6; % Probability of detection of target 1 by detector 4.
prob_detection12 = 0.6; % Probability of detection of target 2 by detector 1.
prob_detection32 = 0.6; % Probability of detection of target 2 by detector 3.
prob_detection42 = 0.6; % Probability of detection of target 2 by detector 4.
prob_detection13 = 0.6; % Probability of detection of target 3 by detector 1.
prob_detection23 = 0.6; % Probability of detection of target 3 by detector 2.
prob_detection43 = 0.6; % Probability of detection of target 3 by detector 4.
prob_detection14 = 0.6; % Probability of detection of target 4 by detector 1.
prob_detection24 = 0.6; % Probability of detection of target 4 by detector 2.
prob_detection34 = 0.6; % Probability of detection of target 4 by detector 3.

mk_k_minus_1_before_prediction = {};% Used in augmenting measurement vector to calculate velocity, for update.

%% Generate intensity of clutter RFS c_k at time k -- used for the update step
% Detected clutter is a Poisson RFS. This clutter intensity is used to model
% detections that are not associated with objects of interest (targets)i.e.
% false detections.

% Clutter intensity to background for target type 1
xrange = [-1500 1500]; % X range of measurements
yrange = [-1500 1500]; % Y range of measurements
A = (xrange(2)-xrange(1)) * (yrange(2)-yrange(1)); % Surveillance region.
c1_z = 1.0/A;   % lambda_c1 = lambda_1/A; % average number of clutter returns per unit area
clutter_intensity1 = @(lambda_1,z_cartesian) lambda_1 * c1_z; % Generate clutter function.

% Clutter intensity due to background for target type 2
xrange = [-1500 1500]; % X range of measurements
yrange = [-1500 1500]; % Y range of measurements

A = (xrange(2)-xrange(1)) * (yrange(2)-yrange(1)); % Surveillance region.
c2_z = 1.0/A;   % lambda_c2 = lambda_2/A; % average number of clutter returns per unit area
clutter_intensity2 = @(lambda_2,z_cartesian) lambda_2 * c2_z; % Generate clutter function.

A = (xrange(2)-xrange(1)) * (yrange(2)-yrange(1)); % Surveillance region.
c3_z = 1.0/A;   % lambda_c2 = lambda_2/A; % average number of clutter returns per unit area
clutter_intensity3 = @(lambda_3,z_cartesian) lambda_3 * c3_z; % Generate clutter function.

A = (xrange(2)-xrange(1)) * (yrange(2)-yrange(1)); % Surveillance region.
c4_z = 1.0/A;   % lambda_c2 = lambda_2/A; % average number of clutter returns per unit area
clutter_intensity4 = @(lambda_4,z_cartesian) lambda_4 * c4_z; % Generate clutter function.

%% Number of targets at each stage
numTargets_Jk_minus_11 = 0; % Number of targets of type 1, previous. J_k-1. Set in end of GM_PHD_Prune
numTargets_Jk_k_minus_11 = 0; % Number of targets given previous, J_k|k-1. Set at end of Step 2 (prediction of existing targets)
numTargets_Jk1 = 0; % Number of targets after update. J_k. Set at end of step 4 (update)

numTargets_Jk_minus_12 = 0; % Number of targets of type 2, previous. J_k-1. Set in end of GM_PHD_Prune
numTargets_Jk_k_minus_12 = 0; % Number of targets given previous, J_k|k-1. Set at end of Step 2 (prediction of existing targets)
numTargets_Jk2 = 0; % Number of targets after update. J_k. Set at end of step 4 (update)

numTargets_Jk_minus_13 = 0; % Number of targets of type 3, previous. J_k-1. Set in end of GM_PHD_Prune
numTargets_Jk_k_minus_13 = 0; % Number of targets given previous, J_k|k-1. Set at end of Step 2 (prediction of existing targets)
numTargets_Jk3 = 0; % Number of targets after update. J_k. Set at end of step 4 (update)

numTargets_Jk_minus_14 = 0; % Number of targets of type 4, previous. J_k-1. Set in end of GM_PHD_Prune
numTargets_Jk_k_minus_14 = 0; % Number of targets given previous, J_k|k-1. Set at end of Step 2 (prediction of existing targets)
numTargets_Jk4 = 0; % Number of targets after update. J_k. Set at end of step 4 (update)

%% Initial intensities: The previous iteration's mean/weight/covariance. Set just after pruning in
% N_type_GM_PHD_PruneAndMerge.m. Used as input for the prediction step. 
prunedIntensity.w1 = {}; % Weights from previous iteration
prunedIntensity.m1 = {}; % Means from previous iteration
prunedIntensity.p1 = {}; % Covariances from previous iteration
prunedIntensity.w2 = {}; % Weights from previous iteration
prunedIntensity.m2 = {}; % Means from previous iteration
prunedIntensity.p2 = {}; % Covariances from previous iteration
prunedIntensity.w3 = {}; % Weights from previous iteration
prunedIntensity.m3 = {}; % Means from previous iteration
prunedIntensity.p3 = {}; % Covariances from previous iteration
prunedIntensity.w4 = {}; % Weights from previous iteration
prunedIntensity.m4 = {}; % Means from previous iteration
prunedIntensity.p4 = {}; % Covariances from previous iteration

%% Prune and Merge parameters
% Merge and prune constants - these numbers have a HUGE impact on the performance and unfortunately need to be
% manually tuned if we make changes. 
T = 10^-5; % Pruning weight threshold. Values of the weights need to be above this value to be considered a target unles it will be deleted (pruned) immediately.
mergeThresholdU = 4; % Merge distance threshold. Points with Mahalanobis distance of less than this value between them will be merged.
weightThresholdToBeExtracted = 0.5; % State extraction weight threshold, values of the weights need to be above this value to be considered a 'real' target i.e. to be extracted.
maxGaussiansJ = 100; % Maximum number of Gaussians after pruning i.e. if there are still too many components after 
                     % merging, only Jmax = 100 components with the heighest weights are saved. NOT USED in this implementation.

%% Store estimated / extracted states
% We store the history of all points X_k (extracted states) for plotting
% purposes. This is updated in the end of the N_type_GM_PHD_ExtractStates
X_k_history1 = [];
X_k_history2 = [];
X_k_history3 = [];
X_k_history4 = [];

%% Calculate performance metric
% The order_p and cutoff_c control filter performance.Parameters given here 
% by default are not particularly well tuned for this problem but they work alright.
order_p = 1; % The order determines how punishing the metric calculation is to larger errors; as p increases, outliers are more heavily penalised
cutoff_c = 100; % Cutoff determines the maximum error for a point.
metric_history = []; % History of the OSPA performance metric

% Cardinality
numTarGt = [];
numTarEstByWeight = [];
numTarEstByStatesSize = [];

%%




