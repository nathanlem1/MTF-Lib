%% ------------------------------------------------------------------------
% N_type_GM_PHD_Simulate_Measurements.m
% -------------------------------------------------------------------------
%
% This file generates simulated measurement data for the simulation
% There will be gaussian noise on the measurement and Poisson-distributed clutter
% in the environment. 
%
% If you want to use this N-type GM-PHD filter implementation for another problem, 
% you will need to replace this script with another one that populates Z, 
% for example from object detection algorithms for visual target tracking.
% Z is used in the update code, zTrue and simMeasurementHistory are used in
% N-type_GM_PHD_Simulate_Plot).
% Note that measurements might not be obtained, for example, if the target is not
% detected.
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

%%

s = sprintf('Step Sim: Simulating measurements.');
disp(s);

Q_size = size(model.Q, 1);  % 4 for CV 
R_size = size(model.R11, 1); % 2 for CV; R_size is the same for all target types 

if simulateTargetStateAndMeasurement == 1   
    %=============== Simulate target movement=============================
    if k == 1
        simTarget1State = simTarget1Start;
    end
    if k > 1 && k <= simTargetEndTime1
        simTarget1State = model.F * simTarget1State + sqrt(model.Q)* randn(Q_size,1); % ground tuth is X_k = F*X_k_1 + W_k_1 -> linear equation for state generation.
    end
    simTarget2State = model.F * simTarget2State + sqrt(model.Q)* randn(Q_size,1); 
    simTarget3State = model.F * simTarget3State + sqrt(model.Q)* randn(Q_size,1); 
    if k == simTargetStartTime
        simTarget4State = simTarget4Start;
    end
    if k > simTargetStartTime
        simTarget4State = model.F * simTarget4State + sqrt(model.Q)* randn(Q_size,1);
    end
    simTarget5State = model.F * simTarget5State + sqrt(model.Q)* randn(Q_size,1); 
    simTarget6State = model.F * simTarget6State + sqrt(model.Q)* randn(Q_size,1);
    simTarget7State = model.F * simTarget7State + sqrt(model.Q)* randn(Q_size,1); 
    if k == simTargetStartTime
        simTarget8State = simTarget8Start;
    end
    if k > simTargetStartTime
        simTarget8State = model.F * simTarget8State + sqrt(model.Q)* randn(Q_size,1);
    end
    simTarget9State = model.F * simTarget9State + sqrt(model.Q)* randn(Q_size,1); 
    simTarget10State = model.F * simTarget10State + sqrt(model.Q)* randn(Q_size,1); 
    simTarget11State = model.F * simTarget11State + sqrt(model.Q)* randn(Q_size,1);
    if k == simTargetStartTime
        simTarget12State = simTarget12Start;
    end 
    if k > simTargetStartTime && k <= simTargetEndTime2
        simTarget12State = model.F * simTarget12State + sqrt(model.Q)* randn(Q_size,1);
    end
    simTarget13State = model.F * simTarget13State + sqrt(model.Q)* randn(Q_size,1);
    simTarget14State = model.F * simTarget14State + sqrt(model.Q)* randn(Q_size,1); 
    simTarget15State = model.F * simTarget15State + sqrt(model.Q)* randn(Q_size,1);
    if k == simTargetStartTime
        simTarget16State = simTarget16Start;
    end 
    if k > simTargetStartTime && k <= simTargetEndTime2
        simTarget16State = model.F * simTarget16State + sqrt(model.Q)* randn(Q_size,1);
    end
    
    % Save target movement for plotting
    simTarget1History = [simTarget1History, simTarget1State];
    simTarget2History = [simTarget2History, simTarget2State];
    simTarget3History = [simTarget3History, simTarget3State];
    simTarget4History = [simTarget4History, simTarget4State];
    simTarget5History = [simTarget5History, simTarget5State];
    simTarget6History = [simTarget6History, simTarget6State];
    simTarget7History = [simTarget7History, simTarget7State];
    simTarget8History = [simTarget8History, simTarget8State];
    simTarget9History = [simTarget9History, simTarget9State];
    simTarget10History = [simTarget10History, simTarget10State];
    simTarget11History = [simTarget11History, simTarget11State];
    simTarget12History = [simTarget12History, simTarget12State];
    simTarget13History = [simTarget13History, simTarget13State];
    simTarget14History = [simTarget14History, simTarget14State];
    simTarget15History = [simTarget15History, simTarget15State];
    simTarget16History = [simTarget16History, simTarget16State];
else    
    % ====== Use already generated trajectories =========================
    % load('simulatedData/simMeasurementHistory03.mat');
%     load('simulatedData/simMeasurementHistory06.mat');
    % load('simulatedData/simMeasurementHistory09.mat');
    load('simulatedData/simTarget1History4.mat');
    load('simulatedData/simTarget2History4.mat');
    load('simulatedData/simTarget3History4.mat');
    load('simulatedData/simTarget4History4.mat');
    load('simulatedData/simTarget5History4.mat');
    load('simulatedData/simTarget6History4.mat');
    load('simulatedData/simTarget7History4.mat');
    load('simulatedData/simTarget8History4.mat');
    load('simulatedData/simTarget9History4.mat');
    load('simulatedData/simTarget10History4.mat');
    load('simulatedData/simTarget11History4.mat');
    load('simulatedData/simTarget12History4.mat');
    load('simulatedData/simTarget13History4.mat');
    load('simulatedData/simTarget14History4.mat');
    load('simulatedData/simTarget15History4.mat');
    load('simulatedData/simTarget16History4.mat');
    
    if k >= 1 && k < simTargetEndTime1
        simTarget1State = simTarget1History(:,k+1);
    end
    simTarget2State = simTarget2History(:,k+1);
    simTarget3State = simTarget3History(:,k+1);
    if k >= simTargetStartTime
        simTarget4State = simTarget4History(:,k+1-simTargetStartTime);
    end
    simTarget5State = simTarget5History(:,k+1);
    simTarget6State = simTarget6History(:,k+1);
    simTarget7State = simTarget7History(:,k+1);
    if k >= simTargetStartTime
        simTarget8State = simTarget8History(:,k+1-simTargetStartTime);
    end
    simTarget9State = simTarget9History(:,k+1);
    simTarget10State = simTarget10History(:,k+1);
    simTarget11State = simTarget11History(:,k+1);
    if k >= simTargetStartTime && k <= simTargetEndTime2
        simTarget12State = simTarget12History(:,k+1-simTargetStartTime);
    end
    simTarget13State = simTarget13History(:,k+1);
    simTarget14State = simTarget14History(:,k+1);
    simTarget15State = simTarget15History(:,k+1);
    if k >= simTargetStartTime && k <= simTargetEndTime2
        simTarget16State = simTarget16History(:,k+1-simTargetStartTime);
    end
end

% First, we generate some clutter in the environment.
% For target type 1
lambda_1 = poissrnd(3); % Number of clutter points, clutter is Poisson distributed though you can assume constant clutter measurements, lambda_1 = 3;  
clutter1 = zeros(2, lambda_1); % The observations are of the form [x; y]
for i = 1:lambda_1
    clutter1X = rand * (xrange(2) - xrange(1)) + xrange(1); % Random number between xrange(1) and xrange(2), uniformly distributed.
    clutter1Y = rand * (yrange(2) - yrange(1)) + yrange(1); % Random number between yrange(1) and yrange(2), uniformly distributed.
    
    clutter1(1,i) = clutter1X;
    clutter1(2,i) = clutter1Y;
end

% For target type 2
lambda_2 = poissrnd(3); % Number of clutter points, clutter is Poisson distributed though you can assume constant clutter measurements, lambda_2 = 3;  
clutter2 = zeros(2, lambda_2); % The observations are of the form [x; y]
for i = 1:lambda_2
    clutter2X = rand * (xrange(2) - xrange(1)) + xrange(1); % Random number between xrange(1) and xrange(2), uniformly distributed.
    clutter2Y = rand * (yrange(2) - yrange(1)) + yrange(1); % Random number between yrange(1) and yrange(2), uniformly distributed.
    
    clutter2(1,i) = clutter2X;
    clutter2(2,i) = clutter2Y;
end

% For target type 3
lambda_3 = poissrnd(3); % Number of clutter points, clutter is Poisson distributed though you can assume constant clutter measurements, lambda_2 = 3;  
clutter3 = zeros(2, lambda_3); % The observations are of the form [x; y]
for i = 1:lambda_3
    clutter3X = rand * (xrange(2) - xrange(1)) + xrange(1); % Random number between xrange(1) and xrange(2), uniformly distributed.
    clutter3Y = rand * (yrange(2) - yrange(1)) + yrange(1); % Random number between yrange(1) and yrange(2), uniformly distributed.
    
    clutter3(1,i) = clutter3X;
    clutter3(2,i) = clutter3Y;
end

% For target type 4
lambda_4 = poissrnd(3); % Number of clutter points, clutter is Poisson distributed though you can assume constant clutter measurements, lambda_2 = 3;  
clutter4 = zeros(2, lambda_4); % The observations are of the form [x; y]
for i = 1:lambda_4
    clutter4X = rand * (xrange(2) - xrange(1)) + xrange(1); % Random number between xrange(1) and xrange(2), uniformly distributed.
    clutter4Y = rand * (yrange(2) - yrange(1)) + yrange(1); % Random number between yrange(1) and yrange(2), uniformly distributed.
    
    clutter4(1,i) = clutter4X;
    clutter4(2,i) = clutter4Y;
end

% We are not guaranteed to detect the target - there is only a probability
detect1 = rand; % Uniform distribution.
detect2 = rand;
detect3 = rand;
detect4 = rand;
detect5 = rand;
detect6 = rand;
detect7 = rand;
detect8 = rand;
detect9 = rand;
detect10 = rand;
detect11 = rand;
detect12 = rand;
detect13 = rand;
detect14 = rand;
detect15 = rand;
detect16 = rand;
 
% A number in the interval (0,prob_detection] will be observed with probability p = prob_detection.
% Detection of target type 1 by detector 1 ( targets 1, 2, 3 and 4) prob_detection11
if(detect1 > prob_detection11)
    meas1 = [];
else
    if k < 1 || k > simTargetEndTime1
        meas1 = [];
    else
        meas1 = model.H11*simTarget1State + sqrt(model.R11) * randn(R_size,1) * noiseScaler1; % simTarget1State(1) is the mean, meas1 is the
                                                                 % true measurement which generated from normal distribution 
    end   
end
% measZ1 = GenerateObservation(model.H, simTarget1State, prob_detection, model.sigma_r, detect1); 
% if ~isempty(measZ1)
%     measX1 = measZ1(1);
%     measY1 = measZ1(2);
% else
%     measX1 = [];
%     measY1 = [];
% end

if(detect2 > prob_detection11)
    meas2 = [];
else
    meas2 = model.H11*simTarget2State + sqrt(model.R11) * randn(R_size,1) * noiseScaler1;  % Z_k = H_k*X_k + V_k, look at GenerateObservation function.
end
if(detect3 > prob_detection11)
    meas3 = [];
else
    meas3 = model.H11*simTarget3State + sqrt(model.R11) * randn(R_size,1) * noiseScaler1;  % Z_k = H_k*X_k + V_k, look at GenerateObservation function.
end
if(detect4 > prob_detection11)
    meas4 = [];
else
    if k < simTargetStartTime
        meas4 = [];
    else
        meas4 = model.H11*simTarget4State + sqrt(model.R11) * randn(R_size,1) * noiseScaler1;  % Z_k = H_k*X_k + V_k, look at GenerateObservation function.
    end
end
% Detection of target type 2 by detector 1 ( target 5) using prob_detection12
if(detect5 > prob_detection12)
    meas_s1t5 = [];
else
    meas_s1t5 = model.H12*simTarget5State + sqrt(model.R12) * randn(R_size,1) * noiseScaler1;
end
% Detection of target type 3 by detector 1 ( target 9) using prob_detection13
if(detect9 > prob_detection13)
    meas_s1t9 = [];
else
    meas_s1t9 = model.H13*simTarget9State + sqrt(model.R13) * randn(R_size,1) * noiseScaler1;
end
% Detection of target type 4 by detector 1 ( target 13) using prob_detection14
if(detect13 > prob_detection14)
    meas_s1t13 = [];
else
    meas_s1t13 = model.H14*simTarget13State + sqrt(model.R14) * randn(R_size,1) * noiseScaler1;
end

% Detection of target type 2 by detector 2 ( targets 5, 6, 7 and 8) prob_detection22
if(detect5 > prob_detection22)
    meas5 = [];
else
    meas5 = model.H22*simTarget5State + sqrt(model.R22) * randn(R_size,1) * noiseScaler2; % simTarget1State(1) is the mean, measX1 is the
                                                                 % true measurement which generated from normal distribution 

if(detect6 > prob_detection22)
    meas6 = [];
else
    meas6 = model.H22*simTarget6State + sqrt(model.R22) * randn(R_size,1) * noiseScaler2;
end
end
if(detect7 > prob_detection22)
    meas7 = [];
else
    meas7 = model.H22*simTarget7State + sqrt(model.R22) * randn(R_size,1) * noiseScaler2;
end
if(detect8 > prob_detection22)
    meas8 = [];
else
    if k < simTargetStartTime
        meas8 = [];
    else
        meas8 = model.H22*simTarget8State + sqrt(model.R22) * randn(R_size,1) * noiseScaler2;
    end
end
% Detection of target type 1 by detector 2 ( target 1) using prob_detection21
if(detect1 > prob_detection21)
    meas_s2t1 = [];
else
    if k < 1 || k > simTargetEndTime1
        meas_s2t1 = [];
    else
        meas_s2t1 = model.H21*simTarget1State + sqrt(model.R21) * randn(R_size,1) * noiseScaler2;  % Z_k = H_k*X_k + V_k, look at GenerateObservation function.
    end  
end
% Detection of target type 3 by detector 2 ( target 10) using prob_detection23
if(detect10 > prob_detection23)
    meas_s2t10 = [];
else
    meas_s2t10 = model.H23*simTarget10State + sqrt(model.R23) * randn(R_size,1) * noiseScaler2;  % Z_k = H_k*X_k + V_k, look at GenerateObservation function.
end
% Detection of target type 4 by detector 2 ( target 14) using prob_detection24
if(detect14 > prob_detection24)
    meas_s2t14 = [];
else
    meas_s2t14 = model.H24*simTarget14State + sqrt(model.R24) * randn(R_size,1) * noiseScaler2;  % Z_k = H_k*X_k + V_k, look at GenerateObservation function.
end

% Detection of target type 3 by detector 3 ( targets 9, 10, 11 and 12) prob_detection33
if(detect9 > prob_detection33)
    meas9 = [];
else
    meas9 = model.H33*simTarget9State + sqrt(model.R33) * randn(R_size,1) * noiseScaler3;  
end
if(detect10 > prob_detection33)
    meas10 = [];
else
    meas10 = model.H33*simTarget10State + sqrt(model.R33) * randn(2,1) * noiseScaler3;
end
if(detect11 > prob_detection33)
    meas11 = [];
else
    meas11 = model.H33*simTarget11State + sqrt(model.R33) * randn(R_size,1) * noiseScaler3;
end
if(detect12 > prob_detection33)
    meas12 = [];
else
    if k < simTargetStartTime || k > simTargetEndTime2
        meas12 = [];
    else
        meas12 = model.H33*simTarget12State + sqrt(model.R33) * randn(R_size,1) * noiseScaler3;
    end
end
% Detection of target type 1 by detector 3 ( target 2) using prob_detection31
if(detect2 > prob_detection31)
    meas_s3t2 = [];
else
    meas_s3t2 = model.H31*simTarget2State + sqrt(model.R31) * randn(R_size,1) * noiseScaler3;  % Z_k = H_k*X_k + V_k, look at GenerateObservation function.
end
% Detection of target type 2 by detector 3 ( target 6) using prob_detection32
if(detect6 > prob_detection32)
    meas_s3t6 = [];
else
    meas_s3t6 = model.H32*simTarget6State + sqrt(model.R32) * randn(R_size,1) * noiseScaler3;  % Z_k = H_k*X_k + V_k, look at GenerateObservation function.
end
% Detection of target type 4 by detector 3 ( target 15) using prob_detection34
if(detect15 > prob_detection34)
    meas_s3t15 = [];
else
    meas_s3t15 = model.H34*simTarget15State + sqrt(model.R34) * randn(R_size,1) * noiseScaler3;  % Z_k = H_k*X_k + V_k, look at GenerateObservation function.
end

% Detection of target type 4 by detector 4 ( targets 13, 14, 15 and 16) prob_detection44
if(detect13 > prob_detection44)
    meas13 = [];
else
    meas13 = model.H44*simTarget13State + sqrt(model.R44) * randn(R_size,1) * noiseScaler4;  
end
if(detect14 > prob_detection44)
    meas14 = [];
else
    meas14 = model.H44*simTarget14State + sqrt(model.R44) * randn(R_size,1) * noiseScaler4;
end
if(detect15 > prob_detection44)
    meas15 = [];
else
    meas15 = model.H44*simTarget15State + sqrt(model.R44) * randn(R_size,1) * noiseScaler4;
end
if(detect16 > prob_detection44)
    meas16 = [];
else
    if k < simTargetStartTime || k > simTargetEndTime2
        meas16 = [];
    else
        meas16 = model.H44*simTarget16State + sqrt(model.R44) * randn(R_size,1) * noiseScaler4;
    end
end
% Detection of target type 1 by detector 4 ( target 3) using prob_detection41
if(detect3 > prob_detection41)
    meas_s4t3 = [];
else
    meas_s4t3 = model.H41*simTarget3State + sqrt(model.R41) * randn(R_size,1) * noiseScaler4;  % Z_k = H_k*X_k + V_k, look at GenerateObservation function.
end
% Detection of target type 2 by detector 4 ( target 7) using prob_detection42
if(detect7 > prob_detection42)
    meas_s4t7 = [];
else
    meas_s4t7 = model.H42*simTarget7State + sqrt(model.R42) * randn(R_size,1) * noiseScaler4;  % Z_k = H_k*X_k + V_k, look at GenerateObservation function.
end
% Detection of target type 3 by detector 4 ( target 11) using prob_detection43
if(detect11 > prob_detection43)
    meas_s4t11 = [];
else
    meas_s4t11 = model.H43*simTarget11State + sqrt(model.R43) * randn(R_size,1) * noiseScaler4;  % Z_k = H_k*X_k + V_k, look at GenerateObservation function.
end

% Generate true measurement (measurement without clutter) 
% (NB. It is a measurement/observation, not the ground truth (mean))
% Expected realistic scenario: 4 target types with confused detections
Z1 = [ meas1 meas2 meas3 meas4 meas_s1t5 meas_s1t9 meas_s1t13 ];
Z2 = [ meas5 meas6 meas7 meas8 meas_s2t1 meas_s2t10 meas_s2t14 ];
Z3 = [ meas9 meas10 meas11 meas12 meas_s3t2 meas_s3t6 meas_s3t15 ];
Z4 = [ meas13 meas14 meas15 meas16 meas_s4t3 meas_s4t7 meas_s4t11 ];

zTrue1 = Z1;  % Store for plotting
zTrue2 = Z2;  % Store for plotting
zTrue3 = Z3;  % Store for plotting
zTrue4 = Z4;  % Store for plotting
simTrueMeasurementHistoryZ1 = [simTrueMeasurementHistoryZ1 zTrue1]; 
simTrueMeasurementHistoryZ2 = [simTrueMeasurementHistoryZ2 zTrue2]; 
simTrueMeasurementHistoryZ3 = [simTrueMeasurementHistoryZ3 zTrue3]; 
simTrueMeasurementHistoryZ4 = [simTrueMeasurementHistoryZ4 zTrue4]; 

% Append clutter (false alarms) to the true measurement to generate cluttered measurement
Z1 = [Z1, clutter1]; % used for the update step
Z2 = [Z2, clutter2]; % used for the update step
Z3 = [Z3, clutter3]; % used for the update step
Z4 = [Z4, clutter4]; % used for the update step

% Store history, should be commented out when already saved measurements are
% used
simMeasurementHistory.Z1{k} =  Z1;
simMeasurementHistory.Z2{k} =  Z2;
simMeasurementHistory.Z3{k} =  Z3;
simMeasurementHistory.Z4{k} =  Z4;

%%


