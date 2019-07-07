%============== N_type_GM_PHD_Filter ===================================

% The N-type GM-PHD Filter is proposed in Nathanael L. Baisa & Andrew Wallace's work:
% Nathanael L. Baisa, Andrew Wallace, "Multiple target, multiple type filtering in the RFS framework",
% Digital Signal Processing, Vol 89, June 2019, Pages 49-59

% This code is for Quad GM-PHD filter which is one example of the N-type GM-PHD
% filter when N = 4 and was used in the above mentioned paper.

% Reference for some estimation code examples is available in: 
% http://www-personal.acfr.usyd.edu.au/tbailey/software/

% The OSPA metric provides a nice way of analysing performance.

% A constant velocity (CV) motion model is used for all target types.

%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

%% Main loop - Quad-GM-PHD Filter
clc;
clear all;
close all;
tic

%Step 0: Initialisation
N_type_GM_PHD_Initialisation;
N_type_GM_PHD_Simulate_Initialise;

k = 0;
while (k < endTime) %k = timestep, endTime = 120
    k = k + 1;
    s = sprintf('======ITERATION %d======', k);
    disp(s);
        
    % Step Sim: Generate sensor Measurements
    % If you want to use this code with your own data or for a different problem,
    % replace this function with your own.
    
    N_type_GM_PHD_Simulate_Measurements;   % Linear N-type GM-PHD filter measurements are simulated direct observations [X; Y] of the target positions   
                                           % Linear N-type GM-PHD filter uses fixed matrices for prediction and update.

    % Step 1: Targets's birth 
    % it uses measurement driven birthing of targets, no spawning is considered as it is handled using birthing only.
    N_type_GM_PHD_Create_Birth; 
    % Step 2: Prediction for existing or persistent targets
    N_type_GM_PHD_Predict_Persistent;
    % Step 3: Construction of PHD update components
    N_type_GM_PHD_Construct_Update_Components;
    % Step 4: Update targets with measurements
    N_type_GM_PHD_Update;
    % Step 5: Prune and merge targets
    N_type_GM_PHD_PruneAndMerge;
    % Step 6: Estimate position of targets
    N_type_GM_PHD_ExtractStates;
    
    % Step Metric: Calculate performance metric
    N_type_GM_PHD_Calculate_Performance_Metric;
    
    % Step Plot: Generate graphs
    N_type_GM_PHD_Simulate_Plot; drawnow;
end
sum(metric_history)/endTime
(sum(numTarGt)- sum(numTarEstByStatesSize))/endTime

% % For Quad-GM-PHD filter
% numTarEstByWeightQuad = numTarEstByWeight;
% numTarEstByStatesSizeQuad  = numTarEstByStatesSize;
% metric_historyQuad = metric_history;
% save('CardinalityByWeightsQuad06', 'numTarEstByWeightQuad'); 
% save('CardinalityByStatesSizeQuad06', 'numTarEstByStatesSizeQuad');
% save('metric_historyQuad06', 'metric_historyQuad'); 
% save('CardinalityGT06', 'numTarGt'); 

% % For four independent GM-PHD filters
% numTarEstByWeightIndep = numTarEstByWeight;
% numTarEstByStatesSizeIndep  = numTarEstByStatesSize;
% metric_historyIndep = metric_history;
% save('CardinalityByWeightsIndep06', 'numTarEstByWeightIndep'); 
% save('CardinalityByStatesSizeIndep06', 'numTarEstByStatesSizeIndep');
% save('metric_historyIndep06', 'metric_historyIndep');

% % Save measurements
% save('simMeasurementHistory09', 'simMeasurementHistory');
% save('simMeasurementHistory06', 'simMeasurementHistory');
% save('simMeasurementHistory03', 'simMeasurementHistory');

% % Save trajectories
% save('simTarget1History4', 'simTarget1History');
% save('simTarget2History4', 'simTarget2History');
% save('simTarget3History4', 'simTarget3History');
% save('simTarget4History4', 'simTarget4History');
% save('simTarget5History4', 'simTarget5History');
% save('simTarget6History4', 'simTarget6History');
% save('simTarget7History4', 'simTarget7History');
% save('simTarget8History4', 'simTarget8History');
% save('simTarget9History4', 'simTarget9History');
% save('simTarget10History4', 'simTarget10History');
% save('simTarget11History4', 'simTarget11History');
% save('simTarget12History4', 'simTarget12History');
% save('simTarget13History4', 'simTarget13History');
% save('simTarget14History4', 'simTarget14History');
% save('simTarget15History4', 'simTarget15History');
% save('simTarget16History4', 'simTarget16History');

timePHD = toc 

% print(gcf,'-dpsc2','dualSim09.eps') % to save in .eps file format which is good for displaying images for paper publications
%%


