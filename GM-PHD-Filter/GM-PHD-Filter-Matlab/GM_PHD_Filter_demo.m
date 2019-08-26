% ============== GM-PHD filter demo =======================================
% 
% This is the GM-PHD-Filter implementation of the algorithm in [1] with the
% assumption of measurement-driven birth of targets which has been used in 
% work [2]. In work [2] i.e. in visual tracking applications, extracting 
% states from the pruned intensity gives better result than extracting them
% from the updated intensity.
% 
% NB. While more restrictive than SMC approaches, Gaussian mixture 
% implementation is much more efficient.
% 
%- Refer for some estimation code examples in: 
% http://www-personal.acfr.usyd.edu.au/tbailey/software/
% 
% References:
% [1]. B.-N. Vo, W.-K. Ma, "The Gaussian Mixture Probability Hypothesis 
% Density Filter", IEEE Transactions on Signal Processing, Vol 54, No. 11, 
% November 2006, pp4091-4104
% [2]. Nathanael L. Baisa, "Online multi-object visual tracking using a 
% GM-PHD filter with deep appearance learning",22nd International Conference 
% on Information Fusion (FUSION), July, 2019.
% 
% -------------------------------------------------------------------------
%  Nathanael L. Baisa: nathanaellmss@gmail.com
%  Original:
%  Modified: March, 2014
% -------------------------------------------------------------------------

clc;
close all;
clear;

% Set model parameters 
model = Set_model; 

% Starting target states (positions , velocities)
m_start1 = [250, 250, 0, 0]';
m_start2 = [10,-10, 5, -5]';
m_start3 = [-250, -250, 0, 0]';

targetStates = [m_start1, m_start2, m_start3];

% Initialize the initial pruned intensity
prunedIntensity.w = {};
prunedIntensity.m = {};
prunedIntensity.p = {};

nScan = 100; % Duration of iterations (simulation)

% Plot figures
figure(1);
hold on;
axis([model.xrange(1) model.xrange(2) model.yrange(1) model.yrange(2)]);
figure(2);
hold on;
figure(3);
hold on;
ospa_all = [];
gtCardinality_all = [];
estimatedCardinality_all = [];

for i = 1:nScan
    disp(i)
    % Generate target states, actual observations, and clutter    
    targetStates = gen_newstates(targetStates, model);
    observations = gen_observations(targetStates, model);
    clutter = gen_clutter(model);
    Z_k = [observations, clutter];
    
    % Apply GM-PHD-Filter
    predictedIntensity = Predict(Z_k, prunedIntensity, model);
    updatedIntensity = Update(Z_k, predictedIntensity, model);
    prunedIntensity = PruneAndMerge(updatedIntensity, model);
    estimates = ExtractStates(prunedIntensity, model); %  extracting estimates from the pruned intensity this gives better result than extracting them from the updated intensity!
    estimatedStates = [];
    weights = [];
    for s =1:length(estimates.m)
        estimatedStates = [estimatedStates, estimates.m{s}];
        weights = [weights, estimates.w{s}];
    end   
    
    % Plot ground-truth states, true observations, estiamted states and clutter
    figure(1), plot_results(targetStates, observations, estimatedStates, clutter)
    
    % Compute OSPA metric, true cardinality, estimated cardinality and plot
    % them 
    ospa_metric = ospa_dist(estimatedStates, targetStates, model.cutoff_c, model.order_p); % Use Ba-Ngu Vo's implementation
    ospa_all = [ospa_all, ospa_metric];
    gtCardinality_all = [gtCardinality_all, size(targetStates,2)];
    estimatedCardinality_all = [estimatedCardinality_all, sum(weights)];
    figure(2), plot_ospa(ospa_all, model, nScan)
    figure(3), plot_cardinality(gtCardinality_all, estimatedCardinality_all);
    
end
%%


