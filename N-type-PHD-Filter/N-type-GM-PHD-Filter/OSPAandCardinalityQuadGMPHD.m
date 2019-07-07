% Loading saved files and plotting the evaluations of OSPA and Cardinality
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

clc;
clear;
close all;

% % Cardinality
% load('CardinalityByStatesSizeQuad03.mat');
% load('CardinalityByStatesSizeIndep03.mat');
% load('CardinalityByWeightsQuad03.mat');
% load('CardinalityByWeightsIndep03.mat');
% load('CardinalityGT03.mat');

load('CardinalityByStatesSizeQuad06.mat');
load('CardinalityByStatesSizeIndep06.mat');
load('CardinalityByWeightsQuad06.mat');
load('CardinalityByWeightsIndep06.mat');
load('CardinalityGT06.mat');

% load('CardinalityByStatesSizeQuad09.mat');
% load('CardinalityByStatesSizeIndep09.mat');
% load('CardinalityByWeightsQuad09.mat');
% load('CardinalityByWeightsIndep09.mat');
% load('CardinalityGT09.mat');

% OSPA metric
% load('metric_historyQuad03.mat');
% load('metric_historyIndep03.mat');

load('metric_historyQuad06.mat');
load('metric_historyIndep06.mat');

% load('metric_historyQuad09.mat');
% load('metric_historyIndep09.mat');

cutoff_c = 100;

figure(1), hold on
plot(numTarGt, 'x-r'); 
plot(numTarEstByStatesSizeQuad , 'x-g'); 
plot(numTarEstByStatesSizeIndep, 'x-b');  
xlabel('Simulation (time) step');
ylabel('Cardinality');
title('Cardinality: ground truth (red), quad-GM-PHD filters(blue), 4 independent GM-PHD filters (green)');
hleg1 = legend('Ground truth cardinality', 'Estimated cardinality by quad-GM-PHD filter','Estimated cardinality by 4 GM-PHD filters');

hold off

figure(2),hold on
axis([0 120 0 cutoff_c]);
plot(metric_historyQuad, 'x-g'); 
plot(metric_historyIndep, 'x-b');
xlabel('Simulation (time) step');
ylabel('OSPA error metric (higher is worse)');
title('OSPA error metric: quad-GM-PHD filters(blue), 4 independent GM-PHD filters (green)');
hleg1 = legend('OSPA error using quad-GM-PHD filter','OSPA error using 4 GM-PHD filters');
hold off
        
   
OSPAtotalQuad = sum(metric_historyQuad)/length(metric_historyQuad)
OSPAtotalIndep = sum(metric_historyIndep)/length(metric_historyIndep)


%% Plot OSPA values of both quad GM-PHD filter and 4 independent GM-PHD filters
% y1 = [32.50, 33.15, 33.32, 34.22]; % old
% y2 = [32.50, 35.57, 48.38, 56.53];

y1 = [28.40, 28.70, 28.81, 29.17];
y2 = [28.40, 32.18, 46.47, 55.86];
x =  [0, 0.3, 0.6, 0.9];
figure, plot(x,y1,'-r'), hold on
plot(x,y2,'-b');
xlabel('Confusion detection probabilities');
ylabel('OSPA error (higher is worse)');
title('OSPA error: quad-GM-PHD filter(red), 4 independent GM-PHD filters (blue)');
hleg1 = legend('OSPA error using quad-GM-PHD filter','OSPA error using 4 GM-PHD filters');
hold off
        
%%



