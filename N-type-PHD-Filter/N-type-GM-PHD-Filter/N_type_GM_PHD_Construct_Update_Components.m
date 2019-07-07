%% ------------------------------------------------------------------------
% N_type_GM_PHD_Construct_Update_Components.m
% -------------------------------------------------------------------------
%
% This file creates the components needed for performing a N-type GM_PHD 
% filter update for each target type using the measurements.
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

%%

s = sprintf('Step 4: Constructing update components for all targets, new and existing.');
disp(s);

% For target type 1
for j = 1:numTargets_Jk_k_minus_11
    
    constructUpdateIntensity.eta1{j} = model.H11*predictedIntensity.m1{j};
    constructUpdateIntensity.S1{j} = model.R11 + model.H11*predictedIntensity.p1{j}*model.H11';
    constructUpdateIntensity.K1{j} = predictedIntensity.p1{j}*model.H11'*inv(constructUpdateIntensity.S1{j});
    constructUpdateIntensity.p1{j} = predictedIntensity.p1{j} - constructUpdateIntensity.K1{j}*model.H11*predictedIntensity.p1{j};
end

% For target type 2
for j = 1:numTargets_Jk_k_minus_12
    
    constructUpdateIntensity.eta2{j} = model.H22*predictedIntensity.m2{j};
    constructUpdateIntensity.S2{j} = model.R22 + model.H22*predictedIntensity.p2{j}*model.H22';
    constructUpdateIntensity.K2{j} = predictedIntensity.p2{j}*model.H22'*inv(constructUpdateIntensity.S2{j});
    constructUpdateIntensity.p2{j} = predictedIntensity.p2{j} - constructUpdateIntensity.K2{j}*model.H22*predictedIntensity.p2{j};
end

% For target type 3
for j = 1:numTargets_Jk_k_minus_13
    
    constructUpdateIntensity.eta3{j} = model.H33*predictedIntensity.m3{j};
    constructUpdateIntensity.S3{j} = model.R33 + model.H33*predictedIntensity.p3{j}*model.H33';
    constructUpdateIntensity.K3{j} = predictedIntensity.p3{j}*model.H33'*inv(constructUpdateIntensity.S3{j});
    constructUpdateIntensity.p3{j} = predictedIntensity.p3{j} - constructUpdateIntensity.K3{j}*model.H33*predictedIntensity.p3{j};
end

% For target type 4
for j = 1:numTargets_Jk_k_minus_14
    
    constructUpdateIntensity.eta4{j} = model.H44*predictedIntensity.m4{j};
    constructUpdateIntensity.S4{j} = model.R44 + model.H44*predictedIntensity.p4{j}*model.H44';
    constructUpdateIntensity.K4{j} = predictedIntensity.p4{j}*model.H44'*inv(constructUpdateIntensity.S4{j});
    constructUpdateIntensity.p4{j} = predictedIntensity.p4{j} - constructUpdateIntensity.K4{j}*model.H44*predictedIntensity.p4{j};
end

%%




