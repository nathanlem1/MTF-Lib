%% ------------------------------------------------------------------------
% N_type_GM_PHD_Predict_Existing.m
% -------------------------------------------------------------------------
%
% This file performs prediction for existing targets.
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

%%

s = sprintf('Step 3: Prediction for existing targets.');
disp(s);

% For target type 1
for j = 1:numTargets_Jk_minus_11
    i1 = i1 + 1;
    predictedIntensity.w1{i1} = prob_survival1 * prunedIntensity.w1{j};
    predictedIntensity.m1{i1} = model.F1 * prunedIntensity.m1{j}; % Assume constant velocity.
    beforePredictionIntensity.m1{i1} = prunedIntensity.m1{j};  % store the position of targets before prediction
    predictedIntensity.p1{i1} = model.Q1 + model.F1 * prunedIntensity.p1{j} * model.F1';
end

numTargets_Jk_k_minus_11 = i1;  % numTargets_Jk_k_minus_11 = numTargets_Jk_minus_11 + numBirthedTargets1;

% For target type 2
for j = 1:numTargets_Jk_minus_12
    i2 = i2 + 1;
    predictedIntensity.w2{i2} = prob_survival2 * prunedIntensity.w2{j};
    predictedIntensity.m2{i2} = model.F2 * prunedIntensity.m2{j}; % Assume constant velocity.
    beforePredictionIntensity.m2{i2} = prunedIntensity.m2{j};  % store the position of targets before prediction
    predictedIntensity.p2{i2} = model.Q2 + model.F2 * prunedIntensity.p2{j} * model.F2';
end

numTargets_Jk_k_minus_12 = i2;  % numTargets_Jk_k_minus_12 = numTargets_Jk_minus_12 + numBirthedTargets2;

% For target type 3
for j = 1:numTargets_Jk_minus_13
    i3 = i3 + 1;
    predictedIntensity.w3{i3} = prob_survival3 * prunedIntensity.w3{j};
    predictedIntensity.m3{i3} = model.F3 * prunedIntensity.m3{j}; % Assume constant velocity.
    beforePredictionIntensity.m3{i3} = prunedIntensity.m3{j};  % store the position of targets before prediction
    predictedIntensity.p3{i3} = model.Q3 + model.F3 * prunedIntensity.p3{j} * model.F3';
end

numTargets_Jk_k_minus_13 = i3;  % numTargets_Jk_k_minus_13 = numTargets_Jk_minus_13 + numBirthedTargets3;

% For target type 4
for j = 1:numTargets_Jk_minus_14
    i4 = i4 + 1;
    predictedIntensity.w4{i4} = prob_survival4 * prunedIntensity.w4{j};
    predictedIntensity.m4{i4} = model.F4 * prunedIntensity.m4{j}; % Assume constant velocity.
    beforePredictionIntensity.m4{i4} = prunedIntensity.m4{j};  % store the position of targets before prediction
    predictedIntensity.p4{i4} = model.Q4 + model.F4 * prunedIntensity.p4{j} * model.F4';
end

numTargets_Jk_k_minus_14 = i4;  % numTargets_Jk_k_minus_13 = numTargets_Jk_minus_13 + numBirthedTargets3;

N_predicted = sprintf('Number of predicted targets by summing weights = %0.0f', round(sum([predictedIntensity.w1{:}])) + round(sum([predictedIntensity.w2{:}])) + round(sum([predictedIntensity.w3{:}])) + round(sum([predictedIntensity.w4{:}])));
disp(N_predicted);

%%

