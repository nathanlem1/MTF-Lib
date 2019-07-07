%% ------------------------------------------------------------------------
% N_type_GM_PHD_ExtactStates.m
% -------------------------------------------------------------------------
%
% This file estimates the positions of targets tracked by the N-type GM-PHD filter.
% We need to extract the likely target positions from the PHD (i.e. we need 
% to find the peaks of the PHD). This is actually fairly tricky. The naive 
% approach is to pull out targets with the highest weights, but this is FAR 
% from the best approach. A large covariance will pull down the peak size, 
% and when targets are close together or have high covariances, there can be
% superposition effects which shift the peak.

% This function implements a method which is pulling out every target with a 
% weight over weightThresholdToBeExtracted (defined in N_type_GM_PHD_Initialisation). 
% There is the option of repeatedly printing out targets with rounded weights 
% greater than 1. This will NOT change filter performance as the extracted 
% state estimate is not fed back into the filter.
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

%%

s = sprintf('Step 8: Estimate target states');
disp(s);

% Storing the extracted states, and their corresponding weights and
% covariances.
X_k1 = [];
X_k_w1 = [];
X_k_P1 = [];
X_k2 = [];
X_k_w2 = [];
X_k_P2 = [];
X_k3 = [];
X_k_w3 = [];
X_k_P3 = [];
X_k4 = [];
X_k_w4 = [];
X_k_P4 = [];

% Extracting states from state extraction step ==> this gives better performance 
% than extracting states from the update step!

% For target type 1
for i = 1:numTargets_J_pruned1
    if (prunedIntensity.w1{i} > weightThresholdToBeExtracted) 
       for j = 1:round(prunedIntensity.w1{i}) % If a target has a rounded weight greater than 1, output it multiple times.          
            X_k1 = [X_k1, prunedIntensity.m1{i}];
            X_k_w1 = [X_k_w1, prunedIntensity.w1{i}];
            X_k_P1 = [X_k_P1, prunedIntensity.p1{i}];
      end
    end
end

% For target type 2 in addition to the target type 1
for i = 1:numTargets_J_pruned2
    if (prunedIntensity.w2{i} > weightThresholdToBeExtracted) 
       for j = 1:round(prunedIntensity.w2{i}) % If a target has a rounded weight greater than 1, output it multiple times.          
            X_k2 = [X_k2, prunedIntensity.m2{i}];
            X_k_w2 = [X_k_w2, prunedIntensity.w2{i}];
            X_k_P2 = [X_k_P2, prunedIntensity.p2{i}];
      end
    end
end

% For target type 3 in addition to the target type 1 and 2
for i = 1:numTargets_J_pruned3
    if (prunedIntensity.w3{i} > weightThresholdToBeExtracted) 
       for j = 1:round(prunedIntensity.w3{i}) % If a target has a rounded weight greater than 1, output it multiple times.          
            X_k3 = [X_k3, prunedIntensity.m3{i}];
            X_k_w3 = [X_k_w3, prunedIntensity.w3{i}];
            X_k_P3 = [X_k_P3, prunedIntensity.p3{i}];
      end
    end
end


% For target type 4 in addition to the target type 1, 2 and 3
for i = 1:numTargets_J_pruned4
    if (prunedIntensity.w4{i} > weightThresholdToBeExtracted) 
       for j = 1:round(prunedIntensity.w4{i}) % If a target has a rounded weight greater than 1, output it multiple times.          
            X_k4 = [X_k4, prunedIntensity.m4{i}];
            X_k_w4 = [X_k_w4, prunedIntensity.w4{i}];
            X_k_P4 = [X_k_P4, prunedIntensity.p4{i}];
      end
    end
end


% % Extracting states from update step
% % For target type 1
% for i = 1:numTargets_Jk11
%     if (updateIntensity2.w1{i} > weightThresholdToBeExtracted) 
%        for j = 1:round(updateIntensity2.w1{i}) % If a target has a rounded weight greater than 1, output it multiple times.
%             X_k1 = [X_k1, updateIntensity2.m1{i}];
%             X_k_w1 = [X_k_w1, updateIntensity2.w1{i}];
%             X_k_P1 = [X_k_P1, updateIntensity2.p1{i}];
%        end
%     end
% end
% 
% % For target type 2
% for i = 1:numTargets_Jk22
%     if (updateIntensity2.w2{i} > weightThresholdToBeExtracted) 
%        for j = 1:round(updateIntensity2.w2{i}) % If a target has a rounded weight greater than 1, output it multiple times.
%             X_k2 = [X_k2, updateIntensity2.m2{i}];
%             X_k_w2 = [X_k_w2, updateIntensity2.w2{i}];
%             X_k_P2 = [X_k_P2, updateIntensity2.p2{i}];
%        end
%     end
% end
% 
% % For target type 3
% for i = 1:numTargets_Jk33
%     if (updateIntensity2.w3{i} > weightThresholdToBeExtracted) 
%        for j = 1:round(updateIntensity2.w3{i}) % If a target has a rounded weight greater than 1, output it multiple times.
%             X_k3 = [X_k3, updateIntensity2.m3{i}];
%             X_k_w3 = [X_k_w3, updateIntensity2.w3{i}];
%             X_k_P3 = [X_k_P3, updateIntensity2.p3{i}];
%        end
%     end
% end
% 
% % For target type 4
% for i = 1:numTargets_Jk44
%     if (updateIntensity2.w4{i} > weightThresholdToBeExtracted) 
%        for j = 1:round(updateIntensity2.w4{i}) % If a target has a rounded weight greater than 1, output it multiple times.
%             X_k4 = [X_k4, updateIntensity2.m4{i}];
%             X_k_w4 = [X_k_w4, updateIntensity2.w4{i}];
%             X_k_P4 = [X_k_P4, updateIntensity2.p4{i}];
%        end
%     end
% end

% Concatenate the results of the N-type GM-PHD filter outputs
X_k = [X_k1 X_k2 X_k3 X_k4];
X_k_w = [X_k_w1 X_k_w2 X_k_w3 X_k_w4];
X_k_P = [X_k_P1 X_k_P2 X_k_P3 X_k_P4];

% Store history for plotting.
X_k_history1 = [X_k_history1, X_k1];
X_k_history2 = [X_k_history2, X_k2];
X_k_history3 = [X_k_history3, X_k3];
X_k_history4 = [X_k_history4, X_k4];

CardinalityByWeight = round(sum(X_k_w));  % Compute the cardinality by summing up the weights.
CardinalityByStatesSize = size(X_k, 2);  % Get Cardinality from the number of extacted states.
cstatesByWt = sprintf('Number of extracted targets by summing weights = %0.0f', CardinalityByWeight);
disp(cstatesByWt);
cstatesByNo = sprintf('Number of extracted targets by counting states = %0.0f', CardinalityByStatesSize);
disp(cstatesByNo);

%%





