%% ------------------------------------------------------------------------
% N-type GM_PHD_PruneAndMerge.m
% -------------------------------------------------------------------------
%
% This file performs pruning and merging after the N-type GM-PHD filter update.
% The PHD update creates a combinatorial explosion in the number of targets.
% We prune the low-weight ones and merge the close-together ones. The weight
% threshold T and distance threshold mergeThresholdU that control this 
% process are set in N_type_GM_PHD_Initialisation.m.
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

%%

s = sprintf('Step 7: Prune and merge targets.');
disp(s);

% Make sure that the previous targets states not be repeated in the state 
% extraction again!
prunedIntensity.w1 = {};
prunedIntensity.m1 = {};
prunedIntensity.p1 = {};
prunedIntensity.w2 = {};
prunedIntensity.m2 = {};
prunedIntensity.p2 = {};
prunedIntensity.w3 = {};
prunedIntensity.m3 = {};
prunedIntensity.p3 = {};
prunedIntensity.w4 = {};
prunedIntensity.m4 = {};
prunedIntensity.p4 = {};

%% For the first target type

%% Prune out the low-weighted targets
I1 = find([updateIntensity.w1{:}] >= T); % Find targets with high enough weights, I is indices of truncated weights

%% Merge the close-together targets
l1 = 0; % counts number of features

while isempty(I1) == false  % We delete from I as we merge
    
    l1 = l1 + 1;
    %Find j, which is i corresponding to highest updateIntensity.w for all i in I
    highWeights = updateIntensity.w1{I1};
    [maxW, j] = max(highWeights);
    j = j(1); % In case of two targets with equal weight
    % j is an index of highWeights (i.e. a position in I)
    % We want the index in updateIntensity.w
    j = I1(j);
    
    % Find all points with Mahalanobis distance less than U from point updateIntensity.m(j)
    L1 = [];  % A vector of indices of merged Gaussians.
    for iterateI = 1:length(I1)
        thisI = I1(iterateI);
        delta_m = updateIntensity.m1{thisI} - updateIntensity.m1{j};
        %mahal_dist = delta_m' * (updateIntensity.p{thisI} \ delta_m); % Invert covariance via left division
        mahal_dist = delta_m' * inv(updateIntensity.p1{thisI}) * delta_m; % inv(A)*b == A\b
        if(mahal_dist <= mergeThresholdU)
            L1 = [L1, thisI]; % Indices of merged Gaussians
        end
    end
    
    % The new weight of the resulted merged Guassian is the summation of the weights of the Gaussian components. 
    w_bar_k_l = sum([updateIntensity.w1{L1}]);
    prunedIntensity.w1{l1} = w_bar_k_l;
    
    % The new mean of the merged Gaussian is the weighted average of the
    % merged means of Gaussian components.
    m_val = 0;
    for i = 1:length(L1)
        thisI = L1(i);
        m_val = m_val + (updateIntensity.w1{thisI} * updateIntensity.m1{thisI});
    end
    m_bar_k_l = m_val/w_bar_k_l;
    prunedIntensity.m1{l1} = m_bar_k_l;
    
    % Calculating covariance P_bar_k is a bit trickier
    P_val = zeros(4,4); 
  
    for i = 1:length(L1)
        thisI = L1(i);
        delta_m = m_bar_k_l - updateIntensity.m1{thisI};  
        P_val = P_val + updateIntensity.w1{thisI} * (updateIntensity.p1{thisI} + delta_m * delta_m');
    end
    P_bar_k_l = P_val / w_bar_k_l;
    prunedIntensity.p1{l1} = P_bar_k_l;

    % Now delete the elements in L from I
    for i = 1:length(L1)
        iToRemove = find(I1 == L1(i));
        I1(iToRemove) = [];
    end
end

numTargets_J_pruned1 = length(prunedIntensity.w1); % The number of targets after pruning

% Here you could do some check to see if numTargets_J_pruned > maxGaussiansJ
% and if needed delete some of the weaker gaussians. I haven't bothered but
% it might be useful for your implementation.

numTargets_Jk_minus_11 = numTargets_J_pruned1; % Number of targets in total, passed into the next filter iteration


%% For the second target type

%% Prune out the low-weighted targets
I2 = find([updateIntensity.w2{:}] >= T); % Find targets with high enough weights, I is indices of truncated weights

%% Merge the close-together targets
l2 = 0; % counts number of features

while isempty(I2) == false  % We delete from I as we merge
    
    l2 = l2 + 1;
    %Find j, which is i corresponding to highest updateIntensity.w for all i in I
    highWeights = updateIntensity.w2{I2};
    [maxW, j] = max(highWeights);
    j = j(1); % In case of two targets with equal weight
    % j is an index of highWeights (i.e. a position in I)
    % We want the index in updateIntensity.w
    j = I2(j);
    
    % Find all points with Mahalanobis distance less than U from point updateIntensity.m(j)
    L2 = [];  % A vector of indices of merged Gaussians.
    for iterateI = 1:length(I2)
        thisI = I2(iterateI);
        delta_m = updateIntensity.m2{thisI} - updateIntensity.m2{j};
        %mahal_dist = delta_m' * (updateIntensity.p{thisI} \ delta_m); % Invert covariance via left division
        mahal_dist = delta_m' * inv(updateIntensity.p2{thisI}) * delta_m; % inv(A)*b == A\b
        if(mahal_dist <= mergeThresholdU)
            L2 = [L2, thisI]; % Indices of merged Gaussians
        end
    end
    
    % The new weight of the resulted merged Guassian is the summation of the weights of the Gaussian components. 
    w_bar_k_l = sum([updateIntensity.w2{L2}]);
    prunedIntensity.w2{l2} = w_bar_k_l;
    
    % The new mean of the merged Gaussian is the weighted average of the
    % merged means of Gaussian components.
    m_val = 0;
    for i = 1:length(L2)
        thisI = L2(i);
        m_val = m_val + (updateIntensity.w2{thisI} * updateIntensity.m2{thisI});
    end
    m_bar_k_l = m_val/w_bar_k_l;
    prunedIntensity.m2{l2} = m_bar_k_l;
    
    % Calculating covariance P_bar_k is a bit trickier
    P_val = zeros(4,4); 
  
    for i = 1:length(L2)
        thisI = L2(i);
        delta_m = m_bar_k_l - updateIntensity.m2{thisI};  
        P_val = P_val + updateIntensity.w2{thisI} * (updateIntensity.p2{thisI} + delta_m * delta_m');
    end
    P_bar_k_l = P_val / w_bar_k_l;
    prunedIntensity.p2{l2} = P_bar_k_l;

    % Now delete the elements in L from I
    for i = 1:length(L2)
        iToRemove = find(I2 == L2(i));
        I2(iToRemove) = [];
    end
end

numTargets_J_pruned2 = length(prunedIntensity.w2); % The number of targets after pruning

% Here you could do some check to see if numTargets_J_pruned > maxGaussiansJ
% and if needed delete some of the weaker gaussians. I haven't bothered but
% it might be useful for your implementation.

numTargets_Jk_minus_12 = numTargets_J_pruned2; % Number of targets in total, passed into the next filter iteration

%% For the third target type

%% Prune out the low-weighted targets
I3 = find([updateIntensity.w3{:}] >= T); % Find targets with high enough weights, I is indices of truncated weights

%% Merge the close-together targets
l3 = 0; % counts number of features

while isempty(I3) == false  % We delete from I as we merge
    
    l3 = l3 + 1;
    %Find j, which is i corresponding to highest updateIntensity.w for all i in I
    highWeights = updateIntensity.w3{I3};
    [maxW, j] = max(highWeights);
    j = j(1); % In case of two targets with equal weight
    % j is an index of highWeights (i.e. a position in I)
    % We want the index in updateIntensity.w
    j = I3(j);
    
    % Find all points with Mahalanobis distance less than U from point updateIntensity.m(j)
    L3 = [];  % A vector of indices of merged Gaussians.
    for iterateI = 1:length(I3)
        thisI = I3(iterateI);
        delta_m = updateIntensity.m3{thisI} - updateIntensity.m3{j};
        %mahal_dist = delta_m' * (updateIntensity.p{thisI} \ delta_m); % Invert covariance via left division
        mahal_dist = delta_m' * inv(updateIntensity.p3{thisI}) * delta_m; % inv(A)*b == A\b
        if(mahal_dist <= mergeThresholdU)
            L3 = [L3, thisI]; % Indices of merged Gaussians
        end
    end
    
    % The new weight of the resulted merged Guassian is the summation of the weights of the Gaussian components. 
    w_bar_k_l = sum([updateIntensity.w3{L3}]);
    prunedIntensity.w3{l3} = w_bar_k_l;
    
    % The new mean of the merged Gaussian is the weighted average of the
    % merged means of Gaussian components.
    m_val = 0;
    for i = 1:length(L3)
        thisI = L3(i);
        m_val = m_val + (updateIntensity.w3{thisI} * updateIntensity.m3{thisI});
    end
    m_bar_k_l = m_val/w_bar_k_l;
    prunedIntensity.m3{l3} = m_bar_k_l;
    
    % Calculating covariance P_bar_k is a bit trickier
    P_val = zeros(4,4); 
  
    for i = 1:length(L3)
        thisI = L3(i);
        delta_m = m_bar_k_l - updateIntensity.m3{thisI};  
        P_val = P_val + updateIntensity.w3{thisI} * (updateIntensity.p3{thisI} + delta_m * delta_m');
    end
    P_bar_k_l = P_val / w_bar_k_l;
    prunedIntensity.p3{l3} = P_bar_k_l;

    % Now delete the elements in L from I
    for i = 1:length(L3)
        iToRemove = find(I3 == L3(i));
        I3(iToRemove) = [];
    end
end

numTargets_J_pruned3 = length(prunedIntensity.w3); % The number of targets after pruning

% Here you could do some check to see if numTargets_J_pruned > maxGaussiansJ
% and if needed delete some of the weaker gaussians. I haven't bothered but
% it might be useful for your implementation.

numTargets_Jk_minus_13 = numTargets_J_pruned3; % Number of targets in total, passed into the next filter iteration

%% For the fourth target type

%% Prune out the low-weighted targets
I4 = find([updateIntensity.w4{:}] >= T); % Find targets with high enough weights, I is indices of truncated weights

%% Merge the close-together targets
l4 = 0; % counts number of features

while isempty(I4) == false  % We delete from I as we merge
    
    l4 = l4 + 1;
    %Find j, which is i corresponding to highest updateIntensity.w for all i in I
    highWeights = updateIntensity.w4{I4};
    [maxW, j] = max(highWeights);
    j = j(1); % In case of two targets with equal weight
    % j is an index of highWeights (i.e. a position in I)
    % We want the index in updateIntensity.w
    j = I4(j);
    
    % Find all points with Mahalanobis distance less than U from point updateIntensity.m(j)
    L4 = [];  % A vector of indices of merged Gaussians.
    for iterateI = 1:length(I4)
        thisI = I4(iterateI);
        delta_m = updateIntensity.m4{thisI} - updateIntensity.m4{j};
        %mahal_dist = delta_m' * (updateIntensity.p{thisI} \ delta_m); % Invert covariance via left division
        mahal_dist = delta_m' * inv(updateIntensity.p4{thisI}) * delta_m; % inv(A)*b == A\b
        if(mahal_dist <= mergeThresholdU)
            L4 = [L4, thisI]; % Indices of merged Gaussians
        end
    end
    
    % The new weight of the resulted merged Guassian is the summation of the weights of the Gaussian components. 
    w_bar_k_l = sum([updateIntensity.w4{L4}]);
    prunedIntensity.w4{l4} = w_bar_k_l;
    
    % The new mean of the merged Gaussian is the weighted average of the
    % merged means of Gaussian components.
    m_val = 0;
    for i = 1:length(L4)
        thisI = L4(i);
        m_val = m_val + (updateIntensity.w4{thisI} * updateIntensity.m4{thisI});
    end
    m_bar_k_l = m_val/w_bar_k_l;
    prunedIntensity.m4{l4} = m_bar_k_l;
    
    % Calculating covariance P_bar_k is a bit trickier
    P_val = zeros(4,4); 
  
    for i = 1:length(L4)
        thisI = L4(i);
        delta_m = m_bar_k_l - updateIntensity.m4{thisI};  
        P_val = P_val + updateIntensity.w4{thisI} * (updateIntensity.p4{thisI} + delta_m * delta_m');
    end
    P_bar_k_l = P_val / w_bar_k_l;
    prunedIntensity.p4{l4} = P_bar_k_l;

    % Now delete the elements in L from I
    for i = 1:length(L4)
        iToRemove = find(I4 == L4(i));
        I4(iToRemove) = [];
    end
end

numTargets_J_pruned4 = length(prunedIntensity.w4); % The number of targets after pruning

% Here you could do some check to see if numTargets_J_pruned > maxGaussiansJ
% and if needed delete some of the weaker gaussians. I haven't bothered but
% it might be useful for your implementation.

numTargets_Jk_minus_14 = numTargets_J_pruned4; % Number of targets in total, passed into the next filter iteration

N_merged = sprintf('Number of targets after merging by summing weights = %0.0f', round(sum([prunedIntensity.w1{:}])) + round(sum([prunedIntensity.w2{:}])) + round(sum([prunedIntensity.w3{:}])) + round(sum([prunedIntensity.w4{:}])));
disp(N_merged);
 
%%
