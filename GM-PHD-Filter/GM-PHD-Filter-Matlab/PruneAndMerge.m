%% ------------------------------------------------------------------------
% pruneAndMerge.m
% -------------------------------------------------------------------------
%
% This file performs pruning and merging after the PHD update.
% The PHD update creates a combinatorial explosion in the number of Gaussian 
% components. We prune the low-weight ones and merge the close-together ones.
% The weight threshold T and distance threshold U that control this process 
% are set in Set_model.m function.
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: 
% Modified: March, 2014
% -------------------------------------------------------------------------

%%

function prunedMergedIntensity = PruneAndMerge(updatedIntensity, model)

prunedMergedIntensity.w = {};
prunedMergedIntensity.m = {};
prunedMergedIntensity.p = {};

%% Prune out the low-weighted targets
I = find([updatedIntensity.w{:}] >= model.T); % Find targets with high enough weights, I is indices of truncated weights

%% Merge the close-together targets
l = 0; % counts number of features

while isempty(I) == false  % We delete from I as we merge
    
    l = l + 1;
    %Find j, which is i corresponding to highest updateIntensity.w for all i in I
    highWeights = updatedIntensity.w{I};
    [maxW, j] = max(highWeights);
    j = j(1); % In case of two targets with equal weight
    % j is an index of highWeights (i.e. a position in I)
    % We want the index in updatedIntensity.w
    j = I(j);
    
    % Find all points with Mahalanobis distance less than U from point updatedIntensity.m(j)
    L = [];  % A vector of indices of merged Gaussians.
    for iterateI = 1:length(I)
        thisI = I(iterateI);
        delta_m = updatedIntensity.m{thisI} - updatedIntensity.m{j};
        %mahal_dist = delta_m' * (updatedIntensity.p{thisI} \ delta_m); % Invert covariance via left division
        mahal_dist = delta_m' * inv(updatedIntensity.p{thisI}) * delta_m; % inv(A)*b == A\b
        if(mahal_dist <= model.U)
            L = [L, thisI]; % Indices of merged Gaussians
        end
    end
    
    % The new weight of the resulted merged Guassian is the summation of the weights of the Gaussian components. 
    w_bar_k_l = sum([updatedIntensity.w{L}]);
    prunedMergedIntensity.w{l} = w_bar_k_l;
    
    % The new mean of the merged Gaussian is the weighted average of the
    % merged means of Gaussian components.
    m_val = 0;
    for i = 1:length(L)
        thisI = L(i);
        m_val = m_val + (updatedIntensity.w{thisI} * updatedIntensity.m{thisI});
    end
    m_bar_k_l = m_val/w_bar_k_l;
    prunedMergedIntensity.m{l} = m_bar_k_l;
    
    % Calculating covariance P_bar_k is a bit trickier
    P_val = zeros(4,4);
    for i = 1:length(L)
        thisI = L(i);
        delta_m = m_bar_k_l - updatedIntensity.m{thisI};  
        P_val = P_val + updatedIntensity.w{thisI} * (updatedIntensity.p{thisI} + delta_m * delta_m');
    end
    P_bar_k_l = P_val / w_bar_k_l;
    prunedMergedIntensity.p{l} = P_bar_k_l;

    % Now delete the elements in L from I
    for i = 1:length(L)
        iToRemove = find(I == L(i));
        I(iToRemove) = [];
    end
end

 
%%
