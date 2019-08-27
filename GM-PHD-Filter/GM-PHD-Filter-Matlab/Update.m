%% ------------------------------------------------------------------------
% Update.m (for GM-PHD filter)
% -------------------------------------------------------------------------
%
% This file performs a PHD filter update on the targets. This does Kalman 
% update of every target with every measurement.
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: 
% Modified: March, 2014
% -------------------------------------------------------------------------

%% Step 3 and 4: (according to the original paper)

function updatedIntensity = Update(Z, predictedIntensity, model)

updatedIntensity.w = {};
updatedIntensity.m = {};
updatedIntensity.p = {};

numTargets_Jk_k_minus_1 = length(predictedIntensity.w); % Number of Gaussian components after the prediction step

for j = 1:numTargets_Jk_k_minus_1
    
    constructUpdateIntensity.eta{j} = model.H*predictedIntensity.m{j};
    constructUpdateIntensity.S{j} = model.R + model.H*predictedIntensity.p{j}*model.H';
    constructUpdateIntensity.K{j} = predictedIntensity.p{j}*model.H'*inv(constructUpdateIntensity.S{j});
    constructUpdateIntensity.p{j} = predictedIntensity.p{j} - constructUpdateIntensity.K{j}*model.H*predictedIntensity.p{j};
end

% -----------Miss-detection part of GM-PHD update -------------------------
% We scale all weights by probability of missed detection (1 - p_D)

for j = 1:numTargets_Jk_k_minus_1
    
    updatedIntensity.w{j} = (1 - model.p_D)*predictedIntensity.w{j};
    updatedIntensity.m{j} = predictedIntensity.m{j};
    updatedIntensity.p{j} = predictedIntensity.p{j};   
end

% ------------Detection part of GM-PHD update -----------------------------
% Every observation updates every target in the map
l = 0;
for zi = 1:size(Z,2)
    l = l + 1; % l is used to calculate an offset from previous updates
    
    for j = 1:numTargets_Jk_k_minus_1
        thisZ = Z(:,zi); % This consists of the observed position.
        %updatedIntensity.w{l*numTargets_Jk_k_minus_1+j} = model.p_D * predictedIntensity.w{j}*mvnpdf(thisZ(1:2),constructUpdateIntensity.eta{j}(1:2),constructUpdateIntensity.S{j}(1:2,1:2)); % Hoping multivariate_normal.pdf is the right one to use; this is for video [x, y, w, h]
        updatedIntensity.w{l*numTargets_Jk_k_minus_1+j} = model.p_D * predictedIntensity.w{j}*mvnpdf(thisZ,constructUpdateIntensity.eta{j},constructUpdateIntensity.S{j}); % Hoping multivariate_normal.pdf is the right one to use; this is for simulation [x, y]         
        updatedIntensity.m{l*numTargets_Jk_k_minus_1+j} = predictedIntensity.m{j} + constructUpdateIntensity.K{j}*(thisZ - constructUpdateIntensity.eta{j});
        updatedIntensity.p{l*numTargets_Jk_k_minus_1+j} = constructUpdateIntensity.p{j};
    end
    
    % Sum up weights for use in reweighting
    totalWeight = 0;
    for i = 1:numTargets_Jk_k_minus_1
        totalWeight = totalWeight + updatedIntensity.w{l*numTargets_Jk_k_minus_1+i};
    end
    
    for j = 1:numTargets_Jk_k_minus_1
        measZ = [Z(1,zi), Z(2,zi)];
        K_k = model.clutterIntensity(measZ);
        updatedIntensity.w{l*numTargets_Jk_k_minus_1+j} = updatedIntensity.w{l*numTargets_Jk_k_minus_1+j}/(K_k + totalWeight);
    end   
end

%%

