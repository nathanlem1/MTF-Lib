% -----------------------------------------------------------------------
% Predict.m (for GM-PHD filter)
% -------------------------------------------------------------------------
% This implements the GM-PHD filter prediction step the GM-PHD filter 
% discussed in [1]. 
% It first implements measurement driven birth intensity 
% scheme which is very important in practice as it removes the need for the
% prior specification or knowledge of birth intensities. 
% This creates new targets (mean, weight, covariance) needed to instantiate 
% them in the next iteration. No spawning is considered as it has no special 
% advantage here, while measurement driven birth intensity is introduced at
% each time step using the position of the current measurements, with zero 
% initial velocity. The assumption of measurement-driven birth of targets
% is very crucial in visual tracking and this code has been used in work [2].
% The prediction of exisiting targets are also implemented in this
% function.
% 
% 
% References:
% [1] Vo, B. & Ma, W. The Gaussian mixture probability hypothesis density
%     filter Signal Processing, IEEE Transactions on, IEEE, 2006, 54
% [2]. Nathanael L. Baisa, "Online multi-object visual tracking using a 
% GM-PHD filter with deep appearance learning", 22nd International Conference 
% on Information Fusion (FUSION), July, 2019.
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: 
% Modified: March, 2014
% -------------------------------------------------------------------------

%% %% Step 1 and 2: (according to the original paper)

function predictedIntensity = Predict(Z_k, prunedIntensity, model)

predictedIntensity.w = {};
predictedIntensity.m = {};
predictedIntensity.p = {};

% Birth new targets
i = 0; % Used for index increment for the prediction step (both birth and existing targets)

for j = 1:length(Z_k)
    i = i + 1;
    m_i = Z_k(:,j);
    m_i(3:4) = [0; 0];
    w_i = model.w_birthsum/(size(Z_k,2));

    P_i = model.P; % Initialise the covariance

    predictedIntensity.w{i} = w_i;
    predictedIntensity.m{i} = m_i;
    predictedIntensity.p{i} = P_i;

end
numBirthedTargets = i;

% Predict existing targets
numTargets_Jk_minus_1 = length(prunedIntensity.w); % Number of Gaussian components after the pruning and merging step

for j = 1:numTargets_Jk_minus_1
    i = i + 1;
    predictedIntensity.w{i} = model.p_S * prunedIntensity.w{j};
    predictedIntensity.m{i} = model.F * prunedIntensity.m{j}; % Assume constant velocity.
    predictedIntensity.p{i} = model.Q + model.F * prunedIntensity.p{j} * model.F';
end

numTargets_Jk_k_minus_1 = i;  % numTargets_Jk_k_minus_1 = numTargets_Jk_minus_1 + numBirthedTargets;

end

%%

