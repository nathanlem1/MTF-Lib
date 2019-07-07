%% ------------------------------------------------------------------------
% N_type_GM_PHD_Update.m
% -------------------------------------------------------------------------
%
% This file performs a N-type GM-PHD filter update on the targets. This is 
% basically a brute-force Kalman update of every target with every measurement 
% and creating a new target from the update results; for each target type 
% handling confusions from the other target types.
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

%%

s = sprintf('Step 6: Performing update.');
disp(s);

% Make sure that the previous targets states not be repeated in the update 
% and state extraction again!
updateIntensity.w1 = {};
updateIntensity.m1 = {};
updateIntensity.p1 = {};
updateIntensity.w2 = {};
updateIntensity.m2 = {};
updateIntensity.p2 = {};
updateIntensity.w3 = {};
updateIntensity.m3 = {};
updateIntensity.p3 = {};
updateIntensity.w4 = {};
updateIntensity.m4 = {};
updateIntensity.p4 = {};

% -----------Miss-detection part of GM-PHD update -------------------------

% First we assume that we failed to detect all targets i.e. miss-detection
% term (1 - P_D) Dk|k-1) of the update step!
% We scale all weights by probability of missed detection
% We already did the prediction step for these so their position &
% covariance will have been updated. What remains is to rescale their weight.

% For target type 1
for j = 1:numTargets_Jk_k_minus_11
    
    updateIntensity.w1{j} = (1 - prob_detection11)*predictedIntensity.w1{j}; % prob_detection = probability(targetType1/targetType1)
    updateIntensity.m1{j} = predictedIntensity.m1{j};
    updateIntensity.p1{j} = predictedIntensity.p1{j};   
end

% For target type 2
for j = 1:numTargets_Jk_k_minus_12
    
    updateIntensity.w2{j} = (1 - prob_detection22)*predictedIntensity.w2{j}; % prob_detection = probability(targetType2/targetType2)
    updateIntensity.m2{j} = predictedIntensity.m2{j};
    updateIntensity.p2{j} = predictedIntensity.p2{j};   
end

% For target type 3
for j = 1:numTargets_Jk_k_minus_13
    
    updateIntensity.w3{j} = (1 - prob_detection33)*predictedIntensity.w3{j}; % prob_detection = probability(targetType2/targetType2)
    updateIntensity.m3{j} = predictedIntensity.m3{j};
    updateIntensity.p3{j} = predictedIntensity.p3{j};   
end

% For target type 4
for j = 1:numTargets_Jk_k_minus_14
    
    updateIntensity.w4{j} = (1 - prob_detection44)*predictedIntensity.w4{j}; % prob_detection = probability(targetType2/targetType2)
    updateIntensity.m4{j} = predictedIntensity.m4{j};
    updateIntensity.p4{j} = predictedIntensity.p4{j};   
end

% ------------Detection part of GM-PHD update -----------------------------
% Now we update all combinations of matching every observation with every
% target in the map i.e. the |Zk| detection terms of the update step!

% For target type 1
l1 = 0;
for zi = 1:size(Z1,2)
    l1 = l1 + 1; % l is used to calculate an offset from previous updates. It maxes out at l = number_of_measurements. A more elegant solution would be to set l = zi but I retain this method for compatibility with Vo&Ma
    
    for j = 1:numTargets_Jk_k_minus_11
        thisZ = Z1(:,zi); % This consists of the observed position.

        % Update weights. Updating weights over all components (position and velocity)
        % produces unacceptably low weights; therefore, only position
        % components will be used
        updateIntensity.w1{l1*numTargets_Jk_k_minus_11+j} = prob_detection11 * predictedIntensity.w1{j}*mvnpdf(thisZ(1:2),constructUpdateIntensity.eta1{j}(1:2),constructUpdateIntensity.S1{j}(1:2,1:2)); % Hoping normpdf is the function to be used!       
        updateIntensity.m1{l1*numTargets_Jk_k_minus_11+j} = predictedIntensity.m1{j} + constructUpdateIntensity.K1{j}*(thisZ - constructUpdateIntensity.eta1{j});
        updateIntensity.p1{l1*numTargets_Jk_k_minus_11+j} = constructUpdateIntensity.p1{j};
    end
    
    % Sum up weights for use in reweighting
    totalWeight1 = 0;
    for i = 1:numTargets_Jk_k_minus_11
        totalWeight1 = totalWeight1 + updateIntensity.w1{l1*numTargets_Jk_k_minus_11+i};
    end
    
    K2 = 0;
    K3 = 0;
    K4 = 0;
    for j = 1:numTargets_Jk_k_minus_11
        measZ = [Z1(1,zi), Z1(2,zi)];
        K_k = clutter_intensity1(lambda_1,measZ);
        if N_type_GM_PHD == 1   % Decide which filter to use: N-type GM-PHD filter (1) or N independent GM-PHD filters (0)
            for i = 1:numTargets_Jk_k_minus_12  % Uses predicted intensity (PHD) (of target type 2)
                K2 =  K2 + predictedIntensity.w2{i} * mvnpdf(measZ', model.H21(1:2,:)*predictedIntensity.m2{i},model.H21(1:2,:)*predictedIntensity.p2{i}*model.H21(1:2,:)'); % ONLY position component is used. Hoping normpdf is the function to be used!       
            end
            K2 = prob_detection21*K2; 
            for i = 1:numTargets_Jk_k_minus_13  % Uses predicted intensity (PHD) (of target type 3)
                K3 =  K3 + predictedIntensity.w3{i} * mvnpdf(measZ', model.H31(1:2,:)*predictedIntensity.m3{i},model.H31(1:2,:)*predictedIntensity.p3{i}*model.H31(1:2,:)'); % ONLY position component is used. Hoping normpdf is the function to be used!       
            end
            K3 = prob_detection31*K3; % Even using prob_detection22, it can differentiate the target types.
            for i = 1:numTargets_Jk_k_minus_14  % Uses predicted intensity (PHD) (of target type 3)
                K4 =  K4 + predictedIntensity.w4{i} * mvnpdf(measZ', model.H41(1:2,:)*predictedIntensity.m4{i},model.H41(1:2,:)*predictedIntensity.p4{i}*model.H41(1:2,:)'); % ONLY position component is used. Hoping normpdf is the function to be used!       
            end
            K4 = prob_detection41*K4; % Even using prob_detection22, it can differentiate the target types.
        end
        K_k = K2 + K3 + K4 + K_k;  % Now the clutter due to the PHD estimator and likelihood of the second target type is included
        updateIntensity.w1{l1*numTargets_Jk_k_minus_11+j} = updateIntensity.w1{l1*numTargets_Jk_k_minus_11+j}/(K_k + totalWeight1);
    end   
end

numTargets_Jk11 = l1*numTargets_Jk_k_minus_11 + numTargets_Jk_k_minus_11;

% For target type 2
l2 = 0;
for zi = 1:size(Z2,2)
    l2 = l2 + 1; % l is used to calculate an offset from previous updates. It maxes out at l = number_of_measurements. A more elegant solution would be to set l = zi but I retain this method for compatibility with Vo&Ma
    
    for j = 1:numTargets_Jk_k_minus_12
        thisZ = Z2(:,zi); % This consists of the observed position.

        % Update weights. Updating weights over all components (position and velocity)
        % produces unacceptably low weights; therefore, only position
        % components will be used
        updateIntensity.w2{l2*numTargets_Jk_k_minus_12+j} = prob_detection22 * predictedIntensity.w2{j}*mvnpdf(thisZ(1:2),constructUpdateIntensity.eta2{j}(1:2),constructUpdateIntensity.S2{j}(1:2,1:2)); % Hoping normpdf is the function to be used!       
        updateIntensity.m2{l2*numTargets_Jk_k_minus_12+j} = predictedIntensity.m2{j} + constructUpdateIntensity.K2{j}*(thisZ - constructUpdateIntensity.eta2{j});
        updateIntensity.p2{l2*numTargets_Jk_k_minus_12+j} = constructUpdateIntensity.p2{j};
    end
    
    % Sum up weights for use in reweighting
    totalWeight2 = 0;
    for i = 1:numTargets_Jk_k_minus_12
        totalWeight2 = totalWeight2 + updateIntensity.w2{l2*numTargets_Jk_k_minus_12+i};
    end
    
    K1 = 0;
    K3 = 0;
    K4 = 0;
    for j = 1:numTargets_Jk_k_minus_12
        measZ = [Z2(1,zi), Z2(2,zi)];
        K_k = clutter_intensity2(lambda_2,measZ);
        if N_type_GM_PHD == 1 % Decide which filter to use: N-type GM-PHD filter (1) or N independent GM-PHD filters (0)
            for i = 1:numTargets_Jk_k_minus_11 % Uses predicted intensity (PHD)
                K1 =  K1 + predictedIntensity.w1{i} * mvnpdf(measZ', model.H12(1:2,:)*predictedIntensity.m1{i},model.H12(1:2,:)*predictedIntensity.p1{i}*model.H12(1:2,:)'); % ONLY position component is used. Hoping normpdf is the function to be used!       
            end
            K1 = prob_detection12*K1; 
            for i = 1:numTargets_Jk_k_minus_13 % Uses predicted intensity (PHD)
                K3 =  K3 + predictedIntensity.w3{i} * mvnpdf(measZ', model.H32(1:2,:)*predictedIntensity.m3{i},model.H32(1:2,:)*predictedIntensity.p3{i}*model.H32(1:2,:)'); % ONLY position component is used. Hoping normpdf is the function to be used!       
            end
            K3 = prob_detection32*K3; 
            for i = 1:numTargets_Jk_k_minus_14 % Uses predicted intensity (PHD)
                K4 =  K4 + predictedIntensity.w4{i} * mvnpdf(measZ', model.H42(1:2,:)*predictedIntensity.m4{i},model.H42(1:2,:)*predictedIntensity.p4{i}*model.H42(1:2,:)'); % ONLY position component is used. Hoping normpdf is the function to be used!       
            end
            K4 = prob_detection42*K4; 
        end
        K_k = K1 + K3 + K4 + K_k; % Now the clutter due to the PHD estimator and likelihood of the first target type is included
        updateIntensity.w2{l2*numTargets_Jk_k_minus_12+j} = updateIntensity.w2{l2*numTargets_Jk_k_minus_12+j}/(K_k + totalWeight2);
    end   
end

numTargets_Jk22 = l2*numTargets_Jk_k_minus_12 + numTargets_Jk_k_minus_12;

% For target type 3
l3 = 0;
for zi = 1:size(Z3,2)
    l3 = l3 + 1; % l is used to calculate an offset from previous updates. It maxes out at l = number_of_measurements. A more elegant solution would be to set l = zi but I retain this method for compatibility with Vo&Ma
    
    for j = 1:numTargets_Jk_k_minus_13
        thisZ = Z3(:,zi); % This consists of the observed position.

        % Update weights. Updating weights over all components (position and velocity)
        % produces unacceptably low weights; therefore, only position
        % components will be used
        updateIntensity.w3{l3*numTargets_Jk_k_minus_13+j} = prob_detection33 * predictedIntensity.w3{j}*mvnpdf(thisZ(1:2),constructUpdateIntensity.eta3{j}(1:2),constructUpdateIntensity.S3{j}(1:2,1:2)); % Hoping normpdf is the function to be used!       
        updateIntensity.m3{l3*numTargets_Jk_k_minus_13+j} = predictedIntensity.m3{j} + constructUpdateIntensity.K3{j}*(thisZ - constructUpdateIntensity.eta3{j});
        updateIntensity.p3{l3*numTargets_Jk_k_minus_13+j} = constructUpdateIntensity.p3{j};
    end
    
    % Sum up weights for use in reweighting
    totalWeight3 = 0;
    for i = 1:numTargets_Jk_k_minus_13
        totalWeight3 = totalWeight3 + updateIntensity.w3{l3*numTargets_Jk_k_minus_13+i};
    end
    
    K1 = 0;
    K2 = 0;
    K4 = 0;
    for j = 1:numTargets_Jk_k_minus_13
        measZ = [Z3(1,zi), Z3(2,zi)];
        K_k = clutter_intensity3(lambda_3,measZ);
        if N_type_GM_PHD == 1  % Decide which filter to use: N-type GM-PHD filter (1) or N independent GM-PHD filters (0)
            for i = 1:numTargets_Jk_k_minus_11 % Uses predicted intensity (PHD)
                K1 =  K1 + predictedIntensity.w1{i} * mvnpdf(measZ', model.H13(1:2,:)*predictedIntensity.m1{i},model.H13(1:2,:)*predictedIntensity.p1{i}*model.H13(1:2,:)'); % ONLY position component is used. Hoping normpdf is the function to be used!       
            end
            K1 = prob_detection13*K1; 
            for i = 1:numTargets_Jk_k_minus_12 % Uses predicted intensity (PHD)
                K2 =  K2 + predictedIntensity.w2{i} * mvnpdf(measZ', model.H23(1:2,:)*predictedIntensity.m2{i},model.H23(1:2,:)*predictedIntensity.p2{i}*model.H23(1:2,:)'); % ONLY position component is used. Hoping normpdf is the function to be used!       
            end
            K2 = prob_detection23*K2; 
            for i = 1:numTargets_Jk_k_minus_14 % Uses predicted intensity (PHD)
                K4 =  K4 + predictedIntensity.w4{i} * mvnpdf(measZ', model.H43(1:2,:)*predictedIntensity.m4{i},model.H43(1:2,:)*predictedIntensity.p4{i}*model.H43(1:2,:)'); % ONLY position component is used. Hoping normpdf is the function to be used!       
            end
            K4 = prob_detection43*K4; 
        end
        K_k = K1 + K2 + K4 + K_k; % Now the clutter due to the PHD estimator and likelihood of the first target type is included
        updateIntensity.w3{l3*numTargets_Jk_k_minus_13+j} = updateIntensity.w3{l3*numTargets_Jk_k_minus_13+j}/(K_k + totalWeight3);
    end   
end

numTargets_Jk33 = l3*numTargets_Jk_k_minus_13 + numTargets_Jk_k_minus_13;

% For target type 4
l4 = 0;
for zi = 1:size(Z4,2)
    l4 = l4 + 1; % l is used to calculate an offset from previous updates. It maxes out at l = number_of_measurements. A more elegant solution would be to set l = zi but I retain this method for compatibility with Vo&Ma
    
    for j = 1:numTargets_Jk_k_minus_14
        thisZ = Z4(:,zi); % This consists of the observed position.

        % Update weights. Updating weights over all components (position and velocity)
        % produces unacceptably low weights; therefore, only position
        % components will be used
        updateIntensity.w4{l4*numTargets_Jk_k_minus_14+j} = prob_detection44 * predictedIntensity.w4{j}*mvnpdf(thisZ(1:2),constructUpdateIntensity.eta4{j}(1:2),constructUpdateIntensity.S4{j}(1:2,1:2)); % Hoping normpdf is the function to be used!       
        updateIntensity.m4{l4*numTargets_Jk_k_minus_14+j} = predictedIntensity.m4{j} + constructUpdateIntensity.K4{j}*(thisZ - constructUpdateIntensity.eta4{j});
        updateIntensity.p4{l4*numTargets_Jk_k_minus_14+j} = constructUpdateIntensity.p4{j};
    end
    
    % Sum up weights for use in reweighting
    totalWeight4 = 0;
    for i = 1:numTargets_Jk_k_minus_14
        totalWeight4 = totalWeight4 + updateIntensity.w4{l4*numTargets_Jk_k_minus_14+i};
    end
    
    K1 = 0;
    K2 = 0;
    K3 = 0;
    for j = 1:numTargets_Jk_k_minus_14
        measZ = [Z4(1,zi), Z4(2,zi)];
        K_k = clutter_intensity4(lambda_4,measZ);
        if N_type_GM_PHD == 1  % Decide which filter to use: N-type GM-PHD filter (1) or N independent GM-PHD filters (0)
            for i = 1:numTargets_Jk_k_minus_11 % Uses predicted intensity (PHD)
                K1 =  K1 + predictedIntensity.w1{i} * mvnpdf(measZ', model.H14(1:2,:)*predictedIntensity.m1{i},model.H14(1:2,:)*predictedIntensity.p1{i}*model.H14(1:2,:)'); % ONLY position component is used. Hoping normpdf is the function to be used!       
            end
            K1 = prob_detection14*K1; 
            for i = 1:numTargets_Jk_k_minus_12 % Uses predicted intensity (PHD)
                K2 =  K2 + predictedIntensity.w2{i} * mvnpdf(measZ', model.H24(1:2,:)*predictedIntensity.m2{i},model.H24(1:2,:)*predictedIntensity.p2{i}*model.H24(1:2,:)'); % ONLY position component is used. Hoping normpdf is the function to be used!       
            end
            K2 = prob_detection24*K2; 
            for i = 1:numTargets_Jk_k_minus_13 % Uses predicted intensity (PHD)
                K3 =  K3 + predictedIntensity.w3{i} * mvnpdf(measZ', model.H34(1:2,:)*predictedIntensity.m3{i},model.H34(1:2,:)*predictedIntensity.p3{i}*model.H34(1:2,:)'); % ONLY position component is used. Hoping normpdf is the function to be used!       
            end
            K3 = prob_detection34*K3; 
        end
        K_k = K1 + K2 + K3 + K_k; % Now the clutter due to the PHD estimator and likelihood of the first target type is included
        updateIntensity.w4{l4*numTargets_Jk_k_minus_14+j} = updateIntensity.w4{l4*numTargets_Jk_k_minus_14+j}/(K_k + totalWeight4);
    end   
end

numTargets_Jk44 = l4*numTargets_Jk_k_minus_14 + numTargets_Jk_k_minus_14;

N_update = sprintf('Number of targets after updating by summing weights = %0.0f', round(sum([updateIntensity.w1{:}])) + round(sum([updateIntensity.w2{:}])) + round(sum([updateIntensity.w3{:}])) + round(sum([updateIntensity.w4{:}])));
disp(N_update);

%%

