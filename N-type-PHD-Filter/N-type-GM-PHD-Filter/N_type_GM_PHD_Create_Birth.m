%% -----------------------------------------------------------------------
% N_type_GM_PHD_Create_Birth.m 
% -------------------------------------------------------------------------
% This implements measurement driven birth intensity which is very
% important in practice as it removes the need for the prior specification
% or knowledge of birth intensities. This populates the birth lists i.e. creates 
% new targets (mean, weight, covariance) needed to instantiate them in the next iteration. 
% No spawning is considered as it has no special advantage here, while
% measurement driven birth intensity is introduced at each time step using 
% the position of the current measurements, with zero initial velocity. 
% The velocity component of the state vector is estimated in subsequent iterations. 
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

%%

disp('Step 1: Creating new targets from measurements');

i1 = 0; % Used for index increment throughout the whole code, for target type 1.
i2 = 0; % Used for index increment throughout the whole code, for target type 2.
i3 = 0; % Used for index increment throughout the whole code, for target type 3.
i4 = 0; % Used for index increment throughout the whole code, for target type 4.
numBirthedTargets1 = 0;
numBirthedTargets2 = 0;
numBirthedTargets3 = 0;
numBirthedTargets4 = 0;

if(k >= 1)
       % For target type 1
        thisMeas1 = simMeasurementHistory.Z1{k};
    
        for j = 1:size(thisMeas1,2)
                    i1 = i1 + 1;
                    m_i = thisMeas1(:,j);
                    m_i(3:4) = [0; 0];
                    w_i = w_birthsum/(size(thisMeas1,2));
                    
                    P_i = covariance_birth; %Initialise the covariance
                    
                    predictedIntensity.w1{i1} = w_i;
                    predictedIntensity.m1{i1} = m_i;
                    birthIntensity.m1{i1} = m_i;  % for observing velocity component.
                    predictedIntensity.p1{i1} = P_i;

        end
        
        numBirthedTargets1 = i1;
        
        % For target type 2        
        thisMeas2 = simMeasurementHistory.Z2{k};
    
        for j = 1:size(thisMeas2,2)
                    i2 = i2 + 1;
                    m_i = thisMeas2(:,j);
                    m_i(3:4) = [0; 0];
                    w_i = w_birthsum/(size(thisMeas2,2));
                    
                    P_i = covariance_birth; %Initialise the covariance
                    
                    predictedIntensity.w2{i2} = w_i;
                    predictedIntensity.m2{i2} = m_i;
                    birthIntensity.m2{i2} = m_i;  % for observing velocity component.
                    predictedIntensity.p2{i2} = P_i;

        end
        
        numBirthedTargets2 = i2;
        
                % For target type 3        
        thisMeas3 = simMeasurementHistory.Z3{k};
    
        for j = 1:size(thisMeas3,2)
                    i3 = i3 + 1;
                    m_i = thisMeas3(:,j);
                    m_i(3:4) = [0; 0];
                    w_i = w_birthsum/(size(thisMeas3,2));
                    
                    P_i = covariance_birth; %Initialise the covariance
                    
                    predictedIntensity.w3{i3} = w_i;
                    predictedIntensity.m3{i3} = m_i;
                    birthIntensity.m3{i3} = m_i;  % for observing velocity component.
                    predictedIntensity.p3{i3} = P_i;

        end
        
        numBirthedTargets3 = i3;
        
                % For target type 4        
        thisMeas4 = simMeasurementHistory.Z4{k};
    
        for j = 1:size(thisMeas4,2)
                    i4 = i4 + 1;
                    m_i = thisMeas4(:,j);
                    m_i(3:4) = [0; 0];
                    w_i = w_birthsum/(size(thisMeas4,2));
                    
                    P_i = covariance_birth; %Initialise the covariance
                    
                    predictedIntensity.w4{i4} = w_i;
                    predictedIntensity.m4{i4} = m_i;
                    birthIntensity.m4{i4} = m_i;  % for observing velocity component.
                    predictedIntensity.p4{i4} = P_i;

        end
        
        numBirthedTargets4 = i4;
end

N_birthed = sprintf('Number of created (measured) targets and clutters altogether = %0.0f', numBirthedTargets1 + numBirthedTargets2 + numBirthedTargets3 + numBirthedTargets4);
disp(N_birthed);

%%


