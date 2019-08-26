%% ------------------------------------------------------------------------
% ExtactStates.m (for GM-PHD filter)
% -------------------------------------------------------------------------
%
% This pulls out every target with a weight over w_thresh (defined in 
% Set_model.m). There is the option of repeatedly printing out targets with
% rounded weights greater than 1. This will NOT change filter performance as
% the extracted state estimate is not fed back into the filter.
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: 
% Modified: March, 2014
% -------------------------------------------------------------------------

%%

function estimates = ExtractStates(prunedIntensity, model)
estimates.w = {};
estimates.m = {};
estimates.p = {};

% Extracting states from state extraction step
l = 0;
for i = 1:length(prunedIntensity.w)
    if (prunedIntensity.w{i} > model.w_thresh) 
       for j = 1:round(prunedIntensity.w{i}) % If a target has a rounded weight greater than 1, output it multiple times.          
           l = l+1;
           estimates.w{l} = prunedIntensity.w{i};
           estimates.m{l} = prunedIntensity.m{i};
           estimates.p{l} = prunedIntensity.p{i};
      end
    end
end

end
%%





