function newStates = gen_newstates(targetStates, model)

% This generates new states in each iteration i.e. this is ground truth
% states of targets that can be observed using a sensor. 
newStates = [];
for i = 1:size(targetStates,2)
    W = model.sigma_v*sqrt(model.Q)*randn(size(model.Q,2),size(targetStates(:,i),2)); % Process noise
    X = model.F*targetStates(:,i)+ W;
    newStates = [newStates, X];
end