function observations = gen_observations(truthStates, model)

% This generates observations from the goundtruth states of targets.
observations = [];
for i = 1:size(truthStates,2)
    detect_i = rand;
    if detect_i <= model.p_D
        V = sqrt(model.R)*randn(size(model.R,2),size(truthStates(:,i),2)); % Observation noise
        Z = model.H*truthStates(:,i) + V;
        observations = [observations, Z];
    end
end