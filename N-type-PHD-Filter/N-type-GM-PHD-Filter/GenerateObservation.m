function Z = GenerateObservation(H, X, prob_detection, sigma_r, detect)

% This function generates measurements or observations using observation
% matrix, state vector X, std. deviation of the measurement, probability of
% detection.

    if (isempty(X) || detect > prob_detection)
        Z = [];
    else

        Z = H*X + sigma_r * randn(size(H,1), 1);
    end
end