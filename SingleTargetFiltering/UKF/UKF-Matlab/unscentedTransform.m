function [X,W]= unscentedTransform(m,P,alpha,kappa)

% standard unscented transformation - note [choose alpha=1 => =lambda=kappa; giving offset of beta for first cov weight]

n_x = length(m); % Dimension of the new state vector which is the concatenation of the state and observation vectors. 
lambda = alpha^2*(n_x + kappa) - n_x;
Psqrtm = chol((n_x + lambda)*P)';% Calculate matrix square root of weighted covariance matrix. This is better than Psqrtm = real(sqrt((n_x + lambda)*P))  
X = repmat(m,[1 2*n_x+1])+ [ zeros(n_x,1) -Psqrtm Psqrtm ];    % Array of the sigma points
W = [ lambda 0.5*ones(1,2*n_x) ]/(n_x+lambda);    % Array of the weights for each sigma point

