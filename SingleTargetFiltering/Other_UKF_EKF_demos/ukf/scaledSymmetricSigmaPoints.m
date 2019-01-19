function [xPts, wPts, nPts] = scaledSymmetricSigmaPoints(x,P,alpha,beta,kappa)

% This function returns the scaled symmetric sigma point distribution.
%
%  [xPts, wPts, nPts] = scaledSymmetricSigmaPoints(x,P,alpha,beta,kappa)  
%
% Inputs:
%	 x	      mean
%	 P	      covariance
%        alpha        scaling parameter 1
%        beta         extra weight on zero'th point
%	 kappa	      scaling parameter 2 (usually set to default 0)
%
% Outputs:
%        xPts	 The sigma points
%        wPts	 The weights on the points
%	 nPts	 The number of points
%
%
%
% (C) 2000      Rudolph van der Merwe  
% (C) 1998-2000 S. J. Julier.

% Number of sigma points and scaling terms
n    = size(x(:),1);
nPts = 2*n+1;            % we're using the symmetric SUT

% Recalculate kappa according to scaling parameters
kappa = alpha^2*(n+kappa)-n;

% Allocate space

wPts=zeros(1,nPts);
xPts=zeros(n,nPts);

% Calculate matrix square root of weighted covariance matrix
Psqrtm=(chol((n+kappa)*P))';  

% Array of the sigma points
xPts=[zeros(size(P,1),1) -Psqrtm Psqrtm];

% Add mean back in
xPts = xPts + repmat(x,1,nPts);  

% Array of the weights for each sigma point
wPts=[kappa 0.5*ones(1,nPts-1) 0]/(n+kappa);

% Now calculate the zero'th covariance term weight
wPts(nPts+1) = wPts(1) + (1-alpha^2) + beta;

