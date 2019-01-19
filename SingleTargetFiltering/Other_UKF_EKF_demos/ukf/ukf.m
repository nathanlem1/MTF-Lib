function [xEst,PEst,xPred,PPred,zPred,inovation,S,K]=ukf(xEst,PEst,U,Q,ffun,z,R,hfun,dt,alpha,beta,kappa);

% TITLE    :  UNSCENTED KALMAN FILTER  
%
% PURPOSE  :  This function performs one complete step of the unscented Kalman filter.
%
% SYNTAX   :  [xEst,PEst,xPred,PPred,zPred,inovation]=ukf(xEst,PEst,U,Q,ffun,z,R,hfun,dt,alpha,beta,kappa)
%
% INPUTS   :  - xEst             : state mean estimate at time k  
%             - PEst             : state covariance at time k
%             - U                : vector of control inputs
%             - Q                : process noise covariance at time k  
%             - z                : observation at k+1  
%             - R                : measurement noise covariance at k+1  
%             - ffun             : process model function  
%             - hfun             : observation model function  
%             - dt               : time step (passed to ffun/hfun)   
%	      - alpha (optional) : sigma point scaling parameter. Defaults to 1.
%             - beta  (optional) : higher order error scaling parameter. Default to 0.  
%             - kappa (optional) : scalar tuning parameter 1. Defaults to 0.  
%
% OUTPUTS  :  - xEst             : updated estimate of state mean at time k+1
%	      - PEst             : updated state covariance at time k+1
%             - xPred            : prediction of state mean at time k+1
%             - PPred            : prediction of state covariance at time k+1
%	      - inovation        : innovation vector
%  
% AUTHORS  :  Simon J. Julier       (sjulier@erols.com)    1998-2000
%             Rudolph van der Merwe (rvdmerwe@ece.ogi.edu) 2000
%
% DATE     :  14 August 2000
%
% NOTES    :  The process model is of the form, x(k+1) = ffun[x(k),v(k),dt,u(k)]
%             where v(k) is the process noise vector. The observation model is 
%             of the form, z(k) = hfun[x(k),w(k),dt,u(k)], where w(k) is the 
%             observation noise vector.
%
%             This code was written to be readable. There is significant
%             scope for optimisation even in Matlab.
%
  

% Process defaults

if (nargin < 10)
  alpha=1;
end;

if (nargin < 11)
  beta=0;
end;

if (nargin < 12)
  kappa=3-size(xEst,1);
end;


% Calculate the dimensions of the problem and a few useful
% scalars

states       = size(xEst(:),1);
observations = size(z(:),1);
vNoise       = size(Q,2);
wNoise       = size(R,2);

noises       = vNoise+wNoise;

% Augment the state vector with the noise vectors.
% Note: For simple, additive noise models this part
% can be done differently to save on computational cost.
% For details, contact Rudolph v.d. Merwe

if (noises)
  N=[Q zeros(vNoise,wNoise); zeros(wNoise,vNoise) R];
  PQ=[PEst zeros(states,noises);zeros(noises,states) N];
  xQ=[xEst;zeros(noises,1)];
else
  PQ=PEst;
  xQ=xEst;
end;

% Calculate the sigma points and there corresponding weights using the Scaled Unscented
% Transformation
[xSigmaPts, wSigmaPts, nsp] = scaledSymmetricSigmaPoints(xQ, PQ, alpha, beta, kappa); 


% Duplicate wSigmaPts into matrix for code speedup
wSigmaPts_xmat = repmat(wSigmaPts(:,2:nsp),states,1);
wSigmaPts_zmat = repmat(wSigmaPts(:,2:nsp),observations,1);


% Work out the projected sigma points and their means
% This routine is fairly generic. The only thing to watch out for are
% angular discontinuities. There is a standard trick for doing this -
% contact me (Julier) for details!

xPredSigmaPts = feval(ffun,xSigmaPts(1:states,:),repmat(U(:),1,nsp),xSigmaPts(states+1:states+vNoise,:),dt);
zPredSigmaPts = feval(hfun,xPredSigmaPts,repmat(U(:),1,nsp),xSigmaPts(states+vNoise+1:states+noises,:),dt);

% Calculate the mean. Based on discussions with C. Schaefer, form
% is chosen to maximise numerical robustness.
% - I vectorized this part of the code for a speed increase : RvdM 2000

xPred = sum(wSigmaPts_xmat .* (xPredSigmaPts(:,2:nsp) - repmat(xPredSigmaPts(:,1),1,nsp-1)),2);
zPred = sum(wSigmaPts_zmat .* (zPredSigmaPts(:,2:nsp) - repmat(zPredSigmaPts(:,1),1,nsp-1)),2);

xPred=xPred+xPredSigmaPts(:,1);
zPred=zPred+zPredSigmaPts(:,1);

% Work out the covariances and the cross correlations. Note that
% the weight on the 0th point is different from the mean
% calculation due to the scaled unscented algorithm.

exSigmaPt = xPredSigmaPts(:,1)-xPred;
ezSigmaPt = zPredSigmaPts(:,1)-zPred;

PPred   = wSigmaPts(nsp+1)*exSigmaPt*exSigmaPt';
PxzPred = wSigmaPts(nsp+1)*exSigmaPt*ezSigmaPt';
S       = wSigmaPts(nsp+1)*ezSigmaPt*ezSigmaPt';

exSigmaPt = xPredSigmaPts(:,2:nsp) - repmat(xPred,1,nsp-1);
ezSigmaPt = zPredSigmaPts(:,2:nsp) - repmat(zPred,1,nsp-1);
PPred     = PPred + (wSigmaPts_xmat .* exSigmaPt) * exSigmaPt';
S         = S + (wSigmaPts_zmat .* ezSigmaPt) * ezSigmaPt';
PxzPred   = PxzPred + exSigmaPt * (wSigmaPts_zmat .* ezSigmaPt)';


%%%%% MEASUREMENT UPDATE

% Calculate Kalman gain
K  = PxzPred / S;

% Calculate Innovation
inovation = z - zPred;

% Update mean
xEst = xPred + K*inovation;

% Update covariance
PEst = PPred - K*S*K';



