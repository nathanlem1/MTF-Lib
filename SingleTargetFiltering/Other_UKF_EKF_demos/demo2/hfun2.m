% Observation model for UKF
%
%  y = hfun1(x,u,n,k)
%
%  INPUT
%         x    :  state vetor at time k
%         u    :  exogenous control input at time k
%         n    :  measurement noise vector at time k
%         k    :  time index
%  OUTPUT
%         y    :  state observation vector at time k
%
%  
% This function assumes a discrete time nonlinear dynamic state space system
% of the following format:
%
%   x(k+1) = ffun[x(k),u(k),v(k),k]
%     y(k) = hfun[x(k),u(k),n(k),k]
%
% Author : Rudolph van der Merwe
% Date   : 24 May 2001
%

function y = hfun1(x,u,n,k)
  
  y = x.^2 + n;