% Process model for UKF
%
%  xout = ffun2(x,u,v,k)
%
%  INPUT
%         x    :  state vetor at time k
%         u    :  exogenous control input at time k
%         v    :  process noise vector at time k
%         k    :  time index
%  OUTPUT
%         xout :  state vector at time k+1
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

function xout = ffun2(x,u,v,k)

  xout = x + sin(x).*x + v;  