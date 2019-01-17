% This function implements Kalman filter (KF) in a clear and understandable
% manner. The tutorial for the KF is given in the following papers:
% 1. Greg Welch and Gary Bishop, 'An Introduction to the Kalman Filter'
% 2. M.S.Arulampalam, S.Maskel, N. Gordon and T Clapp, 'A Tutorial on 
% Particle Filters for Online Nonlinear/Non-Gaussian Bayesian Tracking'.
%
% The state dynamics evolves according to the following equations:
% x(k) = f(k)(x(k-1), u(k), w(k-1)) = F(k)x(k-1) + G(k)u(k) + w(k-1)
%                  f(k)(.) is a linear (in KF case) or non-linear) function.
%                  Meaning the state vector x evolves during one time
%                  step pre-multiplying by the state transition
%                  matrix F. There is optionally (if nonzero) an input
%                  (control) vector u which affects the state linearly, 
%                  and this linear effect on the state is represented 
%                  pre-multiplying by the input matrix G. There is also
%                  gaussian process noise w. In our case, we set u to 0.
% where w ~ N(0,Q(k-1)) meaning w(k-1) is gaussian noise with mean zero and
%                       covariance Q(k-1)
% The observation model is given by:
% z(k) = h(k)(x(k),v(k) = H(k)x(k) + v(k)
%                  h(k)(.) is a linear (in KF case) or non-linear) function.
%                  Meaning the observation vector z is a linear function
%                  of the state vector, and this linear relationship is
%                  represented by pre-multiplication of observation
%                  matrix H. There is also gaussian measurement
%                  noise v.
% where v ~ N(0,R(k)) meaning v(k) is gaussian noise with mean zero and
%                     covariance R(k)
%
% Version 1.0, January, 2019
%
% This function was written by Nathanael L. Baisa
% email: nathanaellmss@gmail.com
%

function [m_update, P_update] = KalmanFilter(z, model, m_update, P_update)
% Inputs are measurement (z), model (F, Q, H, R), and m_update and P_update 
% at time k-1.
% Outputs are m_update and P_update at time k.

% Prediction (time update)
m_predict = model.F * m_update;
P_predict = model.Q + model.F*P_update*model.F';

% Compute the covariance (S) of the innovation term z(k)-H(k)*m(k|k-1) and
% then then Kalman gain (K)
S = model.H*P_predict*model.H' + model.R;  iS = inv(S);
K = P_predict*model.H'*iS;
% S  = model.R+model.H*P_predict*model.H'; Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs);  iS= inv_sqrt_S*inv_sqrt_S';
% K  = P_predict*model.H'*iS;

% Correction (measurement update)
m_update = m_predict + K*(z - model.H*m_predict);
P_update = P_predict - K*model.H*P_predict;

end


