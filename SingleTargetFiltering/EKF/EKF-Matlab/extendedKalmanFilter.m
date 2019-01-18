% This function implements Kalman filter (KF) in a clear and understandable
% manner. The tutorial for the KF is geiven in the following papers:
% 1. Greg Welch and Gary Bishop, 'An Introduction to the Kalman Filter'
%
% The state dynamics evolves according to the following equations:
% x(k) = f(k)(x(k-1), u(k-1), w(k-1)) = f(k)(x(k-1), w(k-1)) [assume u(k-1) =0] 
%                  f(k)(.) is a non-linear function. A local linear 
%                  approximation around the current mean m(k?1) at each 
%                  time should be used to get a Jacobian matrices F(k) 
%                  (partial derivatives of f(k)(.) with respect to x) and 
%                  W(k)(partial derivatives of f(k)(.) with respect to 
%                  Gaussian process noise w).
% where w ~ N(0,Q(k-1)) meaning w(k-1) is gaussian noise with mean zero and
%                       covariance Q(k-1)
% The observation model is given by:
% z(k) = h(k)(x(k),v(k)
%                  h(k)(.) is a non-linear function. A local linear 
%                  approximation around the current mean m(k|k?1) at each 
%                  time should be used to get a Jacobian matrices H(k) 
%                  (partial derivatives of h(k)(.) with respect to x) and 
%                  V(k)(partial derivatives of h(k)(.) with respect to 
%                  Gaussian measurement noise v).
% where v ~ N(0,R(k)) meaning v(k) is gaussian noise with mean zero and
%                     covariance R(k)
%
% Version 1.0, January, 2019
%
% This function was written by Nathanael L. Baisa
% email: nathanaellmss@gmail.com
%

function [m_update, P_update] = extendedKalmanFilter(z, model, m_update, P_update)
    % Inputs are measurement (z), model (F, Q, H, R), and m_update and P_update 
    % at time k-1.
    % Outputs are m_update and P_update at time k.

    % Prediction (time update)
    m_predict = gen_newstate_fn(model,m_update,'noiseless');  
    [F_ekf,W_ekf]= ekf_predict_mat(model,m_update);         % user specified function for application (for Jacobian matrices)
    P_predict = W_ekf*model.Q*W_ekf' + F_ekf*P_update*F_ekf';
    P_predict= (P_predict+P_predict')/2;                    % addition step to avoid numerical problem

    % Compute the covariance (S) of the innovation term z(k)-H(k)*h(k)(m(k|k-1), 0) and
    % then then Kalman gain (K)
    [H_ekf,V_ekf]= ekf_update_mat(model,m_predict);                 % user specified function for application (for Jacobian matrices)
    S =  V_ekf*model.R*V_ekf' + H_ekf*P_predict*H_ekf';  S = (S+ S')/2;   % addition step to avoid numerical problem
    iS = inv(S);
    K = P_predict*H_ekf'*iS;
    % S  = V_ekf*model.R*V_ekf' + H_ekf*P_predict*H_ekf'; S = (S+ S')/2;   % addition step to avoid numerical problem
    % Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs);  iS= inv_sqrt_S*inv_sqrt_S';
    % K  = P_predict*H_ekf'*iS;

    % Correction (measurement update)
    h_eta = gen_observation_fn(model,m_predict,'noiseless');
    m_update = m_predict + K*(z - h_eta);
    P_update = P_predict - K*H_ekf*P_predict;

end


