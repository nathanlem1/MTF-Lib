% This function implements unscented Kalman filter (UKF) in a clear and 
% understandable manner. The tutorial for the UKF is given in the paper
% [1]. A nearly-constant turn (CT) model used in [2] is used in this 
% implementation.
% 1. E.A. Wan and R. Van Der Merwe, 'The unscented Kalman filter for nonlinear estimation'
% 2. Ba-Ngu Vo, and Wing-Kin Ma, 'The Gaussian Mixture Probability Hypothesis
% Density Filter'
%
% The state dynamics evolves according to the following equations:
% x(k) = f(k)(x(k-1), u(k-1), w(k-1)) = f(k)(x(k-1), w(k-1)) [assume u(k-1) =0] 
%                  f(k)(.) is a non-linear function. A local linear 
%                  approximation around the current mean m(k-1) at each 
%                  time should be used to get a Jacobian matrices F(k) 
%                  (partial derivatives of f(k)(.) with respect to x) and 
%                  W(k)(partial derivatives of f(k)(.) with respect to 
%                  Gaussian process noise w).
% where w ~ N(0,Q(k-1)) meaning w(k-1) is gaussian noise with mean zero and
%                       covariance Q(k-1)
% The observation model is given by:
% z(k) = h(k)(x(k),v(k)
%                  h(k)(.) is a non-linear function. A local linear 
%                  approximation around the current mean m(k|k-1) at each 
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

function [m_update, P_update] = unscentedKalmanFilter(z, model, m_update, P_update, alpha, kappa, beta)
    % Inputs are measurement (z), model (F, Q, H, R), and m_update and P_update 
    % at time k-1.
    % Outputs are m_update and P_update at time k.
    
    % Prediction (time update)
    [X_ukf,u]= unscentedTransform( [m_update; zeros(model.v_dim,1) ], blkdiag(P_update,model.Q), alpha, kappa );
    X_pred= gen_newstate_fn( model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.v_dim),:) );

    m_predict = X_pred*u(:);
    X_temp= X_pred- repmat(m_predict,[1 length(u)]);
    u(1)= u(1)+(1-alpha^2+beta);
    P_predict= X_temp*diag(u)*X_temp';    
    
    % Correction (measurement update)
    [X_ukf,u]= unscentedTransform( [m_predict; zeros(model.w_dim,1) ], blkdiag(P_predict,model.R), alpha, kappa );
    Z_pred= gen_observation_fn( model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.w_dim),:) );

    eta= Z_pred*u(:);

    S_temp= Z_pred- repmat(eta,[1 length(u)]);
    u(1)= u(1)+(1-alpha^2+beta);
    S= S_temp*diag(u)*S_temp';
    Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';

    G_temp= X_ukf(1:model.x_dim,:)- repmat(m_predict,[1 length(u)]);
    G= G_temp*diag(u)*S_temp';
    K  = G*iS;

    m_update = repmat(m_predict,[1 size(z,2)]) + K*(z-repmat(eta,[1 size(z,2)]));
    P_update = P_predict- G*iS*G';
   
end


