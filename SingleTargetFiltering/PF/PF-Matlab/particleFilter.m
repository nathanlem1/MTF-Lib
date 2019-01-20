% This function implements particle filter (PF) in a clear and 
% understandable manner.  The tutorial for the PF is given in the paper
% [1]. A nearly-constant turn (CT) model used in [2] is used in this 
% implementation.
% 1. M.S. Arulampalam, S. Maskell, N. Gordon and T. Clapp, 'A tutorial on 
%   particle filters for online nonlinear/non-Gaussian Bayesian tracking'
% 2. Ba-Ngu Vo, and Wing-Kin Ma, 'The Gaussian Mixture Probability Hypothesis
%   Density Filter'
%
% The state dynamics evolves according to the following equations:
% x(k) = f(k)(x(k-1), u(k-1), w(k-1)) = f(k)(x(k-1), w(k-1)) [assume u(k-1) =0] 
%                  f(k)(.) is a non-linear function with non-Gaussian process
%                  noise w.
%
% The observation model is given by:
% z(k) = h(k)(x(k),v(k)
%                  h(k)(.) is a non-linear function with non-Gaussian 
%                  measurement noise v.
%
% This code implements the sampling importance resampling (SIR) particle 
% filter. The key idea is to represent the required posterior density
% function by a set of random samples with associated weights and to 
% compute estimates based on these samples and weights. As the number
% of samples becomes very large, this Monte Carlo (MC) characterization
% becomes an equivalent representation to the usual functional description
% of the posterior pdf, and the PF approaches the optimal Bayesian estimate.
%
% Version 1.0, January, 2019
%
% This function was written by Nathanael L. Baisa
% email: nathanaellmss@gmail.com
%

function [x_update, w_update] = particleFilter(z, model, x_update, w_update, filter)
    % Inputs are measurement (z), model (B, D), and x_update and w_update 
    % at time k-1.
    % Outputs are x_update and w_update at time k.
    
    % Prediction (time update)
    x_predict = gen_newstate_fn(model,x_update,'noise');
    w_predict = w_update; 

    % Correction (measurement update)   
    meas_likelihood =  compute_likelihood(model, z, x_predict)';
    w_update = meas_likelihood.*w_predict;
    x_update = x_predict; 
    
    % Normalize weights
    w_update = w_update/sum(w_update);   
        
    % Resampling
    idx= randsample(length(w_update),filter.J_max,true,w_update); %idx= resample(w_update,filter.J_max);
    w_update= ones(filter.J_max,1)/filter.J_max;
    x_update= x_update(:,idx);
    
end


