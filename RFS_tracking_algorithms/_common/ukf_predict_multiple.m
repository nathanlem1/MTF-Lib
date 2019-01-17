function [m_predict,P_predict] = ukf_predict_multiple(model,m,P,alpha,kappa,beta)      

plength= size(m,2);

m_predict = zeros(size(m));
P_predict = zeros(size(P));

for idxp=1:plength
    [m_temp,P_temp] = ukf_predict_single(model,m(:,idxp),P(:,:,idxp),alpha,kappa,beta);
    m_predict(:,idxp) = m_temp;
    P_predict(:,:,idxp) = P_temp;
end

function [m_predict,P_predict] = ukf_predict_single(model,m,P,alpha,kappa,beta)

[X_ukf,u]= ut( [m; zeros(model.v_dim,1) ], blkdiag(P,model.Q), alpha, kappa );
X_pred= gen_newstate_fn( model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.v_dim),:) );

m_predict = X_pred*u(:);
X_temp= X_pred- repmat(m_predict,[1 length(u)]);
u(1)= u(1)+(1-alpha^2+beta);
P_predict= X_temp*diag(u)*X_temp';
