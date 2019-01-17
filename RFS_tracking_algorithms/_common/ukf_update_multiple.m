function [qz_update,m_update,P_update] = ukf_update_multiple(z,model,m,P,alpha,kappa,beta)

plength= size(m,2);
zlength= size(z,2);

qz_update= zeros(plength,zlength);
m_update = zeros(model.x_dim,plength,zlength);
P_update = zeros(model.x_dim,model.x_dim,plength);

for idxp=1:plength
        [qz_temp,m_temp,P_temp] = ukf_update_single(z,model,m(:,idxp),P(:,:,idxp),alpha,kappa,beta);
       qz_update(idxp,:)   = qz_temp;
        m_update(:,idxp,:) = m_temp;
        P_update(:,:,idxp) = P_temp;
end

function [qz_temp,m_temp,P_temp] = ukf_update_single(z,model,m,P,alpha,kappa,beta)

[X_ukf,u]= ut( [m; zeros(model.w_dim,1) ], blkdiag(P,model.R), alpha, kappa );
Z_pred= gen_observation_fn( model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.w_dim),:) );

eta= Z_pred*u(:);

S_temp= Z_pred- repmat(eta,[1 length(u)]);
u(1)= u(1)+(1-alpha^2+beta);
S= S_temp*diag(u)*S_temp';
Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';

G_temp= X_ukf(1:model.x_dim,:)- repmat(m,[1 length(u)]);
G= G_temp*diag(u)*S_temp';
K  = G*iS;

qz_temp = exp(-0.5*size(z,1)*log(2*pi) - 0.5*log(det_S) - 0.5*dot(z-repmat(eta,[1 size(z,2)]),iS*(z-repmat(eta,[1 size(z,2)]))))';
m_temp = repmat(m,[1 size(z,2)]) + K*(z-repmat(eta,[1 size(z,2)]));
P_temp = P- G*iS*G';



