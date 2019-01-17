function [qz_update,m_update,P_update] = ekf_update_multiple(z,model,m,P)

plength= size(m,2);
zlength= size(z,2);

qz_update= zeros(plength,zlength);
m_update = zeros(model.x_dim,plength,zlength);
P_update = zeros(model.x_dim,model.x_dim,plength);

for idxp=1:plength
        [qz_temp,m_temp,P_temp] = ekf_update_single(z,model,m(:,idxp),P(:,:,idxp));
       qz_update(idxp,:)   = qz_temp;
        m_update(:,idxp,:) = m_temp;
        P_update(:,:,idxp) = P_temp;
end

function [qz_temp,m_temp,P_temp] = ekf_update_single(z,model,m,P)

eta = gen_observation_fn(model,m,'noiseless');

[H_ekf,U_ekf]= ekf_update_mat(model,m);                 % user specified function for application
S= U_ekf*model.R*U_ekf'+H_ekf*P*H_ekf'; S= (S+ S')/2;   % addition step to avoid numerical problem
Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';

K  = P*H_ekf'*iS;

qz_temp = exp(-0.5*size(z,1)*log(2*pi) - 0.5*log(det_S) - 0.5*dot(z-repmat(eta,[1 size(z,2)]),iS*(z-repmat(eta,[1 size(z,2)]))))';
m_temp = repmat(m,[1 size(z,2)]) + K*(z-repmat(eta,[1 size(z,2)]));
P_temp = (eye(size(P))-K*H_ekf)*P;



