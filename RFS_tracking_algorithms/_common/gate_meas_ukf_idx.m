function valid_idx= gate_meas_ukf(z,gamma,model,m,P,alpha,kappa,beta)

valid_idx = [];
zlength = size(z,2); if zlength==0, z_gate= []; return; end
plength = size(m,2);

if zlength>0
for j=1:plength
        [X_ukf,u]= ut( [m(:,j); zeros(model.w_dim,1) ], blkdiag(P(:,:,j),model.R), alpha, kappa );
        Z_pred= gen_observation_fn( model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.w_dim),:) );
        eta= Z_pred*u(:); Sj_temp= Z_pred- repmat(eta,[1 length(u)]); u(1)= u(1)+(1-alpha^2+beta);
        Sj= Sj_temp*diag(u)*Sj_temp';
        Vs= chol(Sj); det_Sj= prod(diag(Vs))^2; inv_sqrt_Sj= inv(Vs);
        iSj= inv_sqrt_Sj*inv_sqrt_Sj'; 
        nu= z- repmat(gen_observation_fn(model,m(:,j),zeros(size(model.D,2),1)),[1 zlength]);
        dist= sum((inv_sqrt_Sj'*nu).^2);
        valid_idx= unique_faster([ valid_idx find( dist < gamma )]);;
end
end
valid_idx=valid_idx(:)';