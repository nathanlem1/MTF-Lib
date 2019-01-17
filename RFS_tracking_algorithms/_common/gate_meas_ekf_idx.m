function valid_idx= gate_meas_ekf(z,gamma,model,m,P)

valid_idx = [];
zlength = size(z,2); if zlength==0, z_gate= []; return; end
plength = size(m,2);

for j=1:plength
        [C_ekf,U_ekf]= ekf_update_mat(model,m(:,j));
        Sj= U_ekf*model.R*U_ekf' + C_ekf*P(:,:,j)*C_ekf'; Sj= (Sj+ Sj')/2; 
        Vs= chol(Sj); det_Sj= prod(diag(Vs))^2; inv_sqrt_Sj= inv(Vs);
        iSj= inv_sqrt_Sj*inv_sqrt_Sj'; 
        nu= z- repmat(gen_observation_fn(model,m(:,j),zeros(size(model.D,2),1)),[1 zlength]);
        dist= sum((inv_sqrt_Sj'*nu).^2);
        valid_idx= unique_faster([ valid_idx find( dist < gamma )]);
end
valid_idx=valid_idx(:)';