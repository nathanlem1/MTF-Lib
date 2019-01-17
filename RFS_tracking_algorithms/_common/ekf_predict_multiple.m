function [m_predict,P_predict] = ekf_predict_multiple(model,m,P)    

plength= size(m,2);

m_predict = zeros(size(m));
P_predict = zeros(size(P));

for idxp=1:plength
    [m_temp,P_temp] = ekf_predict_single(model,m(:,idxp),P(:,:,idxp));
    m_predict(:,idxp) = m_temp;
    P_predict(:,:,idxp) = P_temp;
end

function [m_predict,P_predict] = ekf_predict_single(model,m,P)

m_predict = gen_newstate_fn(model,m,'noiseless');
[F_ekf,G_ekf]= ekf_predict_mat(model,m);                % user specified function for application
P_predict= G_ekf*model.Q*G_ekf' + F_ekf*P*F_ekf';
P_predict= (P_predict+P_predict')/2;                    % addition step to avoid numerical problem