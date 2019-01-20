function X= gen_newstate_fn(model,Xd,V)

%nonlinear state space equation (CT model)

if ~isnumeric(V)
    if strcmp(V,'noise')
        V= model.B*randn(size(model.B,2),size(Xd,2));
    elseif strcmp(V,'noiseless')
        V= zeros(size(model.B,1),size(Xd,2));
    end
end

if isempty(Xd)
    X= [];
else %modify below here for user specified transition model
    X= zeros(size(Xd));
    %-- short hand
    L= size(Xd,2);
    T= model.T; 
    omega= Xd(5,:);
    tol= 1e-10;
    %-- pre calcs
    sin_omega_T= sin(omega*T);
    cos_omega_T= cos(omega*T);
    a= T*ones(1,L); b= zeros(1,L);
    idx= find( abs(omega) > tol );
    a(idx)= sin_omega_T(idx)./omega(idx);
    b(idx)= (1-cos_omega_T(idx))./omega(idx);
    %--  x/y pos/vel
    X(1,:)= Xd(1,:)+ a.*Xd(2,:)- b.*Xd(4,:);
    X(2,:)= cos_omega_T.*Xd(2,:)- sin_omega_T.*Xd(4,:);
    X(3,:)= b.*Xd(2,:) + Xd(3,:)+ a.*Xd(4,:);
    X(4,:)= sin_omega_T.*Xd(2,:)+ cos_omega_T.*Xd(4,:);
    %-- turn rate
    X(5,:)= Xd(5,:);
    %-- add scaled noise 
    X= X+ model.B2*V;
end