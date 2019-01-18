function [F,G]= ekf_predict_mat(model,mu_old)

tol= 1e-6;
omega= mu_old(5);

T= model.T;
sin_omega_T= sin(omega*T);
cos_omega_T= cos(omega*T);

a= T; b= 0;
if abs(omega) > tol
    a= sin_omega_T/omega;
    b= (1-cos_omega_T)/omega;
end;
A= [ 1 a 0 -b; ...
     0 cos_omega_T 0 -sin_omega_T; ...
     0 b 1 a; ...
     0 sin_omega_T 0 cos_omega_T ];

c= 0; d= T^2/2;
if abs(omega) > tol
    c= (omega*T*cos_omega_T - sin_omega_T )/omega^2;
    d= ( omega*T*sin_omega_T - 1 + cos_omega_T )/omega^2;
end;
dA= [ 0 c 0 -d; ...
        0 -T*sin_omega_T 0 -T*cos_omega_T; ...
        0 d 0 c; ...
        0 T*cos_omega_T 0 -T*sin_omega_T ];

F= [ A dA*mu_old(1:4); ...
        zeros(1,4) 1 ];

G= model.B2;
    

    