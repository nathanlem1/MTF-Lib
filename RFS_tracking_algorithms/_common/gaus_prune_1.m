function [w_new,x_new,P_new,u_new,v_new]= gaus_prune_1(w,x,P,u,v,elim_threshold)

idx= find( w > elim_threshold );
w_new= w(idx);
x_new= x(:,idx);
P_new= P(:,:,idx);
u_new= u(idx);
v_new= v(idx);