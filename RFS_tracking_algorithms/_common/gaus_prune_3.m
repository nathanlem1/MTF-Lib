function [w_new,x_new,P_new,u_new,v_new,wc_new,uc_new,vc_new]= gaus_prune_3(w,x,P,u,v,wc,uc,vc,elim_threshold)

idx1= find( w > elim_threshold );
w_new= w(idx1);
x_new= x(:,idx1);
P_new= P(:,:,idx1);
u_new= u(idx1);
v_new= v(idx1);

idx2= find(wc > elim_threshold);
wc_new= wc(idx2);
uc_new= uc(idx2);
vc_new= vc(idx2);