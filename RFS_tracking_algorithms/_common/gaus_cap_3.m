function [w_new,x_new,P_new,u_new,v_new,wc_new,uc_new,vc_new]= gaus_cap_3(w,x,P,u,v,wc,uc,vc,max_number)

if length(w) > max_number
    [notused,idx1]= sort(w,1,'descend');
    w_new= w(idx1(1:max_number)); w_new = w_new * (sum(w)/sum(w_new));
    x_new= x(:,idx1(1:max_number));
    P_new= P(:,:,idx1(1:max_number));
    u_new= u(idx1(1:max_number));
    v_new= v(idx1(1:max_number));
else
    x_new = x;
    P_new = P;
    w_new = w;
    u_new = u;
    v_new = v;
end
if length(wc) > max_number
    [notused,idx2]= sort(wc,1,'descend');
    wc_new= wc(idx2(1:max_number)); wc_new = wc_new * (sum(wc)/sum(wc_new));
    uc_new= uc(idx2(1:max_number));
    vc_new= vc(idx2(1:max_number));
else
    wc_new= wc;
    uc_new = uc;
    vc_new = vc;
end
