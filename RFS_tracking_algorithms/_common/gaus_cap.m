function [w_new,x_new,P_new]= gaus_cap(w,x,P,max_number)

if length(w) > max_number
    [notused,idx]= sort(w,1,'descend');
    w_new= w(idx(1:max_number)); w_new = w_new * (sum(w)/sum(w_new));
    x_new= x(:,idx(1:max_number));
    P_new= P(:,:,idx(1:max_number));
else
    x_new = x;
    P_new = P;
    w_new = w;
end
