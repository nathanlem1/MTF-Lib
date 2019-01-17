function [ x_c, I_idx ]= our_kmeans( X, w, iter )
% custom k-means clustering
% X- particles for clustering
% w- weights of particles
% iter- max no. of adjustment iters
% x_c,I_idx - cluster centres and corresponding particles indices

x_c= []; I_idx= [];
[dim_x,Ns]= size(X);

%--- initialization step
big_I= 1:Ns; el= 0;
while ~isempty(big_I),
    el= el+ 1;
    len= length(big_I); i= round(rand*(len-1))+1;
    xc_cur= X(:,big_I(i));
    d= sum((repmat(xc_cur,[1 len]) - X(:,big_I)).^2);
    [notused,idx_sort]= sort(d);
    w_acc= 0; I_idx_temp= [];
    for k=1:len
        w_acc= w_acc+ w(big_I(idx_sort(k)));
        if w_acc <=1,
            I_idx_temp= [ I_idx_temp big_I(idx_sort(k)) ];
        else
            break;
        end;
    end;
    big_I= setdiff(big_I,I_idx_temp);
    I_idx{el}= I_idx_temp;
    x_c(:,el)= X(:,I_idx{el})*w(I_idx{el}) / sum(w(I_idx{el}));
end;
L= el;  %no. of clusters found

%--- iterative adjustment of the clusters
for m=1:iter
    I_idx= cell(L,1); w_acc= zeros(L,1);
    for i=1:Ns
        d= sum((repmat(X(:,i),[1 L])- x_c).^2);
        [notused,idx_sort]= sort(d);
        k= 1;
        while w_acc(idx_sort(k))+ w(i) > 1,
            k= k+1;
            if k > L, %that is a rare case, increase L
                disp('has to increase size');
                L= L+1;
                x_c(:,L)= X(:,i); 
                idx_sort(k)= L; w_acc(L)= 0; I_idx{L+1}= [];
                break;
            end;
        end;
        w_acc(idx_sort(k))= w_acc(idx_sort(k))+ w(i);
        I_idx{idx_sort(k)}= [ I_idx{idx_sort(k)} i ];
    end;
    for el=1:L
        x_c(:,el)= X(:,I_idx{el})*w(I_idx{el}) / sum(w(I_idx{el}));
    end;
end;
            