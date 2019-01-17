function [X,w]= ut(m,P,alpha,kappa)

% standard unscented transformation - note [choose alpha=1 => =lambda=kappa; giving offset of beta for first cov weight]

n_x= length(m);
lambda= alpha^2*(n_x+kappa) - n_x;
Psqrtm= chol((n_x+lambda)*P)';
X= repmat(m,[1 2*n_x+1])+ [ zeros(n_x,1) -Psqrtm Psqrtm ];
w= [ lambda 0.5*ones(1,2*n_x) ]/(n_x+lambda);
