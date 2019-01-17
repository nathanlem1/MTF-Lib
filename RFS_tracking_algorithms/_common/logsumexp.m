function logsum = logsumexp(w)

%performs log-sum-exp trick to avoid numerical underflow
%input:  w weight vector assumed already log transformed
%output: log(sum(exp(w)))

if all(w==-inf)
    logsum= -inf;
    return;
end

[val,idx] = max(w);
logsum = log(sum(exp(w-val))) + val;

