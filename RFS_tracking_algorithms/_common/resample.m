function resample_idx= resample(w,L)
% particles resampling 
% w- weights with sum(w)= 1
% L- no. of samples desired
% resample_idx- indices for the resampled particles

w=w/sum(w);
resample_idx= [];
[notused,sort_idx]= sort(-w);   %sort in descending order
rv= rand(L,1);
i= 0; threshold= 0;
while ~isempty(rv),
    i= i+1; threshold= threshold+ w(sort_idx(i));
    rv_len= length(rv);
    idx= find(rv>threshold); 
    resample_idx= [ resample_idx; sort_idx(i)*ones(rv_len-length(idx),1) ];
    rv= rv(idx);
end;