function [point_list,value_list,scan_points,pt_list2]= histnd(C,no_of_bins,limit)

[n,data_len]= size(C);
range= limit*[ -1; 1 ];    %range for each dim.
% compute scan points
scan_points= cell(n,1); delta= zeros(n,1);
for i=1:n
    delta(i)= range(i)/(no_of_bins(i)-1);
    scan_points{i}= limit(i,1):delta(i):limit(i,2);
end;
% quantization and truncation with data
Ct= C- limit(:,1)*ones(1,data_len);
Ct= round(diag(1./delta)*Ct)+1;
[I,J]= find(Ct<=0); Ct(I,J)= 1;
for i=1:n
    I= find( Ct(i,:)>no_of_bins(i) );
    Ct(i,I)= no_of_bins(i);
end;
% 
len= 0; point_list= zeros(n,1); value_list= 0;
if n > 1,
    while ~isempty(Ct);
        bt= Ct(:,1); Ct= Ct(:,2:size(Ct,2));
        I= find( sum(abs(Ct-bt*ones(1,size(Ct,2))))~= 0 );
        len= len+1;
        point_list(:,len)= bt;
        value_list(:,len)= 1+ size(Ct,2)- length(I);
        Ct= Ct(:,I);
    end;
else
    while ~isempty(Ct);
        bt= Ct(:,1); Ct= Ct(:,2:size(Ct,2));
        I= find( Ct-bt*ones(1,size(Ct,2))~= 0 );
        len= len+1;
        point_list(:,len)= bt;
        value_list(:,len)= 1+ size(Ct,2)- length(I);
        Ct= Ct(:,I);
    end;
end;
%
if nargout==4,
    pt_list2= diag(delta)*(point_list-1)+ limit(:,1)*ones(1,size(point_list,2));
end;