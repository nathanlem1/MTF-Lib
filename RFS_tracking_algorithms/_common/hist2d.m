function [Z,x,y]= hist2d(C,no_of_bins,limit)
% compute histogram for 2D data
% C- input data with dimension 2 x data_len 
% no_of_bins= [ no of bins for x-axis, no of bins for y-axis ]
% limit= [ lower limit for x-axis, upper limit for x-axis
%          lower limit for y-axis, upper limit for y-axis

data_len= size(C,2);
Z= sparse(no_of_bins(1),no_of_bins(2));

[point_list,value_list,scan_points]= ...
    histnd(C,no_of_bins,limit);
Z= sparse(point_list(1,:),point_list(2,:),value_list,no_of_bins(1),no_of_bins(2));
x= scan_points{1}; y= scan_points{2};
