function [weight count edges mid loc] = histcn(w, X, edges)%varargin)
% function [count edges mid loc] = histcn(w, X, {edge1, edge2, ..., edgeN}) %edge1, edge2, ..., edgeN)
%
% Purpose: compute n-dimensional histogram
%
% INPUT
%   - w: is (M x 1) array, represents weight of each data point
%   - X: is (M x N) array, represents M data points in R^N
%   - edgek: are the bin vectors on dimension k, k=1...N.
%     If it is a scalar (Nk), the bins will be the linear subdivision of
%     the data on the range [min(X(:,k)), max(X(:,k))] into Nk
%     sub-intervals
%     If it's empty, a default of 32 subdivions will be used
% OUTPUT
%   - count: n-dimensional array count of X on the bins, i.e.,
%         count(i1,i2,...,iN) = cardinal of X such that
%                  edge1(i1) <= X(:,i1) < edge1(i1)+1 and
%                       ...
%                  edgeN(iN) <= X(:,iN) < edgeN(iN)+1
%   - edges: (1 x N) cell, each provides the effective edges used in the
%     respective dimension
%   - mid: (1 x N) cell, provides the mid points of the cellpatch used in
%     the respective dimension
%   - loc: (M x N) array, index location of X in the bins. Points have out
%     of range coordinates will have zero at the corresponding dimension.
%
% Usage example:
%       M = 1e5;
%       N = 3;
%       X = randn(M,N);
%       [N edges mid loc] = histcn(X);
%       imagesc(mid{1:2},N(:,:,ceil(end/2)))
% 
% Bruno Luong: <brunoluong@yahoo.com>
% Last update: 25/April/2009

if ndims(X)>2
    error('histcn: X requires to be an (M x N) array of M points in R^N');
end
DEFAULT_NBINS = 32;

% Get the dimension
nd = size(X,2);
%edges = varargin;
if nd<length(edges)
    nd = length(edges); % waisting CPU time warranty
else
    edges(end+1:nd) = {DEFAULT_NBINS};
end

% Allocation of array loc: index location of X in the bins
loc = zeros(size(X));
sz = zeros(1,nd);
% Loop in the dimension
for d=1:nd
    ed = edges{d};
    Xd = X(:,d);
    if isempty(ed)
        ed = DEFAULT_NBINS;
    end
    if isscalar(ed) % automatic linear subdivision
        ed = linspace(min(Xd),max(Xd),ed+1);
    end
    edges{d} = ed;
    % Call histc on this dimension
    [dummy loc(:,d)] = histc(Xd, ed, 1);
    sz(d) = length(ed)-1;
end % for-loop

% Clean
clear dummy

% This is need for seldome points that hit the right border
sz = max([sz; max(loc,[],1)]);

% Compute the mid points
mid = cellfun(@(e) 0.5*(e(1:end-1)+e(2:end)), edges, ...
              'UniformOutput', false);
          
% Count for points where all coordinates are falling in a corresponding bins
if nd==1
    sz = [sz 1]; % Matlab doesn't know what is one-dimensional array!
end
count = accumarray(loc(all(loc>0, 2),:),1,sz);
weight = accumarray(loc(all(loc>0, 2),:),w(all(loc>0, 2)),sz);

return

end % histcn
