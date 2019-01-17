% Copyright (c) 2012, Derek O'Connor
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
%
function [dist,path,pred]= BFMSpathwrap(ncm,source,destination)
% Wrapper to convert output to MATLAB 'graphshortestpath' format (by BT Vo)
[p,D,~]= BFMSpathOT(ncm,source);
dist= D(destination);
pred= p';

if isinf(dist)
    path=[];
else
    path= destination;
    while path(1) ~= source
        path= [pred(path(1)) path];
    end
end
end


% ==================================================================== 
  function [p,D,iter] = BFMSpathOT(G,r)
% --------------------------------------------------------------------
% Basic form of the Bellman-Ford-Moore Shortest Path algorithm
% Assumes G(N,A) is in sparse adjacency matrix form, with |N| = n, 
% |A| = m = nnz(G). It constructs a shortest path tree with root r which 
% is represented by an vector of parent 'pointers' p, along with a vector
% of shortest path lengths D.
% Complexity: O(mn)
% Derek O'Connor, 19 Jan, 11 Sep 2012.  derekroconnor@eircom.net
% -------------------------------------------------------------------- 
% Unlike the original BFM algorithm, this does an optimality test on the
% SP Tree p which may greatly reduce the number of iters to convergence.
% USE: 
% n=10^6; G=sprand(n,n,5/n); r=1; format long g;
% tic; [p,D,iter] = BFMSpathOT(G,r);toc, disp([(1:10)' p(1:10) D(1:10)]);
% WARNING: 
% This algorithm performs well on random graphs but may perform 
% badly on real problems.
% --------------------------------------------------------------------


[m,n,p,D,tail,head,W] = Initialize(G);
p(r)=0; D(r)=0;                  % Set the root of the SP tree (p,D)
for iter = 1:n-1                   % Converges in <= n-1 iters if no 
    optimal = true;                % negative cycles exist in G
    for arc = 1:m                  % O(m) optimality test. 
        u = tail(arc); 
        v = head(arc);
        duv = W(arc);
        if D(v) > D(u) + duv;      %
           D(v) = D(u) + duv;      % Sp Tree not optimal: Update (p,D) 
           p(v) = u;               %
           optimal = false;
        end
    end % for arc
    if optimal
        return                     % SP Tree p is optimal;
    end
end % for iter

%---------------END BFMSpathOT--------------------------------------
 end


function [m,n,p,D,tail,head,W] = Initialize(G)
%
% Transforms the sparse matrix G into the list-of-arcs form
% and intializes the shortest path parent-pointer and distance
% arrays, p and D.
% Derek O'Connor, 21 Jan 2012

[tail,head,W] = find(G); % Get arc list {u,v,duv, 1:m} from G.
[~,n] = size(G); 
m = nnz(G);
p(1:n,1) = 0;            % Shortest path tree of parent pointers
D(1:n,1) = Inf;          % Sh. path distances from node i=1:n to 
                         % the root of the sp tree
                
% NOTE: After a shortest path function calls this function
%       it will need to set D(r) = 0, for a given r, the root
%       of the shortest path tree it is about to construct.
end
    
    


