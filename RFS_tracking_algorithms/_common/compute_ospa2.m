function result = compute_ospa2(X,Y,c,p,wl)

% This is the MATLAB code for the implementation OSPA(2) metric proposed in
% M. Beard, B.-T. Vo, and B.-N. Vo, "Performance Evaluation for Large-Scale Multi-Target Tracking Algorithms," Proc. 21st IEEE Intl. Conf. Information Fusion, Jul. 2018, Cambridge, UK.
% http://ba-ngu.vo-au.com/vo/BVV_OSPA2_FUSION18.pdf
% and
% M. Beard, B.-T. Vo, and B.-N. Vo, "A Solution for Large-Scale Multi-Object Tracking," arXiv preprint, arXiv:1804.06622, Apr. 2018.
% https://arxiv.org/abs/1804.06622
% based on the OSPA metric proposed in
% D. Schuhmacher, B.-T. Vo, and B.-N. Vo, "A consistent metric for performance evaluation in multi-object filtering," IEEE Trans. Signal Processing, Vol. 56, No. 8 Part 1, pp. 3447– 3457, 2008.
% http://ba-ngu.vo-au.com/vo/SVV08_OSPA.pdf
%
% ---BibTeX entry
% @inproceedings{OSPA2,
% author={M. Beard and B.-T. Vo and B.-N. Vo},
% booktitle = {Proc. 21st IEEE Intl. Conf. Information Fusion},
% title={Performance Evaluation for Large-Scale Multi-Target Tracking Algorithms},
% month= {Jul},
% year={2018},
% location= {Cambridge, UK}}
%
% @ARTICLE{LST,
% author={M. Beard and B.-T. Vo and B.-N. Vo},
% title = "{A Solution for Large-scale Multi-object Tracking}",
% journal = {ArXiv e-prints},
% archivePrefix = "arXiv",
% eprint = {1804.06622},
% year = 2018,
% month = apr}
%
% @ARTICLE{OSPA,
% author={D. Schuhmacher and B.-T. Vo and B.-N. Vo},
% journal={IEEE Transactions on Signal Processing},
% title={A Consistent Metric for Performance Evaluation of Multi-Object Filters},
% year={2008},
% month={Aug},
% volume={56},
% number={8},
% pages={3447-3457}}  
% ---
%
% Inputs:
% X    DxTxN array, where D is the state dimension, T is the number of
%      time steps, and N is the number of objects. NaN is used to
%      indicate when an object is not present.
%
% Y    DxTxM array, where D and T must match the dimensions of X, and M
%      is the number of tracks.
%
% c    Cutoff parameter (used for both inner and outer OSPA)
%
% p    Order parameter (used for both inner and outer OSPA)
%
% wl   Size of the moving window. For example, at time t, the metric
%      will be computed over the window [t-wl+1:t].
%
% Output:
% result   3xT array, where result(1,t) is the OSPA(2) at time t,
%          result(2,t) is the localisation component at time t, and
%          result(3,t) is the cardinality component at time t.
%
%
  
  if (size(X,1) ~= size(Y,1)) || (size(X,2) ~= size(Y,2))
    error('Dimensions of X and Y are inconsistent');
  end

  if ~isnumeric(c) || ~isscalar(c) || (c <= 0)
    error('c must be a positive scalar');
  end

  if ~isnumeric(p) || ~isscalar(p) || (p <= 0)
    error('p must be a positive scalar');
  end

  wl = floor(wl);
  if ~isnumeric(wl) || ~isscalar(wl) || (wl <= 0)
    error('Window length must be a positive integer');
  end

  eval_idx = 1:size(X,2);
  truncated = true; 
  win_off = (-wl+1):0;
  
  num_x = size(X,3);
  num_y = size(Y,3);
  num_step = size(X,2);
  num_eval = length(eval_idx);
  
  result = zeros(3,num_eval);
  
  % First, for each time index compute the matrix of inter-point
  % distances and the track existence flags
  
  distances = zeros(num_x,num_y,num_step);
  x_exists = false(num_x,num_step);
  y_exists = false(num_y,num_step);
  
  for i = 1:num_step
      
    % Compute distance between every pair of points
    x = permute(X(:,i,:),[1 3 2]);
    y = permute(Y(:,i,:),[1 3 2]);
    d = permute(sum(abs(bsxfun(@minus,permute(y,[1 3 2]),x)).^p,1),[2 3 1]);
    
    % Compute track existence flags
    x_exists(:,i) = ~isnan(x(1,:));
    y_exists(:,i) = ~isnan(y(1,:));
        
    % Distance between an empty and non-empty state
    one_exists = bsxfun(@xor,x_exists(:,i),y_exists(:,i)');
    d(one_exists) = c^p;     
    
    % Find times when neither object exists
    neither_exists = bsxfun(@and,~x_exists(:,i),~y_exists(:,i)');
    if truncated
      % Truncated window, exclude times when neither objects exists
      d(neither_exists) = NaN;
    else
      % Full window, distance between empty states is zero
      d(neither_exists) = 0;
    end
    
    % Store the distance matrix for this step
    distances(:,:,i) = d;
    
  end
  
  % Cap all inter-point distances at c^p
  if truncated
    % Truncated window
    distances = min(c^p,distances,'includenan');
  else
    % Full window
    distances = min(c^p,distances);
  end
  
  % Compute the OSPA(2) at each evaluation point
  for i = 1:num_eval
    
    % Window indices
    win_idx = eval_idx(i) + win_off;
    idx_val = (win_idx > 0) & (win_idx <= num_step);
    win_idx = win_idx(idx_val);
    
    % Compute the matrix of weighted time-averaged
    % OSPA distances between tracks
    trk_dist = mean(distances(:,:,win_idx),3,'omitnan');
    trk_dist(isnan(trk_dist)) = 0;
    
    % Get the number of objects in X and Y that exist
    % at any time inside the current window
    valid_rows = any(x_exists(:,win_idx),2);
    valid_cols = any(y_exists(:,win_idx),2);
    m = sum(valid_rows);
    n = sum(valid_cols);
    
    % Solve the optimal assignment problem
    trk_dist = trk_dist(valid_rows,valid_cols);
    if isempty(trk_dist)
        cost = 0;
    else
      if m > n
          trk_dist = trk_dist';
      end
      [~,cost] = lapjv(trk_dist);
    end
    
    % Compute the OSPA(2) distance for the current time index
    if max(m,n) == 0
      result(:,i) = 0;
    else
      result(1,i) = ( ( c^p * abs(m-n) + cost ) / max(m,n) ) ^ (1/p);
      result(2,i) = ( cost / max(m,n) ) ^ (1/p);
      result(3,i) = ( c^p * abs(m-n) / max(m,n) ) ^ (1/p);
    end
    
  end
  
end
