function [assignments,costs]= murty(P0,m)
% MURTY   Murty's Algorithm for m-best ranked optimal assignment problem
% Ba Tuong Vo 2015

  % Find the optimal and initial solution
  [S0, C0] = assignmentoptimal(P0);
  S0 = S0';
  num_rows = size(P0,1);
  num_cols = size(P0,2);

  if m == 1
    assignments = S0;
    costs = C0;
    return; 
  end

  % Preallocate a block of memory to hold the queue
  N = 1000;  % Block size
  answer_list_P = zeros(num_rows,num_cols,N);
  answer_list_S = zeros(num_rows,N);
  answer_list_C = NaN*ones(1,N);

  % Initialize answer list
  answer_list_P(:,:,1) = P0;   % problem or cost matrix
  answer_list_S(:,1) = S0';    % solutions or assignemnts
  answer_list_C(1) = C0;       % cost vector for problems/solutions
  answer_index_next = 2;

  assignments = zeros(num_rows,m);
  costs = zeros(m,1);

  for i = 1:m
    
    % If all are cleared, break the loop early 
    if all(isnan(answer_list_C))
        assignments= assignments(:,1:answer_index_next-1);
        costs= costs(1:answer_index_next-1);
        break;
    end
    
    % Grab lowest cost solution index
    [notused,idx_top] = min(answer_list_C(1:answer_index_next-1));
   
    % Copy the current best solution out
    assignments(:,i) = answer_list_S(:,idx_top);
    costs(i,1) = answer_list_C(idx_top);

    % Copy lowest cost problem to temp
    P_now = answer_list_P(:,:,idx_top);
    S_now = answer_list_S(:,idx_top);
   
    % Delete the solution from the queue
    answer_list_C(idx_top) = NaN;
   
    for a = 1:length(S_now)
       
      % Current assignment pair
      aw = a;
      aj = S_now(a);
       
      if aj ~= 0

          % Remove it and calculate new solution
          P_tmp = P_now;
          if aj<=num_cols-num_rows
              P_tmp(aw,aj) = inf;
          else
              P_tmp(aw,num_cols-num_rows+1:end) = inf;
          end
        
        [S_tmp,C_tmp] = assignmentoptimal(P_tmp);
        S_tmp = S_tmp';
           
        % Copy to new list
        if all(S_tmp ~= 0)
          
          % If we have filled the allocated space, allocate more
          if answer_index_next > length(answer_list_C)
            answer_list_P = cat(3,answer_list_P,zeros(num_rows,num_cols,N));
            answer_list_S = cat(2,answer_list_S,zeros(num_rows,N));
            answer_list_C = cat(2,answer_list_C,NaN*ones(1,N));
          end
          
          answer_list_P(:,:,answer_index_next) = P_tmp;
          answer_list_S(:,answer_index_next) = S_tmp;
          answer_list_C(answer_index_next) = C_tmp;
          answer_index_next = answer_index_next + 1;
          
        end
           
        % Enforce current assignment
        v_tmp = P_now(aw,aj);
        P_now(aw,:) = inf;
        P_now(:,aj )= inf;
        P_now(aw,aj) = v_tmp;
           
      end
    end
   
  end
    
  assignments = assignments';
  costs = costs';

end