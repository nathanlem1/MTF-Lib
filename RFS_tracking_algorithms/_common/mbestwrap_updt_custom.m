function [assignments,costs]= murty_wrapper_update(P0,m)

if m==0
    assignments= [];
    costs= [];
    return;
end
    
  n1 = size(P0,1);
  n2 = size(P0,2);

  % Padding blocks for dummy variables
  blk1 = -log(ones(n1,n1));

  P0 = [P0 blk1];

  % Make costs non-negative (required by 'assignmentoptimal')
  x = min(min(P0));
  P0 = P0 - x;
  
  % Murty
  [assignments, costs] = murty_custom(P0,m);

  % Restore correct costs to assignments
  costs = costs + (x.*sum(assignments>0,2))';
      
  % Strip dummy variables
  assignments = assignments(:,1:n1);

  % Dummy assignments are clutter
  assignments(assignments>n2) = 0;

      
end
    