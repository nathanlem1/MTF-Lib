function [assignments,costs]= murty_wrapper_update(P0,m)

if m==0
    assignments= [];
    costs= [];
    return;
end


  % Make costs non-negative (required by 'assignmentoptimal')
  x = min(min(P0));
  P0 = P0 - x;
  
  % Murty
  [assignments, costs] = murty_custom(P0,m);

  % Restore correct costs to assignments
  costs = costs + (x.*sum(assignments>0,2))';

  
end
    