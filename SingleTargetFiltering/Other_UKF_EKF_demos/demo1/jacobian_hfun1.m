% Return Jacobian for hfun1

function Jy = jacobian_hfun1(x,u)
  
  R = u(1);

  Jy = [0 -x(3)/(R^2) x(2)/(R^2)];