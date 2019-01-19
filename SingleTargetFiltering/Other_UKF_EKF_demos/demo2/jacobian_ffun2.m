% Return Jacobian for ffun1

function Jx = jacobian_ffun1(x,u)
  
Jx = 1 + sin(x) + x*cos(x);