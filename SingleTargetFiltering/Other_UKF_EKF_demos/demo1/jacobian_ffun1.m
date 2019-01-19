% Return Jacobian for ffun1

function Jx = jacobian_ffun1(x,u)
  
  R = u(1);
  T = u(2);

  V = x(1);  
  
  theta = atan2(x(3),x(2));
  
  new_theta = theta + (V*T)/R;
  
  
  Jx = [ 1                  0                        0 ;
	-T*sin(new_theta)   x(3)*sin(new_theta)/R   -x(2)*sin(new_theta)/R ;
	 T*cos(new_theta)  -x(3)*cos(new_theta)/R    x(2)*cos(new_theta)/R];