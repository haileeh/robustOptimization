function xkp1 = propagate_dynamics(sys,u,xk)

% persistent w_mag
% if isempty(w_mag)
%     w_mag = 0.005;
% end

xkp1 = xk + (sys.Ax*xk + sys.Bx*u)*sys.h + sys.Cx*sys.w_mag*normrnd(0,1,4,1);