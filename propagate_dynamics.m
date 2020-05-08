function xkp1 = propagate_dynamics(sys,u,xk)

persistent w_mag
if isempty(w_mag)
    w_mag = 0.005;
end

xkp1 = xk + (sys.Ax*xk + sys.Bx*u)*sys.h + sys.Cx*w_mag*normrnd(0,1,4,1);%randn(2,1); 

% xddot(i) = (1/mx)*ux(i) - (bx/mx)*xdot(i) - (cx/mx)*x(i);
% xdot(i) = xdot(i-1) + xddot(i)*h;
% x(i) = x(i-1) + xdot(i)*h;
% yddot(i) = (1/my)*uy(i) - (by/my)*ydot(i) - (cy/my)*y(i);
% ydot(i) = ydot(i-1) + yddot(i)*h;
% y(i) = y(i-1) + ydot(i)*h;