function [xkp1] = propagate_dynamics(sys,u,xk)
% Euler integration
xkp1 = xk + (sys.Ax*xk + sys.Bx*u)*sys.h + sys.Cx*sys.w_mag*normrnd(0,1,4,1);