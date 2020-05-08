%% Set up true parameters
% Time Span
tf = 10;
sys.h = 0.1; % step size
n = tf/sys.h;
tspan = linspace(0,tf,n);
sys.tf = tf;
sys.tspan = tspan;

xk = zeros(4,n);
% Initial conditions
x0 = 1;
xdot0 = 0;
y0 = -1.3;
ydot0 = 0;
sys.IC = [x0; xdot0;y0;ydot0];
xk(:,1) = [x0,xdot0,y0,ydot0];
% theta0 = pi/2;
% thetadot0 = pi/8;
mx = 10; bx=2; cx = 1;
my = 7; by = 1; cy = 0.6;

sys.mx = mx; sys.bx = bx; sys.cx = cx;
sys.my = my; sys.by = by; sys.cy = cy;


%uncertain system parameters
mx_unc = normrnd(mx,mx/10);
bx_unc = normrnd(bx,bx/10);
cx_unc = normrnd(cx,cx/10);
my_unc = normrnd(my,my/10);
by_unc = normrnd(by,by/10);
cy_unc = normrnd(cy,cy/10);

sys.mx_unc = mx_unc; sys.bx_unc = bx_unc; sys.cx_unc = cx_unc;
sys.my_unc = my_unc; sys.by_unc = by_unc; sys.cy_unc = cy_unc;

% I = 0.8;

% [xdot ] = [0 1 0 0][x   ] + [0]
% [xddot]   [c/m b/m 0 0][xdot] + [1/m]u
% [ydot ] = [0 0 0 1][y   ] + [0]
% [yddot] = [0 0 c/m b/m][ydot] + [1/m]u

% y = B xvec
sys.Ax = [0 1 0 0; cx/mx, bx/mx 0 0; 0 0 0 1;0 0 cy/my, by/my];
sys.Bx = [0 0; 1/mx 0; 0 0; 0 1/my];%[0 1/mx 0 1/my]';
sys.Cx = eye(4);
sys.Qx = eye(4);
sys.Rx = 0.1;

sys.Ax_unc = [0 1 0 0; cx_unc/mx_unc, bx_unc/mx_unc 0 0; 0 0 0 1; 0 0 cy_unc/my_unc, by_unc/my_unc];
sys.Bx_unc = [0 0;1/mx_unc 0; 0 0;0 1/my_unc];%[0 1/mx_unc 0 1/my_unc]';
sys.Cx_unc = eye(4);

rng(17);
%% LQR
solveLQR(sys)
keyboard

%% Adaptive Sliding Mode
adaptiveSlidingMode(sys)
keyboard
%% SDP
N = 4;
c=1;
for j=0:sys.h:tf
   u = robust_sdp(sys,N);
   xk(:,c+1) = propagate_dynamics(sys,u(1:2),xk(:,c));
   sys.IC = xk(:,c+1);
   c=c+1;
end

%% Plotting
figure; 
plot(tspan,xk(1,2:end-1));
hold on; grid on;
plot(tspan,xk(3,2:end-1));
plot(tspan,zeros(length(tspan),1),'k');
title('SDP: Time vs. Distance');
xlabel('Time [s]');
ylabel('Distance');
legend('x','y');
