function cost_function = solveLQR(sys)
% solve LQR for different conditions

% % Time Span
tf = 50;
h = sys.h;%0.01; % step size
n = tf/h;
tspan = linspace(0,tf,n);

% n = sys.tf/sys.h;
% h = sys.h;
% tspan = sys.tspan;

xddot = zeros(n,1); yddot = zeros(n,1); thetaddot = zeros(n,1);
xdot = zeros(n,1); ydot = zeros(n,1); thetadot = zeros(n,1);
x = zeros(n,1); y = zeros(n,1); theta = zeros(n,1);
ux = zeros(n,1); uy=zeros(n,1); utheta = zeros(n,1);
x(1) = sys.IC(1);
xdot(1) = sys.IC(2);
y(1) = sys.IC(3);
ydot(1) = sys.IC(4);

theta(1) = pi/2;
thetadot(1) = pi/8;
% mx = 10; bx=2; cx = 1;
% my = 7; by = 1; cy = 0.6;

mx = sys.mx; bx = sys.bx; cx = sys.cx;
my = sys.my; by = sys.by; cy = sys.cy;

I = 0.8;
% [xdot ] = [0 1    ][x   ] + [0]
% [xddot]   [c/m b/m][xdot] + [1]u
% y = B xvec
Ax = [0 1; cx/mx, bx/mx];
Bx = [0 1/mx]';
Qx = [1 0; 0 100];
Rx = 0.1;
Kx = dlqr(Ax,Bx,Qx,Rx,[]);

Ay = [0 1; cy/my, by/my];
By = [0 1/my]';
Qy = [1 0; 0 100];
Ry = 0.1;
Ky = dlqr(Ay,By,Qy,Ry,[]);

% Atheta = [0 1; 0 0];
% Btheta = [0 1/I]';
% Qtheta = [100 0; 0 1];
% Rtheta = 10;
% Ktheta = dlqr(Atheta,Btheta,Qtheta,Rtheta,[]);

% Desired state
xdes = 0;
xdesdot = 0;
ydes = 0;
ydesdot = 0;

if 0
i = 1;
while (i<n)
    i=i+1;
    % Control effort
    ux(i) = -Kx*[x(i-1); xdot(i-1)];
    uy(i) = -Ky*[y(i-1); ydot(i-1)];
%     utheta(i) = -Ktheta*[theta(i-1); thetadot(i-1)];
    % Update dynamics
    xddot(i) = (1/mx)*ux(i) - (bx/mx)*xdot(i-1) - (cx/mx)*x(i-1);
    xdot(i) = xdot(i-1) + xddot(i)*h; 
    x(i) = x(i-1) + xdot(i)*h; 
    yddot(i) = (1/my)*uy(i) - (by/my)*ydot(i-1) - (cy/my)*y(i-1);
    ydot(i) = ydot(i-1) + yddot(i)*h;
    y(i) = y(i-1) + ydot(i)*h;
%     thetaddot(i) = (1/I)*utheta(i);
%     thetadot(i) = thetadot(i-1) + thetaddot(i)*h;
%     theta(i) = theta(i-1) + thetadot(i)*h;

    % Update errors
    x_error = x(i)-xdes;
    x_error_dot = xdot(i)-xdesdot;
    y_error = y(i) - ydes;
    y_error_dot = ydot(i) - ydesdot;
    if (abs(x_error)<1e-3 && abs(y_error)<1e-3)
        break
    end
end

x_nn = x(1:i);
y_nn = y(1:i);
theta_nn = theta(1:i);
xdot_nn = xdot(1:i);
ydot_nn = ydot(1:i);
thetadot_nn = thetadot(1:i);
t = tspan(1:i);
ux = ux(1:i);
uy = uy(1:i);
utheta = utheta(1:i);

mid = 10;

figure;
plot(t(1:mid:end),x_nn(1:mid:end),'-*');
hold on; grid on;
plot(t(1:mid:end),y_nn(1:mid:end),'-.','LineWidth',2);
plot(t,zeros(length(t),1),'k');
title('LQR wo Noise: Time vs. Distance');
xlabel('Time [s]');
ylabel('Distance');
legend('x','y');
axis tight

% figure; 
% subplot(2,1,1)
% plot(t,x_nn); 
% hold on; grid on;
% ylabel('Distance');
% yyaxis right
% plot(t,xdot_nn);
% xlabel('Time [s]');
% ylabel('Velocity');
% title('X-Position and X-Velocity');
% subplot(2,1,2)
% plot(t,y_nn); grid on;
% ylabel('Distance');
% yyaxis right
% plot(t,ydot_nn);
% xlabel('Time [s]');
% ylabel('Velocity');
% title('Y-Position and Y-Velocity');
% subplot(3,1,3)
% plot(t,theta_nn); grid on;
% ylabel('Distance');
% yyaxis right
% plot(t,thetadot_nn);
% xlabel('Time [s]');
% ylabel('Velocity');
% title('\theta-Position and \theta-Velocity');

% figure; 
% plot(tspan,ux);
% hold on; grid on;
% xlabel('Time [s]');
% ylabel('u');
% title('Control Effort');


% xddot = (1/m)*ux - (b/m)*xdot - (c/m)*x;
% yddot = (1/m)*uy - (b/m)*ydot - (c/m)*y;
% thetaddot = (1/I)*tau;

%% add noise
xddot = zeros(n,1); yddot = zeros(n,1); thetaddot = zeros(n,1);
xdot = zeros(n,1); ydot = zeros(n,1); thetadot = zeros(n,1);
x = zeros(n,1); y = zeros(n,1); theta = zeros(n,1);
ux = zeros(n,1); uy=zeros(n,1); utheta = zeros(n,1);
u = zeros(2,n);

x(1) = sys.IC(1);
xdot(1) = sys.IC(2);
y(1) = sys.IC(3);
ydot(1) = sys.IC(4);

theta(1) = pi/2;
thetadot(1) = pi/8;

Q = [1 0 0 0; 0 100 0 0; 0 0 1 0; 0 0 0 100];
Kcom = dlqr(sys.Ax,sys.Bx,Q,sys.Rx,[]);

rng(17);

w_mag = sys.w_mag;
% w = wgn(n,1,10);
% v = wgn(n,1,10);

i = 1; tol = 1e-3;
while (i<n)
    i=i+1;
    % Control effort
%     u(:,i) = -Kcom*[x(i-1); xdot(i-1); y(i-1); ydot(i-1)];
    ux(i) = -Kx*[x(i-1); xdot(i-1)];
    uy(i) = -Ky*[y(i-1); ydot(i-1)];
    % Update dynamics
    xddot(i) = (1/mx)*ux(i) - (bx/mx)*xdot(i-1) - (cx/mx)*x(i-1) - w_mag*normrnd(0,1);
    xdot(i) = xdot(i-1) + xddot(i)*h; 
    x(i) = x(i-1) + xdot(i)*h; 
    yddot(i) = (1/my)*uy(i) - (by/my)*ydot(i-1) - (cy/my)*y(i-1) - w_mag*normrnd(0,1);
    ydot(i) = ydot(i-1) + yddot(i)*h;
    y(i) = y(i-1) + ydot(i)*h;
%     u = [ux(i); uy(i)];
%     xk = [x(i-1);xdot(i-1);y(i-1);ydot(i-1)];
%     xkp1 = propagate_dynamics(sys,u,xk);
%     %spread out xkp1
%     x(i) = xkp1(1);
%     xdot(i) = xkp1(2);
%     y(i) = xkp1(3);
%     ydot(i) = xkp1(4);
   
    % Update errors
    x_error = x(i)-xdes;
    x_error_dot = xdot(i)-xdesdot;
    y_error = y(i) - ydes;
    y_error_dot = ydot(i) - ydesdot;

    if (abs(x_error)<tol && abs(y_error)<tol)
        break
    end
end

x_n = x(1:i);
y_n = y(1:i);
xdot_n = xdot(1:i);
ydot_n = ydot(1:i);
t = tspan(1:i);

mid = 10;

figure;
plot(t(1:mid:end),x_n(1:mid:end),'-*');
hold on; grid on;
plot(t(1:mid:end),y_n(1:mid:end),'-.','LineWidth',2);
plot(t,zeros(length(t),1),'k');
title('LQR w Noise: Time vs. Distance');
xlabel('Time [s]');
ylabel('Distance');
legend('x','y');
axis tight

% figure; 
% subplot(2,1,1)
% plot(t,x_n); 
% hold on; grid on;
% ylabel('Distance');
% yyaxis right
% plot(t,xdot_n);
% xlabel('Time [s]');
% ylabel('Velocity');
% title('X-Position and X-Velocity');
% subplot(2,1,2)
% plot(t,y_n); grid on;
% ylabel('Distance');
% yyaxis right
% plot(t,ydot_n);
% xlabel('Time [s]');
% ylabel('Velocity');
% title('Y-Position and Y-Velocity');

end
%% LQR + uncertain parameters
xddot = zeros(n,1); yddot = zeros(n,1); thetaddot = zeros(n,1);
xdot = zeros(n,1); ydot = zeros(n,1); thetadot = zeros(n,1);
x = zeros(n,1); y = zeros(n,1); theta = zeros(n,1);
ux = zeros(n,1); uy=zeros(n,1); utheta = zeros(n,1);

x(1) = sys.IC(1);
xdot(1) = sys.IC(2);
y(1) = sys.IC(3);
ydot(1) = sys.IC(4);

theta(1) = pi/2;
thetadot(1) = pi/8;

%real system params
mx = 10; bx=2; cx = 1;
my = 7; by = 1; cy = 0.6;
I = 0.8;

w_mag = sys.w_mag;

%randomized system parameters
mx_r = sys.mx_unc;%normrnd(mx, mx/10);
bx_r = sys.bx_unc;%normrnd(bx, bx/10);
cx_r = sys.cx_unc;%normrnd(cx, cx/10);
my_r = sys.my_unc;%normrnd(my, my/10);
by_r = sys.by_unc;%normrnd(by, by/10);
cy_r = sys.cy_unc;%normrnd(cy, cy/10);
% [xdot ] = [0 1    ][x   ] + [0]
% [xddot]   [c/m b/m][xdot] + [1]u
% y = B xvec
Ax = [0 1; cx_r/mx_r, bx_r/mx_r];
Bx = [0 1/mx_r]';
Qx = [1 0; 0 100];
Rx = 0.1;
Kx = dlqr(Ax,Bx,Qx,Rx,[]);

Ay = [0 1; cy_r/my_r, by_r/my_r];
By = [0 1/my_r]';
Qy = [1 0; 0 100];
Ry = 0.1;
Ky = dlqr(Ay,By,Qy,Ry,[]);

Q = [1 0 0 0; 0 100 0 0; 0 0 1 0; 0 0 0 100];
R = [0.1 0; 0 0.1];

i = 1;
while (i<n) 
    i=i+1;
    % Control effort
    ux(i) = -Kx*[x(i-1); xdot(i-1)];
    uy(i) = -Ky*[y(i-1); ydot(i-1)];
    % Update dynamics
    xddot(i) = (1/mx)*ux(i) - (bx/mx)*xdot(i-1) - (cx/mx)*x(i-1) - w_mag*normrnd(0,1);
    xdot(i) = xdot(i-1) + xddot(i)*h; 
    x(i) = x(i-1) + xdot(i)*h; 
    yddot(i) = (1/my)*uy(i) - (by/my)*ydot(i-1) - (cy/my)*y(i-1) - w_mag*normrnd(0,1);
    ydot(i) = ydot(i-1) + yddot(i)*h;
    y(i) = y(i-1) + ydot(i)*h;
   
    % Update errors
    x_error = x(i)-xdes;
    x_error_dot = xdot(i)-xdesdot;
    y_error = y(i) - ydes;
    y_error_dot = ydot(i) - ydesdot;
    if (abs(x_error)<1e-3  && abs(y_error)<1e-3)
        break
    end
end

x_r = x(1:i);
y_r = y(1:i);
xdot_r = xdot(1:i);
ydot_r = ydot(1:i);
t = tspan(1:i);
ux_r = ux(1:i);
uy_r = uy(1:i);

% figure; 
% subplot(2,1,1)
% plot(t,x_r); 
% hold on; grid on;
% ylabel('Distance');
% yyaxis right
% plot(t,xdot_r);
% xlabel('Time [s]');
% ylabel('Velocity');
% title('X-Position and X-Velocity');
% subplot(2,1,2)
% plot(t,y_r); grid on;
% ylabel('Distance');
% yyaxis right
% plot(t,ydot_r);
% xlabel('Time [s]');
% ylabel('Velocity');
% title('Y-Position and Y-Velocity');

mid = 10;

figure;
plot(t(1:mid:end),x_r(1:mid:end),'-*');
hold on; grid on;
plot(t(1:mid:end),y_r(1:mid:end),'-.','LineWidth',2);
plot(t,zeros(length(t),1),'k');
title('Uncertain Params, Noisy, LQR: Time vs. Distance');
xlabel('Time [s]');
ylabel('Distance');
legend('x','y');
axis tight

% solve for cost_function
cost_function = 0;
for j=1:length(x_r)
    xVec = [x_r(j); xdot_r(j); y_r(j); ydot_r(j)];
    xTerm = xVec'*Q*xVec;
    uVec = [ux_r(j); uy_r(j)];
    uTerm = uVec'*R*uVec;
    cost_function = cost_function + xTerm + uTerm;
end

cost_function