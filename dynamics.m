clear all
% Time Span
tf = 100;
h = 0.01; % step size
n = tf/h;
tspan = linspace(0,tf,n);

xddot = zeros(n,1); yddot = zeros(n,1); thetaddot = zeros(n,1);
xdot = zeros(n,1); ydot = zeros(n,1); thetadot = zeros(n,1);
x = zeros(n,1); y = zeros(n,1); theta = zeros(n,1);
ux = zeros(n,1); uy=zeros(n,1); utheta = zeros(n,1);
x(1) = 1;
y(1) = -1.3;
theta(1) = pi/2;
thetadot(1) = pi/8;
mx = 10; bx=2; cx = 1;
my = 7; by = 1; cy = 0.6;
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

Atheta = [0 1; 0 0];
Btheta = [0 1/I]';
Qtheta = [100 0; 0 1];
Rtheta = 10;
Ktheta = dlqr(Atheta,Btheta,Qtheta,Rtheta,[]);

% Desired state
xdes = 0;
xdesdot = 0;
ydes = 0;
ydesdot = 0;

i = 1;
while (i<n)
    i=i+1;
    % Control effort
    ux(i) = -Kx*[x(i-1); xdot(i-1)];
    uy(i) = -Ky*[y(i-1); ydot(i-1)];
%     utheta(i) = -Ktheta*[theta(i-1); thetadot(i-1)];
    % Update dynamics
    xddot(i) = (1/mx)*ux(i) - (bx/mx)*xdot(i) - (cx/mx)*x(i);
    xdot(i) = xdot(i-1) + xddot(i)*h; 
    x(i) = x(i-1) + xdot(i)*h; 
    yddot(i) = (1/my)*uy(i) - (by/my)*ydot(i) - (cy/my)*y(i);
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
    if (abs(x_error)<1e-3 && abs(x_error_dot)<1e-3 && abs(y_error)<1e-3 && abs(y_error_dot)<1e-3)
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

figure; 
subplot(2,1,1)
plot(t,x_nn); 
hold on; grid on;
ylabel('Distance');
yyaxis right
plot(t,xdot_nn);
xlabel('Time [s]');
ylabel('Velocity');
title('X-Position and X-Velocity');
subplot(2,1,2)
plot(t,y_nn); grid on;
ylabel('Distance');
yyaxis right
plot(t,ydot_nn);
xlabel('Time [s]');
ylabel('Velocity');
title('Y-Position and Y-Velocity');
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

x(1) = 1;
y(1) = -1.3;
theta(1) = pi/2;
thetadot(1) = pi/8;

rng(17);
w = wgn(n,1,10);
v = wgn(n,1,10);

i = 1;
while (i<n)
    i=i+1;
    % Control effort
    ux(i) = -Kx*[x(i-1); xdot(i-1)];
    uy(i) = -Ky*[y(i-1); ydot(i-1)];
    % Update dynamics
    xddot(i) = (1/mx)*ux(i) - (bx/mx)*xdot(i) - (cx/mx)*x(i) - w(i);
    %xddot(i) = awgn(xddot(i),20);
    xdot(i) = xdot(i-1) + xddot(i)*h; 
    x(i) = x(i-1) + xdot(i)*h; 
    yddot(i) = (1/my)*uy(i) - (by/my)*ydot(i) - (cy/my)*y(i) - v(i);
    %yddot(i) = awgn(yddot(i),20);
    ydot(i) = ydot(i-1) + yddot(i)*h;
    y(i) = y(i-1) + ydot(i)*h;
   
    % Update errors
    x_error = x(i)-xdes;
    x_error_dot = xdot(i)-xdesdot;
    y_error = y(i) - ydes;
    y_error_dot = ydot(i) - ydesdot;
    if (abs(x_error)<1e-3 && abs(x_error_dot)<1e-3 && abs(y_error)<1e-3 && abs(y_error_dot)<1e-3)
        break
    end
end

x_n = x(1:i);
y_n = y(1:i);
xdot_n = xdot(1:i);
ydot_n = ydot(1:i);
t = tspan(1:i);

figure; 
subplot(2,1,1)
plot(t,x_n); 
hold on; grid on;
ylabel('Distance');
yyaxis right
plot(t,xdot_n);
xlabel('Time [s]');
ylabel('Velocity');
title('X-Position and X-Velocity');
subplot(2,1,2)
plot(t,y_n); grid on;
ylabel('Distance');
yyaxis right
plot(t,ydot_n);
xlabel('Time [s]');
ylabel('Velocity');
title('Y-Position and Y-Velocity');

%% randomize system parameters
xddot = zeros(n,1); yddot = zeros(n,1); thetaddot = zeros(n,1);
xdot = zeros(n,1); ydot = zeros(n,1); thetadot = zeros(n,1);
x = zeros(n,1); y = zeros(n,1); theta = zeros(n,1);
ux = zeros(n,1); uy=zeros(n,1); utheta = zeros(n,1);

x(1) = 1;
y(1) = -1.3;
theta(1) = pi/2;
thetadot(1) = pi/8;

%real system params
mx = 10; bx=2; cx = 1;
my = 7; by = 1; cy = 0.6;
I = 0.8;
%randomized system parameters
mx_r = normrnd(mx, mx/10);
bx_r = normrnd(bx, bx/10);
cx_r = normrnd(cx, cx/10);
my_r = normrnd(my, my/10);
by_r = normrnd(by, by/10);
cy_r = normrnd(cy, cy/10);
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

i = 1;
while (i<n) 
    i=i+1;
    % Control effort
    ux(i) = -Kx*[x(i-1); xdot(i-1)];
    uy(i) = -Ky*[y(i-1); ydot(i-1)];
    % Update dynamics
    xddot(i) = (1/mx)*ux(i) - (bx/mx)*xdot(i) - (cx/mx)*x(i);
    xdot(i) = xdot(i-1) + xddot(i)*h; 
    x(i) = x(i-1) + xdot(i)*h; 
    yddot(i) = (1/my)*uy(i) - (by/my)*ydot(i) - (cy/my)*y(i);
    ydot(i) = ydot(i-1) + yddot(i)*h;
    y(i) = y(i-1) + ydot(i)*h;
   
    % Update errors
    x_error = x(i)-xdes;
    x_error_dot = xdot(i)-xdesdot;
    y_error = y(i) - ydes;
    y_error_dot = ydot(i) - ydesdot;
    if (abs(x_error)<1e-3 && abs(x_error_dot)<1e-3 && abs(y_error)<1e-3 && abs(y_error_dot)<1e-3)
        break
    end
end

x_r = x(1:i);
y_r = y(1:i);
xdot_r = xdot(1:i);
ydot_r = ydot(1:i);
t = tspan(1:i);

figure; 
subplot(2,1,1)
plot(t,x_r); 
hold on; grid on;
ylabel('Distance');
yyaxis right
plot(t,xdot_r);
xlabel('Time [s]');
ylabel('Velocity');
title('X-Position and X-Velocity');
subplot(2,1,2)
plot(t,y_r); grid on;
ylabel('Distance');
yyaxis right
plot(t,ydot_r);
xlabel('Time [s]');
ylabel('Velocity');
title('Y-Position and Y-Velocity');