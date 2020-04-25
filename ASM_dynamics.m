% adaptive sliding mode
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

% new vectors
sx = zeros(n,1);
ax_hat_dot = zeros(2,n);
Kx_hat_dot = zeros(n,1);
ax_hat = zeros(2,n);
Kx_hat = zeros(n,1);
sy = zeros(n,1);
ay_hat_dot = zeros(2,n);
Ky_hat_dot = zeros(n,1);
ay_hat = zeros(2,n);
Ky_hat = zeros(n,1);

% parameters
lambdax = 1;
gammax = 4;
Gammax = 5;
lambday = 0.5;
gammay = 3;
Gammay = 4;

% initialize error terms
xdes = 0; ydes = 0;
xdesdot = 0; ydesdot = 0;
x_error = x(1)-xdes;
x_error_dot = xdot(1)-xdesdot;
y_error = y(1)-ydes;
y_error_dot = ydot(1)-ydesdot;

rng(15);
w = wgn(n,1,0.01);
v = wgn(n,1,0.01);

i=1;
while ( i<n ) 
    i=i+1;
    % compute s
    sx(i) = x_error_dot + lambdax*x_error;
    sy(i) = y_error_dot + lambday*y_error;
    % Construct Y
    Yx = [xdot(i-1), x(i-1)];
    Yy = [ydot(i-1), y(i-1)];
    % adaptation laws
    ax_hat_dot(:,i) = -Gammax*Yx'*sx(i);
    ay_hat_dot(:,i) = -Gammay*Yy'*sy(i);
    Kx_hat_dot(i) = -gammax*sx(i)^2;
    Ky_hat_dot(i) = -gammay*sy(i)^2;
    % integrate the adaptation laws
    ax_hat(:,i) = ax_hat(:,i-1) + ax_hat_dot(:,i)*h;
    ay_hat(:,i) = ay_hat(:,i-1) + ay_hat_dot(:,i)*h;
    Kx_hat(i) = Kx_hat(i-1) + Kx_hat_dot(i)*h;
    Ky_hat(i) = Ky_hat(i-1) + Ky_hat_dot(i)*h;
    % solve for u_hat
    xrddot = -lambdax*x_error_dot;
    yrddot = -lambday*y_error_dot;
    ux_hat = Yx*ax_hat(:,i) + xrddot;
    uy_hat = Yy*ay_hat(:,i) + yrddot;
    % Solve for u_p
    ux_p = -Kx_hat(i)*sx(i);
    uy_p = -Ky_hat(i)*sy(i);
    % control effort
    ux(i) = ux_hat - ux_p;
    uy(i) = uy_hat - uy_p;

    % Update dynamics
    xddot(i) = (1/mx)*ux(i) - (bx/mx)*xdot(i) - (cx/mx)*x(i);% -w(i);
    xddot(i) = awgn(xddot(i),20);
    xdot(i) = xdot(i-1) + xddot(i)*h; 
    x(i) = x(i-1) + xdot(i)*h; 
    yddot(i) = (1/my)*uy(i) - (by/my)*ydot(i) - (cy/my)*y(i);%-v(i);
    yddot(i) = awgn(yddot(i),20);
    ydot(i) = ydot(i-1) + yddot(i)*h;
    y(i) = y(i-1) + ydot(i)*h;

    % Update errors
    x_error = x(i)-xdes;
    x_error_dot = xdot(i)-xdesdot;
    y_error = y(i)-ydes;
    y_error_dot = ydot(i)-ydesdot;

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
sx = sx(1:i);
sy = sy(1:i);

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

figure; 
plot(t,ux,t,uy);
hold on; grid on;
legend('u_x','u_y');
xlabel('Time [s]');
ylabel('u');
title('Control Effort');

figure;
plot(t,sx);
hold on; grid on;
plot(t,sy);
legend('s_x','s_y');
xlabel('Time [s]');
ylabel('s');
title('Sliding Variable');