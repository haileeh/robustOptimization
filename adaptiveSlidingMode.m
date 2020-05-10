function adaptiveSlidingMode(sys,iter)
% adaptive sliding mode

% % Time Span
tf = sys.tf;%100;
h = 0.01; % step size
n = tf/h;
tspan = linspace(0,tf,n);

% n = sys.tf/sys.h;
% h = sys.h;
% tspan = sys.tspan;

xddot = zeros(n,1); yddot = zeros(n,1);
xdot = zeros(n,1); ydot = zeros(n,1);
x = zeros(n,1); y = zeros(n,1);
ux = zeros(n,1); uy=zeros(n,1);
x(1) = sys.IC(1);
xdot(1) = sys.IC(2);
y(1) = sys.IC(3);
ydot(1) = sys.IC(4);

mx = 10; bx=2; cx = 1;
my = 7; by = 1; cy = 0.6;

%using uncertain parameters: dont need them here: maybe to initialize
%ahat?? no it made it worse
% mx = sys.mx_unc; bx = sys.bx_unc; cx = sys.cx_unc;
% my = sys.my_unc; by = sys.by_unc; cy = sys.cy_unc;

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

tol = sys.tol; 

rng(17);

% uncertainty
w_mag = sys.w_mag;

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
    xddot(i) = (1/mx)*ux(i) - (bx/mx)*xdot(i-1) - (cx/mx)*x(i-1) - w_mag*normrnd(0,1);
    xdot(i) = xdot(i-1) + xddot(i)*h; 
    x(i) = x(i-1) + xdot(i)*h; 
    yddot(i) = (1/my)*uy(i) - (by/my)*ydot(i-1) - (cy/my)*y(i-1) - w_mag*normrnd(0,1);
    ydot(i) = ydot(i-1) + yddot(i)*h;
    y(i) = y(i-1) + ydot(i)*h;

    % Update errors
    x_error = x(i)-xdes;
    x_error_dot = xdot(i)-xdesdot;
    y_error = y(i)-ydes;
    y_error_dot = ydot(i)-ydesdot;

    if (abs(x_error)<tol && abs(y_error)<tol)
        break
    end
end

x_nn = x(1:i);
y_nn = y(1:i);

xdot_nn = xdot(1:i);
ydot_nn = ydot(1:i);

t = tspan(1:i);
ux = ux(1:i);
uy = uy(1:i);

sx = sx(1:i);
sy = sy(1:i);

mid = 20;

figure;
plot(t(1:mid:end),x_nn(1:mid:end),'-*');
hold on; grid on;
plot(t(1:mid:end),y_nn(1:mid:end),'-.','LineWidth',2);
plot(t,zeros(length(t),1),'k');
title('ASM: Time vs. Distance');
xlabel('Time [s]');
ylabel('Distance');
legend('x','y');
axis tight

baseFileName = sprintf('ASM#%d.png', iter);
fullFileName = fullfile('Plots', baseFileName);
saveas(gcf, fullFileName);
close(gcf)

if 0
    figure;
    plot(t(1:mid:end),ux(1:mid:end),'-*');
    hold on; grid on;
    plot(t(1:mid:end),uy(1:mid:end),'-.','LineWidth',2);
    plot(t,zeros(length(t),1),'k');
    legend('u_x','u_y');
    xlabel('Time [s]');
    ylabel('u');
    title('Control Effort');
    axis tight
    
    baseFileName = sprintf('ASM_ctrl#%d.png', iter);
    fullFileName = fullfile('Plots', baseFileName);
    saveas(gcf, fullFileName);
    close(gcf)
    
    figure;
    plot(t(1:mid:end),sx(1:mid:end),'-*');
    hold on; grid on;
    plot(t(1:mid:end),sy(1:mid:end),'-.','LineWidth',2);
    plot(t,zeros(length(t),1),'k');
    legend('s_x','s_y');
    xlabel('Time [s]');
    ylabel('s');
    title('Sliding Variable');
    axis tight
    
    baseFileName = sprintf('ASM_sliding#%d.png', iter);
    fullFileName = fullfile('Plots', baseFileName);
    saveas(gcf, fullFileName);
    close(gcf)
end

% wasn't formulated with LQR so makes sense the cost function would be high
Q = [1 0 0 0; 0 100 0 0; 0 0 1 0; 0 0 0 100];
R = [0.1 0; 0 0.1];

%% solve for cost_function
cost_function = 0;
for j=1:length(x_nn)
    xVec = [x_nn(j); xdot_nn(j); y_nn(j); ydot_nn(j)];
    xTerm = xVec'*Q*xVec;
    uVec = [ux(j); uy(j)];
    uTerm = uVec'*R*uVec;
    cost_function = cost_function + xTerm + uTerm;
end

%% save cost_function and time of completion to file
fid = fopen('data.txt','a');
fprintf(fid,'ASM Data for Iter #%d with w_mag=%0.4f\n', iter,sys.w_mag);
fprintf(fid,'Time: %f [sec] \n', t(end));
fprintf(fid,'Cost function: %f \n\n',cost_function);
fclose(fid);
