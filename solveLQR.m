function solveLQR(sys,iter,fix_tf)
% solve LQR for different conditions
rng(17);
% % Time Span
tf = sys.tf;
h = sys.h;%0.01; % step size
n = tf/h;
tspan = linspace(0,tf,n);

% Desired state
xdes = 0;
xdesdot = 0;
ydes = 0;
ydesdot = 0;

%% LQR + uncertain parameters
xddot = zeros(n,1); yddot = zeros(n,1);
xdot = zeros(n,1); ydot = zeros(n,1);
x = zeros(n,1); y = zeros(n,1); 
ux = zeros(n,1); uy=zeros(n,1);

x(1) = sys.IC(1);
xdot(1) = sys.IC(2);
y(1) = sys.IC(3);
ydot(1) = sys.IC(4);

%real system params
mx = 10; bx=2; cx = 1;
my = 7; by = 1; cy = 0.6;

w_mag = sys.w_mag;

%randomized system parameters
mx_r = sys.mx_unc;%normrnd(mx, mx/10);
bx_r = sys.bx_unc;%normrnd(bx, bx/10);
cx_r = sys.cx_unc;%normrnd(cx, cx/10);
my_r = sys.my_unc;%normrnd(my, my/10);
by_r = sys.by_unc;%normrnd(by, by/10);
cy_r = sys.cy_unc;%normrnd(cy, cy/10);

Ax = [0 1; -cx_r/mx_r, -bx_r/mx_r];
Bx = [0 1/mx_r]';
Qx = [1 0; 0 1]; % changed from 100 to 10 then to 1
Rx = 0.1;
Kx = lqr(Ax,Bx,Qx,Rx,[]);

Ay = [0 1; -cy_r/my_r, -by_r/my_r];
By = [0 1/my_r]';
Qy = [1 0; 0 1]; % changed from 100 to 10 then to 1
Ry = 0.1;
Ky = lqr(Ay,By,Qy,Ry,[]);

Q = [1 0 0 0; 0 100 0 0; 0 0 1 0; 0 0 0 100];
R = [0.1 0; 0 0.1];

tol = sys.tol; %1e-3;

sys.Ax(2,1:2) = -sys.Ax(2,1:2);
sys.Ax(4,3:4) = -sys.Ax(4,3:4);
sys.Cx = sys.Cx;

i = 1;
while (i<n)
    i=i+1;
    % Control effort
    ux(i) = -Kx*[x(i-1); xdot(i-1)];
    uy(i) = -Ky*[y(i-1); ydot(i-1)];
    % Update dynamics
    u = [ux(i); uy(i)];
    xk = [x(i-1);xdot(i-1);y(i-1);ydot(i-1)];
    [xkp1] = propagate_dynamics(sys,u,xk);
    x(i) = xkp1(1);
    xdot(i) = xkp1(2);
    y(i) = xkp1(3);
    ydot(i) = xkp1(4);
%     xddot(i) = (1/mx)*ux(i) - (bx/mx)*xdot(i-1) - (cx/mx)*x(i-1) - w_mag*normrnd(0,1);
%     xdot(i) = xdot(i-1) + xddot(i)*h;
%     x(i) = x(i-1) + xdot(i)*h;
%     yddot(i) = (1/my)*uy(i) - (by/my)*ydot(i-1) - (cy/my)*y(i-1) - w_mag*normrnd(0,1);
%     ydot(i) = ydot(i-1) + yddot(i)*h;
%     y(i) = y(i-1) + ydot(i)*h;
    
    % Update errors
    x_error = x(i)-xdes;
    y_error = y(i) - ydes;
    if ~fix_tf
        if (abs(x_error)<tol  && abs(y_error)<tol)
            break
        end
    end
end

x_r = x(1:i);
y_r = y(1:i);
xdot_r = xdot(1:i);
ydot_r = ydot(1:i);
t = tspan(1:i);
ux_r = ux(1:i);
uy_r = uy(1:i);


%% Plot
mid = 5;

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

baseFileName = sprintf('LQR#%d.png', iter);
fullFileName = fullfile('Plots', baseFileName);
saveas(gcf, fullFileName);

close(gcf)

%% solve for cost_function
cost_function = 0;
for j=1:length(x_r)
    xVec = [x_r(j); xdot_r(j); y_r(j); ydot_r(j)];
    xTerm = xVec'*Q*xVec;
    uVec = [ux_r(j); uy_r(j)];
    uTerm = uVec'*R*uVec;
    cost_function = cost_function + xTerm + uTerm;
end

%% save cost_function and time of completion to file
fid = fopen('data.txt','a');
fprintf(fid,'LQR Data for Iter #%d with w_mag=%0.4f\n', iter,sys.w_mag);
fprintf(fid,'Time: %f [sec] \n', t(end));
fprintf(fid,'Cost function: %f \n\n',cost_function);
fclose(fid);


