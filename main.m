%% Set up true parameters
% Time Span
tf = 15;%10;
sys.h = 0.05;%0.1; % step size
n = tf/sys.h;
tspan = linspace(0,tf,n);
sys.tf = tf;
sys.tspan = tspan;

% Initial conditions
x0 = 1;
xdot0 = 0;
y0 = -1.3;
ydot0 = 0;

% certain system parameters
mx = 10; bx=2; cx = 1;
my = 7; by = 1; cy = 0.6;

sys.mx = mx; sys.bx = bx; sys.cx = cx;
sys.my = my; sys.by = by; sys.cy = cy;

rng(17);

%uncertain system parameters
varArray = [10; 25; 50];
% mx_unc = normrnd(mx,mx/10);
% bx_unc = normrnd(bx,bx/10);
% cx_unc = normrnd(cx,cx/10);
% my_unc = normrnd(my,my/10);
% by_unc = normrnd(by,by/10);
% cy_unc = normrnd(cy,cy/10);
%
% sys.mx_unc = mx_unc; sys.bx_unc = bx_unc; sys.cx_unc = cx_unc;
% sys.my_unc = my_unc; sys.by_unc = by_unc; sys.cy_unc = cy_unc;

sys.tol = 1e-2;

% sys.w_mag = 0.005;
wArray = [0.005; 0.05; 0.5; 1];
sys.gamma = 10; % 1 or 10 only for adaptive; sdp can do full range

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

% for loop should go here, vary iter, gamma, w_mag
iter = 1;
for varIter =1:3
    varDen = varArray(varIter);
    mx_unc = normrnd(mx,mx/varDen);
    bx_unc = normrnd(bx,bx/varDen);
    cx_unc = normrnd(cx,cx/varDen);
    my_unc = normrnd(my,my/varDen);
    by_unc = normrnd(by,by/varDen);
    cy_unc = normrnd(cy,cy/varDen);
    
    sys.mx_unc = mx_unc; sys.bx_unc = bx_unc; sys.cx_unc = cx_unc;
    sys.my_unc = my_unc; sys.by_unc = by_unc; sys.cy_unc = cy_unc;
    for counter=1:4
        sys.w_mag = wArray(counter);
        sys.IC = [x0; xdot0;y0;ydot0];
        xk = zeros(4,n);
        xk(:,1) = [x0,xdot0,y0,ydot0];
        
        %% LQR
        solveLQR(sys,iter)
        
        %% Adaptive Sliding Mode
        adaptiveSlidingMode(sys,iter)
        
        %% SDP
        N = 4;
        c=1;
        tol = 1e-2;
        uspan = zeros(8,n);
        for j=0:sys.h:tf
            u = robust_sdp(sys,N);
            [xk(:,c+1)] = propagate_dynamics(sys,u(1:2),xk(:,c));
            sys.IC = xk(:,c+1);
            uspan(:,c) = u;
            c=c+1;
            if (abs(xk(1,c))<tol && abs(xk(3,c))<tol)
                break
            end
        end
        
        % Plotting
        if c < length(tspan)
            t=tspan(1:c);
            xk = xk(:,1:c);
            uspan = uspan(:,1:c);
        else
            t = tspan;
            xk = xk(:,2:end-1);
        end
        comp_xk = xk;
        comp_t = t;
        
        %% Adaptive SDP
        N=4;
        c=1;
        uprev=zeros(8,1);
        xk = zeros(4,n);
        u_asdp = zeros(8,n);
        sys.IC = [x0; xdot0;y0;ydot0];
        xk(:,1) = [x0,xdot0,y0,ydot0];
        for j=0:sys.h:tf
            u = adaptive_robust_sdp(sys,N,c-1,uprev);
            [xk(:,c+1)] = propagate_dynamics(sys,u(1:2),xk(:,c));
            sys.IC = xk(:,c+1);
            u_asdp(:,c) = u;
            c=c+1;
            uprev = [u(1:2);u(1:2);u(1:2);u(1:2)];
            if (abs(xk(1,c))<tol && abs(xk(3,c))<tol)
                break
            end
        end
        
        if c <= length(tspan)
            t=tspan(1:c);
            xk = xk(:,1:c);
            u_asdp = u_asdp(:,1:c);
        else
            t = tspan;
            xk = xk(:,2:end-1);
        end
        
        %% Plotting
        figure;
        plot(comp_t,comp_xk(1,:),'-*');
        hold on; grid on;
        plot(comp_t,comp_xk(3,:),'-.','LineWidth',2);
        plot(comp_t,zeros(length(comp_t),1),'k');
        title('SDP: Time vs. Distance');
        xlabel('Time [s]');
        ylabel('Distance');
        legend('x','y');
        axis tight
        
        baseFileName = sprintf('SDP#%d.png', iter);
        fullFileName = fullfile('Plots', baseFileName);
        saveas(gcf, fullFileName);
        close(gcf)
        
        figure;
        plot(t,xk(1,:),'-*');
        hold on; grid on;
        plot(t,xk(3,:),'-.','LineWidth',2);
        plot(t,zeros(length(t),1),'k');
        title('A-SDP: Time vs. Distance');
        xlabel('Time [s]');
        ylabel('Distance');
        legend('x','y');
        axis tight
        
        baseFileName = sprintf('ASDP#%d.png', iter);
        fullFileName = fullfile('Plots', baseFileName);
        saveas(gcf, fullFileName);
        close(gcf)
        
        n1 = length(comp_xk);
        n2 = length(xk);
        if n1>n2
            len = n2;
        else
            len = n1;
        end
        
        figure;
        plot(tspan(1:len),xk(1,1:len)-comp_xk(1,1:len),'-*');
        hold on; grid on;
        plot(tspan(1:len),xk(3,1:len)-comp_xk(3,1:len),'-.','LineWidth',2);
        plot(tspan(1:len),zeros(length(tspan(1:len)),1),'k');
        title('Difference in Positions');
        legend('X','Y');
        xlabel('Time [s]');
        ylabel('A-SDP - SDP: Distance Difference');
        axis tight
        
        %% solve for cost_function for SDP
%         Q = [1 0 0 0; 0 100 0 0; 0 0 1 0; 0 0 0 100];
        Q = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
        R = [0.1 0; 0 0.1];
        
        cost_function_sdp = zeros(length(comp_xk),1);
        for j=1:length(comp_xk)
            xVec = comp_xk(:,j);
            xTerm = xVec'*Q*xVec;
            uVec = uspan(1:2,j);
            uTerm = uVec'*R*uVec;
            cost_function_sdp(j) = xTerm + uTerm;
        end
        
        %% solve for cost_function for A-SDP        
        cost_function_asdp = zeros(length(xk),1);
        for j=1:length(xk)
            xVec = xk(:,j);
            xTerm = xVec'*Q*xVec;
            uVec = u_asdp(1:2,j);
            uTerm = uVec'*R*uVec;
            cost_function_asdp(j) = xTerm + uTerm;
        end
        
        %% save cost_function and time of completion to file
        fid = fopen('data.txt','a');
        fprintf(fid,'SDP Data for Iter #%d with w_mag=%0.4f\n', iter,sys.w_mag);
        fprintf(fid,'Time: %f [sec] \n', comp_t(end));
        fprintf(fid,'Cost function: %f \n\n',sum(cost_function_sdp));%cost_function_sdp(1));
        
        fprintf(fid,'A-SDP Data for Iter #%d with w_mag=%0.4f\n', iter,sys.w_mag);
        fprintf(fid,'Time: %f [sec] \n', t(end));
        fprintf(fid,'Cost function: %f \n\n',sum(cost_function_asdp));%cost_function_asdp(1));
        fclose(fid);
        
        iter = iter+1;
    end
end