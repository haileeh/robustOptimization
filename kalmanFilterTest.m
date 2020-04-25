% Problem 3
% Part 1
A = 1; B = 0; G = 1;
C = 1; D = 0;H = 0;

sys = ss(A,[B G], C, [D H],-1,'inputname',{'u' 'w'},'outputname','y');

Q = 1;
R = 1;
N = 0;
[kest, L, P] = kalman(sys, Q, R, N);

%% Part 2
% initiate
A_k = A; W_k = Q; C_k = C; V_k = R;
n = 11000;
v = wgn(n,1,1);
w = wgn(n,1,1);
x_hat(1) = 0;
x(1) = 0;
Q(1) = P;

% initial state propagation
x_plus(1) = A_k*x_hat(1);
Q_plus = A_k*Q(1)*A_k' + W_k;
y(1) = x(1) + v(1);
x(2) = A_k*x(1) + w(1);

res(1) = y(1) - C_k*x_plus(1);

for k = 2:n
    % regular system update
    x(k+1) = A_k*x(k) + w(k);
    y(k) = x(k) + v(k);
    
    % measurement update
    L(k) = Q_plus(k-1)*C_k'*inv(C_k*Q_plus(k-1)*C_k' + V_k);
    % could just use L_ss
    x_hat(k) = x_plus(k-1) + L(k)*(y(k)-C_k*x_plus(k-1));
    Q(k) = (1 - L(k)*C_k)*Q_plus(k-1)*(1 - L(k)*C_k)' + L(k)*V_k*L(k)';
    %state propagation
    x_plus(k) = A_k*x_hat(k);
    Q_plus(k) = A_k*Q(k)*A_k' + W_k;
    
    res(k) = y(k) - C_k*x_plus(k-1);
end
    
figure
plot(x,'-k','LineWidth',2); hold on;
plot(x_hat,'--r');
ylabel('Data');
xlabel('Time Steps');
legend({'$x$','$\hat{x}$'},'interpreter','latex');
title('Actual and Estimate');

figure
plot(res)
title('Residual');
xlabel('Time Steps');

figure
autocorr(res);

%% Part 3
% initiate
A_k = A; W_k = Q; C_k = C; V_k = R;
n = 11000;
v = wgn(n,1,1);
w = wgn(n,1,1);
x_hat(1) = 0;
x(1) = 0;
Q(1) = P;

% initial state propagation
x_plus(1) = A_k*x_hat(1);
Q_plus = A_k*Q(1)*A_k' + W_k;
y(1) = x(1) + v(1);
x(2) = A_k*x(1) + w(1);

res(1) = y(1) - C_k*x_plus(1);

for k = 2:n
    % regular system update
    x(k+1) = A_k*x(k) + w(k);
    y(k) = x(k) + v(k);
    
    % measurement update
    L(k) = Q_plus(k-1)*C_k'*inv(C_k*Q_plus(k-1)*C_k' + V_k);
    L(k) = 0.8*L(k);
    % could just use L_ss
    x_hat(k) = x_plus(k-1) + L(k)*(y(k)-C_k*x_plus(k-1));
    Q(k) = (1 - L(k)*C_k)*Q_plus(k-1)*(1 - L(k)*C_k)' + L(k)*V_k*L(k)';
    %state propagation
    x_plus(k) = A_k*x_hat(k);
    Q_plus(k) = A_k*Q(k)*A_k' + W_k;
    
    res(k) = y(k) - C_k*x_plus(k-1);
end
    
figure
plot(x); hold on;
plot(x_hat);
ylabel('Data');
xlabel('Time Steps');
legend({'$x$','$\hat{x}$'},'interpreter','latex');
title('Actual and Estimate');

figure
plot(res)
title('Residual');
xlabel('Time Steps');

figure
autocorr(res);  
