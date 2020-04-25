%SDP approach
clear all
% Time Span
tf = 1;%100;
h = .1;%0.01; % step size
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

% x-direction
Ax = [0 1; cx/mx, bx/mx];
Bx = [0 1/mx]';
Qx = [1 0; 0 100];
Rx = 0.1;
% Kx = dlqr(Ax,Bx,Qx,Rx,[]);
Cx = eye(2);
Atilde = zeros(2,2,n);
Btilde = zeros(2,n,n);
Ctilde = zeros(2,2*n,n);
Atilde(:,:,1) = Ax;
Atilde(:,:,2) = Ax*Ax;
Btilde(:,:,1) = [Ax*Bx, zeros(2,n-1)];
Btilde(:,:,2) = [Ax*Ax*Bx, Bx, zeros(2,n-2)];
Ctilde(:,1:2,1) = Ax*Cx;
Ctilde(:,1:2,2) = Ax*Ax*Cx;
Ctilde(:,3:4,2) = Cx;
i = 2;
while(i<n)
    i=i+1;

    Ctilde(:,:,i) = Ctilde(:,:,i-1);
    Ctilde(:,1:2,i) = Atilde(:,:,i-1)*Cx;
    Ctilde(:,2*i-1:2*i,i) = Cx;
    Btilde(:,:,i) = Btilde(:,:,i-1);
    Btilde(:,1,i) = Atilde(:,:,i-1)*Bx;
    Btilde(:,i,i) = Bx;
    Atilde(:,:,i) = Atilde(:,:,i-1)*Ax;
end

C = zeros(2*n); D=zeros(n,2*n); Btemp =zeros(n,n); BA=zeros(n,2);
for j=1:n-1
    C = C + Ctilde(:,:,j)'*Qx*Ctilde(:,:,j);
    D = D + Btilde(:,:,j)'*Qx*Ctilde(:,:,j);
    Btemp = Btemp + Btilde(:,:,j)'*Qx*Btilde(:,:,j);
    BA = BA + Btilde(:,:,j)'*Qx*Atilde(:,:,j);
end
c = C*x(1); %this should have the xdot initial velocity too
Rhat = Rx*eye(n);
B = Rhat + Btemp;
b = BA*x(1);