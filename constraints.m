function [constraint_matrix,B,b,Atilde,Btilde,Ctilde] = constraints(n,sys)
% Form the large constraint matrix for the SDP

Ax = sys.Ax_unc;
Bx = sys.Bx_unc;
Cx = sys.Cx; % same as Cx_unc
Qx = sys.Qx;
Rx = sys.Rx;
gamma = sys.gamma;

[~,n_u] = size(Bx);
n_x = length(Ax);
n_w = n_x;

%% Solve for the tilde matrices 
% Initialize the tilde matrices
Atilde = zeros(n_x,n_x,n);
Btilde = zeros(n_x,n*n_u,n);
Ctilde = zeros(n_x,n_w*n,n);
Atilde(:,:,1) = Ax;
Atilde(:,:,2) = Ax*Ax;
Btilde(:,:,1) = [Ax*Bx, zeros(n_x,n*n_u-2)];
Btilde(:,:,2) = [Ax*Ax*Bx, Bx, zeros(n_x,n*n_u-4)];
Ctilde(:,1:n_x,1) = Ax*Cx;
Ctilde(:,1:n_x,2) = Ax*Ax*Cx;
Ctilde(:,(n_x+1):2*n_x,2) = Cx;

% Form remaining columns of the tilde matrices
i = 2;
while(i<n)
    i=i+1;

    Ctilde(:,:,i) = Ctilde(:,:,i-1);
    Ctilde(:,1:n_x,i) = Atilde(:,:,i-1)*Cx; %was 1:2 
    Ctilde(:,4*i-3:4*i,i) = Cx;         % was 5:6 (2*i-1:2*i)
    Btilde(:,:,i) = Btilde(:,:,i-1);
    Btilde(:,1:2,i) = Atilde(:,:,i-1)*Bx;
    Btilde(:,2*i-1:2*i,i) = Bx;
    Atilde(:,:,i) = Atilde(:,:,i-1)*Ax;
end

%% Solve for the desired matrices in terms of tilde matrices
C = zeros(n_x*n); 
D=zeros(n_u*n,n_w*n); 
Btemp =zeros(n*n_u,n*n_u); 
BA=zeros(n*n_u,n_x);
CA = zeros(n_x*n,n_x);
for j=1:n
    C = C + Ctilde(:,:,j)'*Qx*Ctilde(:,:,j); 
    D = D + Btilde(:,:,j)'*Qx*Ctilde(:,:,j);
    Btemp = Btemp + Btilde(:,:,j)'*Qx*Btilde(:,:,j);
    BA = BA + Btilde(:,:,j)'*Qx*Atilde(:,:,j);
    CA = CA + Ctilde(:,:,j)'*Qx*Atilde(:,:,j);
end

c = CA*sys.IC;

Rhat = zeros(n*n_u);
for k = 1:n
    Rhat((k-1)*n_u +1:(k)*n_u, (k-1)*n_u +1:(k)*n_u) = Rx;
end
B = Rhat + Btemp;
b = BA*sys.IC;

h = c-D'*inv(B)*b;
F = B^(-1/2)*D;
% placeholder values for decision variables
yvec = ones(n*n_u,1); % decision variable
z = 5; %decision variable
lambda = 0; %decision variable
% Bring everything together in the desired constraint
constraint_matrix = ...  
    [eye(n*n_u), yvec, F;...    
     yvec', z-gamma^2*lambda, -h';...   
     F', -h, lambda*eye(n_w*n)-C+F'*F]; 

