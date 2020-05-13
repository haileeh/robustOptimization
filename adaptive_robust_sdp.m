function u = adaptive_robust_sdp(sys,N,iter,input_u)

import mosek.fusion.*;

n_u = 2;
n_w = 4;
n_x = 4;

%sys.Ax = sys.Ax + eye(size(sys.Ax,2));
sys.Ax_unc = sys.Ax_unc + eye(size(sys.Ax_unc,2));

gamma = sys.gamma;%10; % user selected: can't solve 0.001 or 0.01 or 0.1
[constraint_matrix,B,b,Atilde,Btilde,Ctilde] = constraints(gamma,N,sys);
dim = size(constraint_matrix,1);

% SDP
mod = Model('SDO_horizon');
% constraint
psdMat = mod.variable('psdMat', Domain.inPSDCone(dim));
sx = mod.variable('sx', Domain.inRange(-10.,3.)); % might want to change this
sy = mod.variable('sy', Domain.inRange(-10.,3.));

% decision variables
lambda = mod.variable('lambda', Domain.greaterThan(0.)); 
z = mod.variable('z'); 

% define s constraint
worst_w = sys.w_mag*ones(16,1); %equal to gamma? or some relationship to gamma
if iter==0
    old_u = [1; -0.7];%[0.1; -0.07]; % this is my guess
    old_u = [old_u; old_u; old_u; old_u];
else
    old_u = input_u;
end
temp = Atilde(:,:,end)*sys.IC+Btilde(:,:,end)*old_u+Ctilde(:,:,end)*worst_w;
sfactor = 0.5;
% conSx = temp(2)-sfactor*temp(1);
conSx = Expr.sub(temp(2), Expr.mul(0.1,Expr.mul(lambda,temp(1))));
conSy = Expr.sub(temp(4), Expr.mul(0.1,Expr.mul(lambda,temp(3))));


% PSD Constraint
% (1,1) = (1:8,1:8)
for i=1:N*n_u
    for j=1:N*n_u
        psd_i = psdMat.index([i-1,j-1]);
        mod.constraint(psd_i, Domain.equalsTo(constraint_matrix(i,j)));
    end
end

% (1,2)

% (1,3) = (1:8,10:25)
for i = 1:N*n_u
    for j = N*n_u+2 : N*n_u+1 + N*n_w
        psd_i = psdMat.index([i-1,j-1]);
        mod.constraint(psd_i, Domain.equalsTo(constraint_matrix(i, j)) );  
    end
end

% (2,1)

% (2,2)
i = N*n_u + 1;
j = N*n_u + 1;
psd_i = psdMat.index([i-1,j-1]);
constraint = Expr.sub(z, Expr.mul(gamma^2, lambda));
constraint = Expr.sub(constraint, conSx);
constraint = Expr.sub(constraint, conSy);
mod.constraint(Expr.sub ( psd_i, constraint), Domain.equalsTo(0.0));


% (2,3)
i = N*n_u + 1;
for j = N*n_u + 2 : N*n_u+1 + N*n_w
    psd_i = psdMat.index([i-1,j-1]);
    mod.constraint( psd_i, Domain.equalsTo(constraint_matrix(i, j)) ); 
end

% (3,1)
for i = N*n_u+2 : N*n_u+1 + N*n_w
    for j = 1:N*n_u
        psd_i = psdMat.index([i-1,j-1]);
        mod.constraint( psd_i, Domain.equalsTo(constraint_matrix(i, j)) ); 
    end
end

% (3,2)
for i = N*n_u + 2 : N*n_u+1 + N*n_w
    j = N*n_u + 1;
    psd_i = psdMat.index([i-1,j-1]);
    mod.constraint( psd_i, Domain.equalsTo(constraint_matrix(i, j)) );
end

% (3,3)
for i = N*n_u + 2 : N*n_u+1 + N*n_w
    for j = N*n_u + 2 : N*n_u+1 + N*n_w
        psd_i = psdMat.index([i-1,j-1]);
        if i==j
            constraint = Expr.add(lambda, constraint_matrix(i, j));
            mod.constraint(Expr.sub( psd_i, constraint), Domain.equalsTo( 0.0 ));
        else
            mod.constraint( psd_i, Domain.equalsTo(constraint_matrix(i, j)) );
        end
    end
end

mod.objective(ObjectiveSense.Minimize, z);

mod.solve();

PSD_answer = psdMat.level();

y = PSD_answer(dim*N*n_u+1 : dim*N*n_u + N*n_u);

% convert to u
u = B^(-1/2)*y - B\b;
