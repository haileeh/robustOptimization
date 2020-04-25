function [out,tspan] = AggregatedDynamicsDiscreteSolver(Target,Chaser,tspan,term_type)

n = length(tspan);
% initial conditions + setup
qw = zeros(2,1); qw(1) = 1;
qx = zeros(2,1); qy = zeros(2,1); qz = zeros(2,1);
qwdot = zeros(2,1); qxdot = zeros(2,1);
qydot = zeros(2,1); qzdot = zeros(2,1);
yaw1 = zeros(2,1); yaw1d = zeros(2,1);
yaw1d(1) = Target.yaw1d_0; yaw1dd = zeros(2,1); yaw1dd(1) = Target.yaw1dd_0;
pitch1 = zeros(2,1); pitch1d = zeros(2,1);
pitch1d(1) = Target.pitch1d_0; pitch1dd = zeros(2,1); pitch1dd(1) = Target.pitch1dd_0;
roll1 = zeros(2,1); roll1d = zeros(2,1);
roll1d(1) = Target.roll1d_0; roll1dd = zeros(2,1); roll1dd(1) = Target.roll1dd_0;
x1 = zeros(2,1); x1(1) = Target.x1_0; x1d = zeros(2,1);
x1d(1) = Target.x1d_0; x1dd = zeros(2,1); x1dd(1) = 0;
y1 = zeros(2,1); y1(1) = Target.y1_0; y1d = zeros(2,1);
y1d(1) = Target.y1d_0; y1dd = zeros(2,1); y1dd(1) = 0;
z1 = zeros(2,1); z1(1) = Target.z1_0; z1d = zeros(2,1);
z1d(1) = Target.z1d_0; z1dd = zeros(2,1); z1dd(1) = 0;
x2 = zeros(2,1); x2(1) = Chaser.x2_0; x2d = zeros(2,1);
x2d(1) = Chaser.x2d_0; x2dd = zeros(2,1); x2dd(1) = 0;
y2 = zeros(2,1); y2(1) = Chaser.y2_0; y2d = zeros(2,1);
y2d(1) = Chaser.y2d_0; y2dd = zeros(2,1); y2dd(1) = 0;
z2 = zeros(2,1); z2(1) = Chaser.z2_0; z2d = zeros(2,1);
z2d(1) = Chaser.z2d_0; z2dd = zeros(2,1); z2dd(1) = 0;
BodyToInertialRot = zeros(3,3,n); BodyToInertialRot(:,:,1) = FindRotMat(0, 0, 0);
CurrentB2IRot = BodyToInertialRot;
CentAccInBody_mpsps = zeros(3,2);
AngAccInBody_mpsps = zeros(3,2);
Mx = zeros(2,1); Fz = zeros(2,1); Fth = zeros(2,1);
if isfield(Chaser,'roll1dd_user')
    roll1dd_user = Chaser.roll1dd_user;
    pitch1dd_user = Chaser.pitch1dd_user;
    yaw1dd_user = Chaser.yaw1dd_user;
end

t_wait = 1; 
i = 1;
tol = 0.001;
h = tspan(2)-tspan(1); %assumes linspace
R = (Target.CM2DPinX_m + Chaser.CM2DPinX_m);
Ixx = Target.MOI_kgm2(1,1); Iyy = Target.MOI_kgm2(2,2); Izz = Target.MOI_kgm2(3,3);
rollstopped = 0;
RInBody = R*[1;0;0];
    
while (abs(yaw1d(i)) >= tol || abs(pitch1d(i)) >= tol || abs(roll1d(i)) >= tol)
    
    % Determine double dot terms
    if i~=1
        roll1dd(i) = (Iyy-Izz)*pitch1d(i)*yaw1d(i)/Ixx;
        pitch1dd(i) = (Izz-Ixx)*roll1d(i)*yaw1d(i)/Iyy;
        yaw1dd(i) = (Ixx-Iyy)*pitch1d(i)*roll1d(i)/Izz;
    end
    xdd_trans = 0; ydd_trans = 0; zdd_trans = 0;
    if tspan(i) >= t_wait % start the process of stopping the Target
        %If nutating, we want to stop roll1d first before stopping yaw1d,
        %pitch1d
        if strcmp(term_type, 'user_defined')
                roll1dd(i) = roll1dd_user(i) + roll1dd(i);
                pitch1dd(i) = pitch1dd_user(i) + pitch1dd(i);
                yaw1dd(i) = yaw1dd_user(i) + yaw1dd(i);
        else
            if (Ixx == Iyy) && (Iyy == Izz)
                % Target is a sphere; no nutation
                if strcmp(term_type, 'NegQuadratic')
                    yaw1dd_slowdown = NegativeQuadraticTerm(Target,...
                        Target.yaw1d_0, t_wait, tspan(i));
                    roll1dd_slowdown = NegativeQuadraticTerm(Target,...
                        Target.roll1d_0, t_wait, tspan(i));
                    pitch1dd_slowdown = NegativeQuadraticTerm(Target,...
                        Target.pitch1d_0, t_wait, tspan(i));
                elseif strcmp(term_type, 'PosLinear')
                    yaw1dd_slowdown = PositiveLinearTerm(Target,...
                        Target.yaw1d_0, t_wait, tspan(i));
                    roll1dd_slowdown = PositiveLinearTerm(Target,...
                        Target.roll1d_0, t_wait, tspan(i));
                    pitch1dd_slowdown = PositiveLinearTerm(Target,...
                        Target.pitch1d_0, t_wait, tspan(i));
                elseif strcmp(term_type, 'NegLinear')
                    yaw1dd_slowdown = NegativeLinearTerm(Target,...
                        Target.yaw1d_0, t_wait, tspan(i));
                    roll1dd_slowdown = NegativeLinearTerm(Target,...
                        Target.roll1d_0, t_wait, tspan(i));
                    pitch1dd_slowdown = NegativeLinearTerm(Target,...
                        Target.pitch1d_0, t_wait, tspan(i));
                elseif strcmp(term_type, 'Constant')
                    yaw1dd_slowdown = ConstantAccTerm(Target,...
                        Target.yaw1d_0);
                    roll1dd_slowdown = ConstantAccTerm(Target,...
                        Target.roll1d_0);
                    pitch1dd_slowdown = ConstantAccTerm(Target,...
                        Target.pitch1d_0);
                end
                roll1dd(i) = roll1dd_slowdown + roll1dd(i);
                pitch1dd(i) = pitch1dd_slowdown + pitch1dd(i);
                yaw1dd(i) = yaw1dd_slowdown + yaw1dd(i);
            else
                %Target may be nutating
                if abs(roll1d(i)) > tol/10
                    if strcmp(term_type, 'NegQuadratic')
                        roll1dd_slowdown = NegativeQuadraticTerm(Target,...
                            Target.roll1d_0, t_wait, tspan(i));
                    elseif strcmp(term_type, 'PosLinear')
                        roll1dd_slowdown = PositiveLinearTerm(Target,...
                            Target.roll1d_0, t_wait, tspan(i));
                    elseif strcmp(term_type, 'NegLinear')
                        roll1dd_slowdown = NegativeLinearTerm(Target,...
                            Target.roll1d_0, t_wait, tspan(i));
                    elseif strcmp(term_type, 'Constant')
                        roll1dd_slowdown = ConstantAccTerm(Target,...
                            Target.roll1d_0);
                    end
                    roll1dd(i) = roll1dd_slowdown + roll1dd(i);
                end
                if abs(roll1d(i)) <= tol/10
                    % slow down ang vel in other axes now
                    if rollstopped == 0
                        % save off these ang vels
                        Target.yaw1d_0 = yaw1d(i);
                        Target.pitch1d_0 = pitch1d(i);
                    end
                    if abs(pitch1d(i)) > tol
                        if strcmp(term_type, 'NegQuadratic')
                            pitch1dd_slowdown = NegativeQuadraticTerm(Target,...
                                Target.pitch1d_0, t_wait, tspan(i));
                        elseif strcmp(term_type, 'PosLinear')
                            pitch1dd_slowdown = PositiveLinearTerm(Target,...
                                Target.pitch1d_0, t_wait, tspan(i));
                        elseif strcmp(term_type, 'NegLinear')
                            pitch1dd_slowdown = NegativeLinearTerm(Target,...
                                Target.pitch1d_0, t_wait, tspan(i));
                        elseif strcmp(term_type, 'Constant')
                            pitch1dd_slowdown = ConstantAccTerm(Target,...
                                Target.pitch1d_0);
                        end
                        pitch1dd(i) = pitch1dd_slowdown + pitch1dd(i);
                    end
                    if abs(yaw1d(i)) > tol
                        if strcmp(term_type, 'NegQuadratic')
                            yaw1dd_slowdown = NegativeQuadraticTerm(Target,...
                                Target.yaw1d_0, t_wait, tspan(i));
                        elseif strcmp(term_type, 'PosLinear')
                            yaw1dd_slowdown = PositiveLinearTerm(Target,...
                                Target.yaw1d_0, t_wait, tspan(i));
                        elseif strcmp(term_type, 'NegLinear')
                            yaw1dd_slowdown = NegativeLinearTerm(Target,...
                                Target.yaw1d_0, t_wait, tspan(i));
                        elseif strcmp(term_type, 'Constant')
                            yaw1dd_slowdown = ConstantAccTerm(Target,...
                                Target.yaw1d_0);
                        end
                        yaw1dd(i) = yaw1dd_slowdown + yaw1dd(i);
                    end
                    rollstopped = 1;
                end
            end
        end
        
        % Always occur post t_wait
        if abs(y1d(i)) > tol
            ydd_trans = ConstantAccTerm(Target, Target.y1d_0);
        end
        if abs(x1d(i)) > tol
            xdd_trans = ConstantAccTerm(Target, Target.x1d_0);
        end
        if abs(z1d(i)) > tol
            zdd_trans = ConstantAccTerm(Target, Target.z1d_0);
        end
    end

    % Quaternion to DCM
    if i~=1
        BodyToInertialRot(:,:,i) = quat2dcm([qw(i), qx(i), qy(i), qz(i)]);
        BodyToInertialRot(:,:,i) = BodyToInertialRot(:,:,i)';
        CurrentB2IRot(:,:,i) =  BodyToInertialRot(:,:,i);
    end
    omega_123 = [roll1d(i); pitch1d(i); yaw1d(i)];
    alpha_123 = [roll1dd(i); pitch1dd(i); yaw1dd(i)];
    
    % Necessary for maintaining synchronous motion
    CentAccInBody_mpsps(:,i) = cross(omega_123, cross(omega_123, RInBody));
    AngAccInBody_mpsps(:,i) = cross(alpha_123, RInBody);
    
    CentAccInInertial_mpsps = CurrentB2IRot(:,:,i)*CentAccInBody_mpsps(:,i);
    AngAccInInertial_mpsps = CurrentB2IRot(:,:,i)*AngAccInBody_mpsps(:,i);
    TotalAccInInertial_mpsps = CentAccInInertial_mpsps + AngAccInInertial_mpsps;
    
    x2dd(i) = TotalAccInInertial_mpsps(1) + xdd_trans;
    y2dd(i) = TotalAccInInertial_mpsps(2) + ydd_trans;
    z2dd(i) = TotalAccInInertial_mpsps(3) + zdd_trans;
    
    x1dd(i) = xdd_trans;
    y1dd(i) = ydd_trans;
    z1dd(i) = zdd_trans;
    
    qwdot(i) = -1/2 * (roll1d(i).*qx(i) + pitch1d(i).*qy(i) + yaw1d(i).*qz(i));
    qxdot(i) = 1/2 * (roll1d(i).*qw(i) + pitch1d(i).*qz(i) - yaw1d(i).*qy(i));
    qydot(i) = 1/2 * (pitch1d(i).*qw(i) + yaw1d(i).*qx(i) - roll1d(i).*qz(i));
    qzdot(i) = 1/2 * (yaw1d(i).*qw(i) + roll1d(i).*qy(i) - pitch1d(i).*qx(i));
    
    % results from TARGET AMB (wrt aggregated system)
    Hdot = Target.MOICombowrtTargetCM_kgm2*alpha_123 + ...
        cross(omega_123, Target.MOICombowrtTargetCM_kgm2*omega_123);
    Mx(i) = Hdot(1);
    Fz(i) = -Hdot(2)/R;
    Fth(i) = Hdot(3)/R;
    
    % Update timer
    i = i+1;
    if i>n
        break;
    end
    
    % update velocity terms
    x1d(i) = x1d(i-1) + x1dd(i-1)*h; y1d(i) = y1d(i-1) + y1dd(i-1)*h;
    z1d(i) = z1d(i-1) + z1dd(i-1)*h;
    nextstep = CurrentB2IRot(:,:,i-1)*cross(omega_123, RInBody);
    x2d(i) = x1d(i) + nextstep(1); y2d(i) = y1d(i) + nextstep(2);
    z2d(i) = z1d(i) + nextstep(3);
    
    % Update position terms
    x1(i) = x1(i-1) + x1d(i)*h; y1(i) = y1(i-1) + y1d(i)*h;
    z1(i) = z1(i-1) + z1d(i)*h;
    nextstep = CurrentB2IRot(:,:,i-1)*RInBody;
    x2(i) = x1(i) + nextstep(1); y2(i) = y1(i) + nextstep(2);
    z2(i) = z1(i) + nextstep(3);
    
    % update angular velocity
    yaw1d(i) = yaw1d(i-1) + yaw1dd(i-1)*h;
    pitch1d(i) = pitch1d(i-1) + pitch1dd(i-1)*h;
    roll1d(i) = roll1d(i-1) + roll1dd(i-1)*h;
    
    %update angles, delta angle
    yaw1(i) = yaw1(i-1) + yaw1d(i)*h; 
    pitch1(i) = pitch1(i-1) + pitch1d(i)*h;
    roll1(i) = roll1(i-1) + roll1d(i)*h;
    
    %update quaternion
    qw(i) = qw(i-1) + qwdot(i-1)*h; qx(i) = qx(i-1) + qxdot(i-1)*h;
    qy(i) = qy(i-1) + qydot(i-1)*h; qz(i) = qz(i-1) + qzdot(i-1)*h;
    
end
if i<=n
    tspan = tspan(1:i);
end

% Pack output
out.roll1 = roll1;      out.pitch1 = pitch1;        out.yaw1 = yaw1;
out.roll1d = roll1d;    out.pitch1d = pitch1d;      out.yaw1d = yaw1d;
out.roll1dd = roll1dd;  out.pitch1dd = pitch1dd;    out.yaw1dd = yaw1dd;
out.x1 = x1;            out.y1 = y1;                out.z1 = z1;
out.x2 = x2;            out.y2 = y2;                out.z2 = z2;
out.x1d = x1d;          out.y1d = y1d;              out.z1d = z1d;
out.x2d = x2d;          out.y2d = y2d;              out.z2d = z2d;
out.x1dd = x1dd;        out.y1dd = y1dd;            out.z1dd = z1dd;
out.x2dd = x2dd;        out.y2dd = y2dd;            out.z2dd = z2dd;
out.CurrentB2IRot = CurrentB2IRot;
out.qw = qw; out.qx = qx; out.qy = qy; out.qz = qz;

%% Calculate loads/moments
TotalForceInBody_N = Chaser.mass_kg*(CentAccInBody_mpsps + ...
    AngAccInBody_mpsps);
out.ShearLoadY_N = TotalForceInBody_N(2,:) - Fth';
out.ShearLoadZ_N = TotalForceInBody_N(3,:) - Fz';

out.Torsion_Nm = Chaser.MOI_kgm2(1,1)*roll1dd - ...
    (Chaser.MOI_kgm2(2,2)-Chaser.MOI_kgm2(3,3)).*...
    pitch1d(1:length(roll1dd)).*yaw1d(1:length(roll1dd)) - Mx;
out.BendingMomY_Nm = Chaser.MOI_kgm2(2,2)*pitch1dd - ...
    (Chaser.MOI_kgm2(3,3)-Chaser.MOI_kgm2(1,1)).*...
    roll1d(1:length(pitch1dd)).*yaw1d(1:length(pitch1dd));
out.BendingMomZ_Nm = Chaser.MOI_kgm2(3,3)*yaw1dd - ...
    (Chaser.MOI_kgm2(1,1)-Chaser.MOI_kgm2(2,2)).*...
    pitch1d(1:length(yaw1dd)).*roll1d(1:length(yaw1dd));

p = length(out.ShearLoadY_N);
out.int.ShearLoadY_N = cumtrapz(tspan(1:p), abs(out.ShearLoadY_N));
out.int.ShearLoadZ_N = cumtrapz(tspan(1:p), abs(out.ShearLoadZ_N));
out.int.Torsion_Nm = cumtrapz(tspan(1:p), abs(out.Torsion_Nm));
out.int.BendingMomY_Nm = cumtrapz(tspan(1:p), abs(out.BendingMomY_Nm));
out.int.BendingMomZ_Nm = cumtrapz(tspan(1:p), abs(out.BendingMomZ_Nm));
out.int.TotalLM = out.int.ShearLoadY_N(end) + out.int.ShearLoadZ_N(end) + ...
    out.int.Torsion_Nm(end) + out.int.BendingMomY_Nm(end)...
    + out.int.BendingMomZ_Nm(end);

%% Calculate deltaV
n = length(yaw1dd);
t = tspan(1:n);

out.DeltaVRoll_mps = cumtrapz(t, roll1dd);
out.DeltaVPitch_mps = cumtrapz(t, pitch1dd);
out.DeltaVYaw_mps = cumtrapz(t, yaw1dd);
out.TotalDeltaV_mps = out.DeltaVRoll_mps(end) + ...
    out.DeltaVPitch_mps(end) + out.DeltaVYaw_mps(end);

end

function current_acc = ConstantAccTerm(p,init_vel)
current_acc = (0 - init_vel) / p.deltaT;
end

function current_acc = NegativeLinearTerm(p,init_vel,t_wait,current_time)
a = init_vel/(p.deltaT)^2;
current_acc = -2*a*(current_time-t_wait);
end

function current_acc = PositiveLinearTerm(p,init_vel,t_wait,current_time)
a = init_vel/(p.deltaT)^2;
current_acc = 2*a*(current_time-(p.deltaT+t_wait));
end

function current_acc = NegativeQuadraticTerm(p,init_vel,t_wait,current_time)
a = init_vel/(p.deltaT)^3;
current_acc = -3*a*(current_time-t_wait)^2;
end
