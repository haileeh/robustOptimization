%% script for finding Target, Chaser EOM
% only applies forces using Chaser to change Target's dynamics
% except for roll

anim = 1;

TargetShapeType = 'rect prism';
Target = SetupTargetParams(5, 0.10, 0.35, 0.5, TargetShapeType, deg2rad([0, 0, 5]));
ChaserShapeType = TargetShapeType;
Chaser = SetupChaserParams(5, 0.10, 0.35, 0.5, ChaserShapeType, Target);

disp_vec = [-Target.CM2DPinX_m; 0; 0]; %from CMC to new point
Target.MOICombowrtTargetCM_kgm2 = Target.MOI_kgm2 + Chaser.MOI_kgm2 + ...
    (Target.mass_kg+Chaser.mass_kg)*(dot(disp_vec,disp_vec)*eye(3) - disp_vec*disp_vec');

Target.deltaT = 100;

%%%
% Target Nutation
Ixx = Target.MOI_kgm2(1,1); Iyy = Target.MOI_kgm2(2,2); Izz = Target.MOI_kgm2(3,3);
wx = Target.roll1d_0; wy = Target.pitch1d_0; wz = Target.yaw1d_0;
HT = sqrt((Ixx*wx)^2 + (Iyy*wy)^2 + (Izz*wz)^2);
gamma = acos(Ixx*wx/HT); % nutation angle
mu = atan2(Izz*wz,(Iyy*wy)); % precession angle

wxdot = (Iyy-Izz)*wy*wz/Ixx;
wydot = (Izz-Ixx)*wx*wz/Iyy;
wzdot = (Ixx-Iyy)*wx*wy/Izz;

mudot = (Izz*Iyy*(wzdot*wy-wydot*wz))/(Iyy^2*wy^2 + Izz^2*wz^2);
%%%

% Time Span
tf = 10*Target.deltaT;
h = 0.1; % step size
n = tf/h;
tspan = linspace(0,tf,n);

%% Discrete solution
%term_type = 'Constant'; %Options = Constant, PosLinear, NegLinear,...
% NegQuadratic
term_type = 'PosLinear';

[out,tspan] = AggregatedDynamicsDiscreteSolver(Target, Chaser, tspan, term_type);

%% Plot Angles, Angular Velocities, Angular Accelerations
plotAngleStates(tspan,out)

%% Plot Moments and Loads
plotLoadsMoments(tspan,out)

%% check quaternions
qw = out.qw; qx = out.qx; qy = out.qy; qz = out.qz;
n = length(qw);
for i=1:n
    q = [qw(i), qx(i), qy(i), qz(i)];
   [roll(i), pitch(i), yaw(i)] = quat2angle(q, 'XYZ');
end

%% Animation 3D
if anim
    animateDetumble(out, Chaser, Target);
end

