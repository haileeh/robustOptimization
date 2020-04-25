function Target = SetupTargetParams(mass, CM2DPDistance, CMToEdgeY_m, CMToEdgeZ_m, TargetShapeType, AngVelVec_rad)

% Coordinates
i = [1 0 0];
j = [0 1 0];
k = [0 0 1];

% Target mass
Target.mass_kg = mass;
% Distance from Target's CM to the docking port interface
Target.CM2DPinX_m = CM2DPDistance; 

% Distances from CM to edges
Target.CMToEdgeZ_m = CMToEdgeZ_m;
Target.CMToEdgeY_m = CMToEdgeY_m;
Target.CMToEdgeX_m = Target.CM2DPinX_m;

% Set up inertia matrix based on shape type
if strcmp(TargetShapeType, 'cylinder')
    Target.MOI_kgm2 = [1/2*Target.mass_kg*Target.CMToEdgeY_m^2, 0, 0;...
        0, 1/12*Target.mass_kg*(3*Target.CMToEdgeY_m^2 + Target.CMToEdgeX_m^2), 0;...
        0, 0, 1/12*Target.mass_kg*(3*Target.CMToEdgeY_m^2 + Target.CMToEdgeX_m^2)];
elseif strcmp(TargetShapeType, 'sphere')
    Target.MOI_kgm2 = (2/5)*Target.mass_kg*Target.CMToEdgeX_m^2*eye(3);
elseif strcmp(TargetShapeType, 'rect prism')
    Target.MOI_kgm2 = (1/12)*Target.mass_kg*...
        [((2*Target.CMToEdgeY_m)^2+(2*Target.CMToEdgeZ_m)^2), 0, 0;...
        0, ((2*Target.CMToEdgeZ_m)^2+(2*Target.CMToEdgeX_m)^2), 0;...
        0, 0, ((2*Target.CMToEdgeY_m)^2+(2*Target.CMToEdgeX_m)^2)];
end

% Simplified notation for following calcultions
Ixx = Target.MOI_kgm2(1,1); Iyy = Target.MOI_kgm2(2,2); Izz = Target.MOI_kgm2(3,3);

% Initial Conditions
% Initial angles should always be equal to 0 since the assumption is that
% the global frame = the body frame at the start of this phase; so they are
% not included as an option for ICs 
Target.roll1d_0 = AngVelVec_rad(1); %wx 
Target.pitch1d_0 = AngVelVec_rad(2); %wy
Target.yaw1d_0 = AngVelVec_rad(3); %wz 

% Initial angular accelerations
Target.roll1dd_0 = (Iyy-Izz)*Target.pitch1d_0*Target.yaw1d_0/Ixx;
Target.pitch1dd_0 = (Izz-Ixx)*Target.roll1d_0*Target.yaw1d_0/Iyy;
Target.yaw1dd_0 = (Ixx-Iyy)*Target.pitch1d_0*Target.roll1d_0/Izz;

% Angular velocity construction in Target
Target.BodyToInertialRot = FindRotMat(0, 0, 0);
omega_123 = [Target.roll1d_0; Target.pitch1d_0; Target.yaw1d_0];
Target.omega_xyz = Target.BodyToInertialRot*omega_123;

% Notation in inertial frame for determining ICs
%r1_xyz = [(-Target.CM2DPinX_m);0;0];
r1_xyz = [0;0;0];
v1_xyz = cross(Target.omega_xyz, r1_xyz);

Target.x1_0 = dot(r1_xyz,i);
Target.x1d_0 = dot(v1_xyz,i); % Change to 0 for no translation
Target.y1_0 = dot(r1_xyz,j);
Target.y1d_0 = dot(v1_xyz,j); % Change to 0 for no translation
Target.z1_0 = dot(r1_xyz,k);
Target.z1d_0 = dot(v1_xyz,k); % Change to 0 for no translation

Target.shape = TargetShapeType;

end