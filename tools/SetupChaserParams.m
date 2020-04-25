function Chaser = SetupChaserParams(mass, CM2DPDistance, CMToEdgeY_m, CMToEdgeZ_m, ChaserShapeType, Target)

% Coordinate frame
i = [1 0 0];
j = [0 1 0];
k = [0 0 1];

% Chaser mass
Chaser.mass_kg = mass;
% Distance from Chaser's CM to docking port interface
Chaser.CM2DPinX_m = CM2DPDistance;

% Distances from CM to edges
Chaser.CMToEdgeZ_m = CMToEdgeZ_m;
Chaser.CMToEdgeY_m = CMToEdgeY_m;
Chaser.CMToEdgeX_m = Chaser.CM2DPinX_m;

% Set up inertia matrix
if strcmp(ChaserShapeType, 'cylinder')
    Chaser.MOI_kgm2 = [1/2*Chaser.mass_kg*Chaser.CMToEdgeY_m^2, 0, 0;...
        0, 1/12*Chaser.mass_kg*(3*Chaser.CMToEdgeY_m^2 + Chaser.CMToEdgeX_m^2), 0;...
        0, 0, 1/12*Chaser.mass_kg*(3*Chaser.CMToEdgeY_m^2 + Chaser.CMToEdgeX_m^2)];
elseif strcmp(ChaserShapeType, 'sphere')
    Chaser.MOI_kgm2 = (2/5)*Chaser.mass_kg*Chaser.CMToEdgeX_m^2*eye(3);
elseif strcmp(ChaserShapeType, 'rect prism')
    Chaser.MOI_kgm2 = (1/12)*Chaser.mass_kg*...
        [((2*Chaser.CMToEdgeY_m)^2+(2*Chaser.CMToEdgeZ_m)^2), 0, 0;...
        0, ((2*Chaser.CMToEdgeZ_m)^2+(2*Chaser.CMToEdgeX_m)^2), 0;...
        0, 0, ((2*Chaser.CMToEdgeY_m)^2+(2*Chaser.CMToEdgeX_m)^2)];
end

% Set up initial conditions
R = (Target.CM2DPinX_m+Chaser.CM2DPinX_m);
r2_xyz = (Target.BodyToInertialRot*[R;0;0]);
v2_xyz = cross(Target.omega_xyz, r2_xyz);

Chaser.x2_0 = dot(r2_xyz,i) + Target.x1_0;
Chaser.x2d_0 = dot(v2_xyz,i) + Target.x1d_0;
Chaser.y2_0 = dot(r2_xyz,j) + Target.y1_0;
Chaser.y2d_0 = dot(v2_xyz,j)+ Target.y1d_0;
Chaser.z2_0 = dot(r2_xyz,k) + Target.z1_0;
Chaser.z2d_0 = dot(v2_xyz,k) + Target.z1d_0;
Chaser.shape = ChaserShapeType;

end