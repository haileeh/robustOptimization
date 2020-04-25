function animateDetumble(out, Chaser, Target)

% unpack array
yaw1 = (out.yaw1); pitch1 = (out.pitch1); roll1 = (out.roll1);
x1 = out.x1; y1 = out.y1; z1 = out.z1;
x2 = out.x2; y2 = out.y2; z2 = out.z2;
ChaserShapeType = Chaser.shape;
TargetShapeType = Target.shape;

i = [1 0 0];
j = [0 1 0];
k = [0 0 1];

CurrentB2IRot(:,:,1) = eye(3);

if strcmp(ChaserShapeType, 'sphere')
    [x,y,z] = sphere;
    x_Chaser = x2(1) + Chaser.CM2DPinX_m.*x;
    y_Chaser = y2(1) + Chaser.CM2DPinX_m.*y;
    z_Chaser = z2(1) + Chaser.CM2DPinX_m.*z;
elseif strcmp(ChaserShapeType, 'cylinder')
    [zc,yc,xc] = cylinder(Chaser.CMToEdgeY_m);
    x_Chaser = [x2(1)-Chaser.CMToEdgeX_m*ones(1,21); x2(1)+Chaser.CMToEdgeX_m*ones(1,21)];
    y_Chaser = y2(1) + yc;
    z_Chaser = z2(1) + zc;
    h_c = figure;
    S = surf(x_Chaser, y_Chaser, z_Chaser);
    fvc2 = surf2patch(S.XData,S.YData,S.ZData);
    close(h_c);
elseif strcmp(ChaserShapeType, 'rect prism')
    % axisymmetric rectangular prism; assumes center of mass is in the middle
    [x_Chaser, y_Chaser, z_Chaser]=plot_rect3D(Chaser.CMToEdgeX_m, Chaser.CMToEdgeY_m, CurrentB2IRot(:,:,1),[x2(1) y2(1) z2(1)]);
end

if strcmp(TargetShapeType, 'sphere')
    [x,y,z]=sphere;
    x_Target = x1(1) + Target.CM2DPinX_m.*x;
    y_Target = y1(1) + Target.CM2DPinX_m.*y;
    z_Target = z1(1) + Target.CM2DPinX_m.*z;
elseif strcmp(TargetShapeType, 'cylinder')
    [zt,yt,xt] = cylinder(Target.CMToEdgeY_m);
    x_Target = [x1(1)-Target.CMToEdgeX_m*ones(1,21); x1(1)+Target.CMToEdgeX_m*ones(1,21)];
    y_Target = y1(1) + yt;
    z_Target = z1(1) + zt;
    h_t = figure;
    S = surf(x_Target, y_Target, z_Target);
    fvc1 = surf2patch(S.XData,S.YData,S.ZData);
    close(h_t);
elseif strcmp(TargetShapeType, 'rect prism')
    [x_Target, y_Target, z_Target]=plot_rect3D(Target.CMToEdgeX_m, Target.CMToEdgeY_m, CurrentB2IRot(:,:,1),[x1(1) y1(1) z1(1)]);
end

% Initialize figure
figure; grid on;
set(gcf, 'Position', [100, 100, 1000, 600])
xl = xlabel('$x(t) [km]$','Interpreter','LaTeX');
yl = ylabel('$y(t) [km]$','Interpreter','LaTeX');
zl = zlabel('$z(t) [km]$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(zl,'FontSize',18);
tl = title('Trajectory');
set(tl,'FontSize',18);

axis([-2 2 -2 2 -2 2]);

% Initialize patches for Target and Chaser
if strcmp(TargetShapeType, 'sphere') || strcmp(TargetShapeType, 'rect prism')
    TargetPatch = patch(x_Target,y_Target, z_Target, [ 0.9100 0.4100 0.1700]);
    TargetPatch.EdgeColor = 'k';
elseif strcmp(TargetShapeType, 'cylinder')
    TargetPatch = patch('Faces',fvc1.faces,'Vertices', fvc1.vertices,...
        'FaceColor', [ 0.9100 0.4100 0.1700]);
    TargetPatch.EdgeColor = 'k';
    % Create lids
    [TargetLid1,TargetLid2]=AddLid2Cyl([x1(1)-Target.CMToEdgeX_m; x1(1)+Target.CMToEdgeX_m],y1(1),z1(1),Target.CMToEdgeY_m,'Target');
    TargetLid1OG = TargetLid1; TargetLid2OG = TargetLid2;
end

if strcmp(ChaserShapeType, 'sphere') || strcmp(ChaserShapeType, 'rect prism')
    ChaserPatch = patch(x_Chaser, y_Chaser, z_Chaser,'b');
    ChaserPatch.EdgeColor = 'k';
elseif strcmp(ChaserShapeType, 'cylinder')
    ChaserPatch = patch('Faces',fvc2.faces,'Vertices', fvc2.vertices, 'FaceColor','b');
    ChaserPatch.EdgeColor = 'k';
    % Create lids
    [ChaserLid1,ChaserLid2]=AddLid2Cyl([x2(1)-Chaser.CMToEdgeX_m; x2(1)+Chaser.CMToEdgeX_m],y2(1),z2(1),Chaser.CMToEdgeY_m,'Chaser');
    ChaserLid1OG = ChaserLid1; ChaserLid2OG = ChaserLid2;
end

alpha(0.5)

% Initialize patches for docking ports
dp_123 = [(Target.CM2DPinX_m+Chaser.CM2DPinX_m); 0; 0];
dp_xyz(:,1) = CurrentB2IRot(:,:,1)*dp_123;
DP1 = patch([x1(1) (x1(1) + dot(dp_xyz(:,1),i))], ...
    [y1(1) (y1(1)+dot(dp_xyz(:,1),j))],...
    [z1(1) (z1(1)+dot(dp_xyz(:,1),k))], 'r', 'LineWidth', 3);
DP1.EdgeColor = 'r';
DP2 = patch([x2(1) (x2(1) - dot(dp_xyz(:,1),i))], ...
    [y2(1) (y2(1)-dot(dp_xyz(:,1),j))],...
    [z2(1) (z2(1)-dot(dp_xyz(:,1),k))], 'r', 'LineWidth', 3);
DP2.EdgeColor = 'r';

waitforbuttonpress
BodyToInertialRot(:,:,1) = CurrentB2IRot(:,:,1);

q = 10;
for l=(q+1):q:length(yaw1)
    
    BodyToInertialRot(:,:,l) = FindRotMat(roll1(l)-roll1(l-q),...
        pitch1(l)-pitch1(l-q), yaw1(l)-yaw1(l-q));
    CurrentB2IRot(:,:,l) = BodyToInertialRot(:,:,l)*CurrentB2IRot(:,:,l-q);
    
    if strcmp(TargetShapeType, 'sphere')
        TargetPatch.XData = x1(l) + Target.CM2DPinX_m.*x;
        TargetPatch.YData = y1(l) + Target.CM2DPinX_m.*y;
        TargetPatch.ZData = z1(l) + Target.CM2DPinX_m.*z;
    elseif strcmp(TargetShapeType, 'cylinder')
        if l==q+1
            x_Target = [-Target.CMToEdgeX_m*ones(1,21); Target.CMToEdgeX_m*ones(1,21)];
            h_t = figure;
            S = surf(x_Target, y_Target, z_Target);
            fvc1 = surf2patch(S.XData,S.YData,S.ZData);
            close(h_t);
        end
        
        for b=1:length(fvc1.vertices)
            vert_out_t(b,:) =  CurrentB2IRot(:,:,l)*(fvc1.vertices(b,:))';
        end
        TargetPatch.Vertices = vert_out_t + [x1(l), y1(l), z1(l)];
        
        % Update Lid 1 - Interpolate patch data
        [r,~]=size(TargetLid1OG.Vertices);
        v = TargetPatch.Vertices(1:2:end,1);
        xq = linspace(1,length(v),r);
        for pp=1:3
            TargetLid1.Vertices(:,pp) = interp1(TargetPatch.Vertices(1:2:end,pp),xq);
        end
        
        % Update Lid 2 - Interpolate patch data
        [r,c]=size(TargetLid2OG.Vertices);
        xq = linspace(1,length(v),r);
        for pp=1:3
            TargetLid2.Vertices(:,pp) = interp1(TargetPatch.Vertices(2:2:end,pp),xq);
        end
        
    elseif strcmp(TargetShapeType, 'rect prism')
        [TargetPatch.XData, TargetPatch.YData, TargetPatch.ZData]=plot_rect3D(Target.CMToEdgeX_m, Target.CMToEdgeY_m, CurrentB2IRot(:,:,l),[x1(l) y1(l) z1(l)]);
    end
    hold on;
    plot3(x1(l), y1(l), z1(l), '.k');
    
    if strcmp(ChaserShapeType, 'sphere')
        ChaserPatch.XData = x2(l) + Chaser.CM2DPinX_m.*x;
        ChaserPatch.YData = y2(l) + Chaser.CM2DPinX_m.*y;
        ChaserPatch.ZData = z2(l) + Chaser.CM2DPinX_m.*z;
    elseif strcmp(ChaserShapeType, 'cylinder')
        if l==q+1
            x_Chaser = [-Chaser.CMToEdgeX_m*ones(1,21); Chaser.CMToEdgeX_m*ones(1,21)];
            h_c = figure;
            S = surf(x_Chaser, y_Chaser, z_Chaser);
            fvc2 = surf2patch(S.XData,S.YData,S.ZData);
            close(h_c);
        end
        
        for b=1:length(fvc2.vertices)
            vert_out_c(b,:) =  CurrentB2IRot(:,:,l)*(fvc2.vertices(b,:))';
        end
        ChaserPatch.Vertices = vert_out_c + [x2(l), y2(l), z2(l)];
        
        % Update Lid 1 - Interpolate patch data
        [r,~]=size(ChaserLid1OG.Vertices);
        v = ChaserPatch.Vertices(1:2:end,1);
        xq = linspace(1,length(v),r);
        for pp=1:3
            ChaserLid1.Vertices(:,pp) = interp1(ChaserPatch.Vertices(1:2:end,pp),xq);
        end
        
        % Update Lid 2 - Interpolate patch data
        [r,c]=size(ChaserLid2OG.Vertices);
        xq = linspace(1,length(v),r);
        for pp=1:3
            ChaserLid2.Vertices(:,pp) = interp1(ChaserPatch.Vertices(2:2:end,pp),xq);
        end
        
    elseif strcmp(ChaserShapeType, 'rect prism')
        [ChaserPatch.XData, ChaserPatch.YData, ChaserPatch.ZData]=plot_rect3D(Chaser.CMToEdgeX_m, Chaser.CMToEdgeY_m, CurrentB2IRot(:,:,l),[x2(l) y2(l) z2(l)]);
    end
    plot3(x2(l), y2(l), z2(l), '.g');
    
    dp_xyz(:,l) = CurrentB2IRot(:,:,l)*dp_123;
    DP1.XData = [x1(l) (x1(l) + dot(dp_xyz(:,l),i))];
    DP1.YData = [y1(l) (y1(l) + dot(dp_xyz(:,l),j))];
    DP1.ZData = [z1(l) (z1(l) + dot(dp_xyz(:,l),k))];
    hold on;
    %plot3((x1(l) + dot(dp_xyz(:,l),i)), (y1(l) + dot(dp_xyz(:,l),j)), (z1(l) + dot(dp_xyz(:,l),k)),'.r');
    
    DP2.XData = [x2(l) (x2(l) - dot(dp_xyz(:,l),i))];
    DP2.YData = [y2(l) (y2(l) - dot(dp_xyz(:,l),j))];
    DP2.ZData = [z2(l) (z2(l) - dot(dp_xyz(:,l),k))];
    
    drawnow;
end