function [XData, YData, ZData]=plot_rect3D(height, radius, RotMat,r0)
% Plotting a 3D rectangular prism about its center of mass

htdiv2 = height;%height/2;
raddiv2 = radius;%radius/2;

Face1 = [htdiv2, htdiv2, htdiv2, htdiv2;...
    -raddiv2, raddiv2, raddiv2, -raddiv2;...
    -raddiv2, -raddiv2, raddiv2, raddiv2];

Face2 = [htdiv2, -htdiv2, -htdiv2, htdiv2;...
    raddiv2, raddiv2, -raddiv2, -raddiv2;...
    raddiv2, raddiv2, raddiv2, raddiv2];

Face3 = [-htdiv2,  -htdiv2,  -htdiv2,  -htdiv2;...
    raddiv2, raddiv2,  -raddiv2,  -raddiv2;...
    -raddiv2, raddiv2, raddiv2, -raddiv2];

Face4 = [htdiv2,  -htdiv2,  -htdiv2, htdiv2;...
    raddiv2, raddiv2,  -raddiv2,  -raddiv2;...
    -raddiv2,  -raddiv2,  -raddiv2,  -raddiv2];

Face5 = [htdiv2,  -htdiv2,  -htdiv2,  htdiv2;...
    raddiv2, raddiv2,  raddiv2, raddiv2;...
    -raddiv2,  -raddiv2, raddiv2, raddiv2];

Face6 = [htdiv2, -htdiv2, -htdiv2, htdiv2;...
    -raddiv2, -raddiv2, -raddiv2, -raddiv2;...
    -raddiv2, -raddiv2, raddiv2, raddiv2];

for pp=1:4
    Face1(:,pp) = RotMat*Face1(:,pp);
    Face2(:,pp) = RotMat*Face2(:,pp);
    Face3(:,pp) = RotMat*Face3(:,pp);
    Face4(:,pp) = RotMat*Face4(:,pp);
    Face5(:,pp) = RotMat*Face5(:,pp);
    Face6(:,pp) = RotMat*Face6(:,pp);
end

XData = [Face1(1,:); Face2(1,:); Face3(1,:);...
    Face4(1,:); Face5(1,:); Face6(1,:)]';

YData = [Face1(2,:); Face2(2,:); Face3(2,:);...
    Face4(2,:); Face5(2,:); Face6(2,:)]';

ZData = [Face1(3,:); Face2(3,:); Face3(3,:);...
    Face4(3,:); Face5(3,:); Face6(3,:)]';

XData = r0(1) + XData;
YData = r0(2) + YData;
ZData = r0(3) + ZData;


end