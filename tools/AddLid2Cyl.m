function [Lid1,Lid2]=AddLid2Cyl(x,y0,z0,CMToEdgeY_m,SatChoice)

theta = 0:0.05:2*pi;

y = y0+CMToEdgeY_m*cos(theta);
z = z0+CMToEdgeY_m*sin(theta);

z(end) = 0;

x1 = x(1);%0;

x2 = x(2);%h;

if strcmp(SatChoice, 'Chaser')
    col = 'b';
else
    col = [ 0.9100 0.4100 0.1700];
end

Lid1=patch(x1*ones(size(y)),y,z,col);

%set(gca,'NextPlot','Add');

Lid2=patch(x2*ones(size(y)),y,z,col);

%surf([y;y],[z;z],[x1*ones(size(y));x2*ones(size(y))],'parent',gca)

end