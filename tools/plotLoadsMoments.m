function plotLoadsMoments(tspan,out)

m = length(out.ShearLoadY_N);
figure()
shear_y = plot(tspan(1:m),out.ShearLoadY_N); hold on;
shear_z = plot(tspan(1:m), out.ShearLoadZ_N);
tor_x = plot(tspan(1:m),out.Torsion_Nm); 
bend_y = plot(tspan(1:m),out.BendingMomY_Nm);
bend_z = plot(tspan(1:m),out.BendingMomZ_Nm);
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$Shear Loads \& Moments$','Interpreter','LaTeX');
ll = legend([shear_y, shear_z, tor_x, bend_y, bend_z],...
    {'$S_{\theta}(t)$','$S_z(t)$','$T_x(t)$','$B_\theta(t)$','$B_z(t)$'},'Location','Northeast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18,'Interpreter','LaTeX');
set(gca,'FontSize',16,'FontName','Times');
title('Chaser Shear Loads & Moments')
grid on