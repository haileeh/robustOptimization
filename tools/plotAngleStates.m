function plotAngleStates(tspan, out)

yaw1 = (out.yaw1); pitch1 = (out.pitch1); roll1 = (out.roll1);
yaw1d = out.yaw1d; pitch1d = out.pitch1d; roll1d = out.roll1d;
yaw1dd = out.yaw1dd; pitch1dd = out.pitch1dd; roll1dd = out.roll1dd;

figure;
subplot(3,1,1); hold on; grid on;
plot(tspan, yaw1, tspan, roll1, tspan, pitch1);
legend('yaw','roll','pitch','Location','bestoutside');
title('Angles');
xlabel('Time [s]'); ylabel('Angle [rad]');
xlim([tspan(1), tspan(end)]);

subplot(3,1,2); hold on; grid on;
plot(tspan, yaw1d, tspan, roll1d, tspan, pitch1d);
legend('d/dt yaw','d/dt roll','d/dt pitch','Location','bestoutside');
title('Anglular Velocity');
xlabel('Time [s]'); ylabel('Angular Velocity [rps]');
xlim([tspan(1), tspan(end)]);

subplot(3,1,3); hold on; grid on;
plot(tspan(1:length(yaw1dd)), yaw1dd, tspan(1:length(roll1dd)), roll1dd, tspan(1:length(pitch1dd)), pitch1dd);
legend('d^2/dt^2 yaw','d^2/dt^2 roll','d^2/dt^2 pitch','Location','bestoutside');
title('Anglular Acceleration');
xlabel('Time [s]'); ylabel('Angular Acceleration [rpsps]');
xlim([tspan(1), tspan(length(yaw1dd))]);