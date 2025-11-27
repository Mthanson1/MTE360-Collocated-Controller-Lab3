
save('SimResults_Part3_2','xr','x1sim','x2sim','usim');

figure(1)
subplot(4,1,1)
plot(x1sim.time,x1sim.signals.values)
xlabel('time [sec]')
ylabel('position [mm]')
legend('x1')

subplot(4,1,2)
plot(xr.time,xr.signals.values,x2sim.time,x2sim.signals.values)
xlabel('time [sec]')
ylabel('position [mm]')
legend('xr','x2')


subplot(4,1,3)
plot(xr.time,xr.signals.values-x2sim.signals.values)
xlabel('time [sec]')
ylabel('position [mm]')
legend('e')

subplot(4,1,4)
plot(usim.time,usim.signals.values)
ylabel('control signal [V]')
xlabel('time [sec]')
legend('u')