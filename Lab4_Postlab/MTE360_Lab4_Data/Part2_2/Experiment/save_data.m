
save('ExpResults_Part3_2','xr','x1exp','x2exp','uexp');

figure(1)
subplot(4,1,1)
plot(x1exp.time,x1exp.signals.values)
xlabel('time [sec]')
ylabel('position [mm]')
legend('x1')

subplot(4,1,2)
plot(xr.time,xr.signals.values,x2exp.time,x2exp.signals.values)
xlabel('time [sec]')
ylabel('position [mm]')
legend('xr','x2')


subplot(4,1,3)
plot(xr.time,xr.signals.values-x2exp.signals.values)
xlabel('time [sec]')
ylabel('position [mm]')
legend('e')

subplot(4,1,4)
plot(uexp.time,uexp.signals.values)
ylabel('control signal [V]')
xlabel('time [sec]')
legend('u')