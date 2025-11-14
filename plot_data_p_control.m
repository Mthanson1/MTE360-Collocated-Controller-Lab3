clc;
%
t_exp=ScopeData_x1.time;
xr=ScopeData_xr.signals.values;
x1=ScopeData_x1.signals.values;
x2=ScopeData_x2.signals.values;
u=ScopeData_u.signals.values;
%%
 
figure(2)
f(1)=subplot(4,1,1);
plot(t_exp,xr,t_exp,x1)
ylabel('Position x1 [mm]')
grid on

f(2)=subplot(4,1,2);
plot(t_exp,xr , t_exp, x2)
ylabel('Position x2 [mm]')
grid on

f(3)=subplot(4,1,3);
plot(t_exp,xr-x2);
ylabel('Error [mm]')
grid on;

f(3)=subplot(4,1,4);
plot(t_exp,u)
ylabel('Input [V]')
xlabel('Time [sec]')
grid on
linkaxes(f,'x');

%% Modify the .mat file name to be representitive of the data you are saving
save('exp_data','xr','x1','x2','u','t_exp');