%% ------------------------------------------------------------------------
% 2) Identified model parameters (note lower-case c)
% -------------------------------------------------------------------------

b1 = 0.0257; %[Vs/mm]
b2 = 0.0026; %[Vs/mm]

d1 = 0.5431;  %[V]       
d2 = 0.0543; %[V] 

m1 = 0.00031231;   %[Vs^2/mm]   
m2 = 0.00016865; %[Vs^2/mm]

k  = 0.025; %[V/mm]
C  = 0.000; %[Vs/mm]

a1= (m1*b2 + m2*b1) / (m1*m2);
a2= (b1*b2 + (m1+m2)*k) / (m1*m2);
a3= ((b1+b2) * k) / (m1*m2);

num = [C k];
num = (1/(m1*m2)) .* num;
den = [1 a1 a2 a3 0];

G = tf(num,den);
disp("X2 Transfer Function: "); G

sisotool(G);


%% Notch/Lead/Lag Filter

G_notch = tf([0.066^2 0.093 1],[0.066^2 0.132 1]);
G_lead = tf([1 7.8824697],[1 79.28987]);

omega_c

K_p = 1/abs(freqresp(G*G_notch*G_lead^2,omega_c));

