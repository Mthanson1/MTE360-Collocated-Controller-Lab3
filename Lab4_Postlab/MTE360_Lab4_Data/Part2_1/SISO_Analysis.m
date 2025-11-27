clc; close all;

%% System parameters
mT = 8.5E-4;
bT = 0.0033;
dT = 0.44;
k  = 0.07;
C  = 0;
m1 = mT/(1+0.54);
m2 = mT - m1;

b1 = bT/(1+0.1);
b2 = bT - b1;

d1 = dT/(1+0.1);
d2 = dT - d1;

%% Transfer function

a1 = (m1*b2 + m2*b1 + (m1+m2)*C) / (m1*m2);
a2 = (b1*b2 + (b1+b2)*C + (m1+m2)*k)/(m1*m2);
a3 = (b1+b2)*k/(m1*m2);

s  = tf('s');
G2 = 1/(m1*m2)*(C*s+k) / (s*(s^3+a1*s^2+a2*s+a3));
[Wn,Zeta] = damp(G2);
wn = Wn(3)
zeta = Zeta(3)

sisotool(G2)