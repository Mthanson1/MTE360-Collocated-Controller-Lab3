clc;clear;close all
%
data=load('input_trajectory.txt');
time=data(:,1);
s=data(:,2);
Ts = 0.001;
%
clear data
%
% YOU NEED TO PUT YOUR OWN VALUES FOR THE MODEL PARAMETERS


k  = 0.025;
C = 0;
% m2/m1 = 0.54, b2/b1 = 0.1, d2/d1 = 0.1;
% m2 = 0.54m1, 0.54m1+m1 = mT, 
m1 = 0.00031231;
m2 = 0.00016865; %m2 = 0.54*m1;

b1 = 0.0257;
b2 = 0.0026; %b2 = 0.1*b1;

d1 = 0.5431;
d2 = 0.0543; %d2 = 0.1*d1;


mT =m1+m2;
bT = b1+b2;
dT = d1+d2;

q = 1000;  % derivative gain