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

mT = 8.5E-4;
bT = 0.0033;
dT = 0.44;
k  = 0.07;
C = 0;
% m2/m1 = 0.54, b2/b1 = 0.1, d2/d1 = 0.1;
% m2 = 0.54m1, 0.54m1+m1 = mT, 
m1 = mT/(1+0.54);
m2 = mT - m1; %m2 = 0.54*m1;

b1 = bT/(1+0.1);
b2 = bT - b1; %b2 = 0.1*b1;

d1 = dT/(1+0.1);
d2 = dT - d1; %d2 = 0.1*d1;

q = 1000;  % derivative gain