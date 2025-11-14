% m function to fit a sine wave on a noisy signal
%
% function [Amplitude , Phase] = fit_sine_wave(t,y,w)
%
% Least Squares parameter estimation is used to find
% a function in the form: 
%
% y_fit = Amplitude * sin( w*t + Phase )
% 
% which matches the data set t[time] and y[signal]
% as closely as possible 
% 
% inputs:
% t [sec]    : time array
% y [   ]    : noisy signal that a sine wave has to be fit
% w [rad/sec]: frequency of sine wave to be fit
% 
% outputs:
% Amplitude [ ]  : amplitude of sine wave
% Phase     [rad]: phase of sine wave

function [Amplitude , Phase] = fit_sine_wave(t,y,w)

Phi = [sin(w*t) cos(w*t)];                  % construct regressor matrix
Y = y;                                      % noisy output signal
theta = pinv(Phi)*Y;                        % find parameters [A*cos(pha) , A*sin(pha)]

Amplitude = sqrt(theta(1)^2 + theta(2)^2);  % amplitude
Phase = atan2(theta(2),theta(1));           % phase


