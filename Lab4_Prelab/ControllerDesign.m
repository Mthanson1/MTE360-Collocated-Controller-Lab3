clc; clear; close all;

%% ------------------------------------------------------------------------
% 0) Paths / model
% -------------------------------------------------------------------------
baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir,'../Lab3_Group6_Data','2-4');   % Section 2.5, U = ±2.0 V
mdl    = {'NotchDLead','NotchDLeadLag','NotchDLeadLagFeedForward'};

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

%uncomment following line to do SISO work
%sisotool(G);


%% Notch/Lead/Lag Filters
wn  = 15.245;  %rad/s
zeta = sqrt(2);

a=7.8824697;
b=79.28987;

c=10;
d=0.1;

K_Val = {0.352, 32, 32};
q = 1000;

%G_notch = tf([1 2*zeta*wns wns^2],[1 2*wns wns^2]);
%G_lead = tf([1 c],[1 d]);
%omega_c = 25;
%K_p = 1/abs(freqresp(G*G_notch*G_lead^2,omega_c));

matFiles = dir(fullfile(dataDir,'**','*.mat'));   % recursive search
if isempty(matFiles)
    error('No MAT files found under %s', dataDir);
end

Data = load(fullfile(matFiles(1).folder, matFiles(1).name));

%Load Data from .mat file
t = Data.t_exp;
xr_exp = Data.xr;

%% ------------------------------------------------------------------------
%    Input signal for Simulink – Inport/From-Workspace named 'usim'
%    Set that block's data variable to 'usim' and format to Timeseries.
% -------------------------------------------------------------------------
usim = timeseries(xr_exp, t);
usim.TimeInfo.Units = 'seconds';
assignin('base','usim',usim);   % ensure base workspace has it

%% ------------------------------------------------------------------------
% 4) Run simulation
% -------------------------------------------------------------------------
stopTime = num2str(t(end));

for i = 1:1%numel(mdl)
    Kp=K_Val{i};
    load_system(mdl{i});

    simOut = sim(mdl{i}, ...
              'StopTime', stopTime, ...
               'ReturnWorkspaceOutputs','on');

    %% ------------------------------------------------------------------------
    % 5) Plot overlay: black background, white experiment, red simulation
    % -------------------------------------------------------------------------
    fig = figure('Color','k', 'Position',[100 100 1100 700]);
    
    darken = @(ax) set(ax, 'Color','k', 'XColor','w','YColor','w', ...
                        'GridColor',[0.35 0.35 0.35], ...
                        'MinorGridColor',[0.25 0.25 0.25]);

    % ---- v1 ----
    ax1 = subplot(2,1,1); hold(ax1,'on');
    plot(t, v1_exp,   'w',  'LineWidth', 1.4);
    plot(t, v1_sim_i, 'r','LineWidth', 2.0);
    grid on;
    xlabel('Time [s]', 'Color','w');
    ylabel('v_1 [mm/s]', 'Color','w');
    title('Section 3.3, U = \pm2.0 V: v_1', 'Color','w');
    legend({'Experiment','Simulation'}, 'TextColor','w', 'Location','Best', 'Box','off');
    %ylim([-150 150])
    darken(ax1);
    
    % ---- v2 ----
    ax2 = subplot(2,1,2); hold(ax2,'on');
    plot(t, v2_exp,   'w',  'LineWidth', 1.4);
    plot(t, v2_sim_i, 'r','LineWidth', 2.0);
    grid on;
    xlabel('Time [s]', 'Color','w');
    ylabel('v_2 [mm/s]', 'Color','w');
    title('Section 3.3, U = \pm2.0 V: v_2', 'Color','w');
    legend({'Experiment','Simulation'}, 'TextColor','w', 'Location','Best', 'Box','off');
    %ylim([-100 100])
    darken(ax2);
    
    %% ------------------------------------------------------------------------
    % 6) Save figure
    % -------------------------------------------------------------------------
    outDir = fullfile(baseDir,'Section_3_Plots');
    if ~exist(outDir,'dir'); mkdir(outDir); end
    outPng = fullfile(outDir,'MTE360_3-3_2_0V_sim_overlay.png');
    exportgraphics(fig, outPng, 'Resolution',300,'BackgroundColor','black');
end

