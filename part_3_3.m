clc; clear; close all;

%% ------------------------------------------------------------------------
% 0) Paths / model
% -------------------------------------------------------------------------
baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir,'Lab3_Group6_Data','2-5','2_0V');   % Section 2.5, U = ±2.0 V
mdl    = 'flexible_drive';

%% ------------------------------------------------------------------------
% 1) Load ONLY Section 2.5, U = ±2.0 V data (average over all trials)
% -------------------------------------------------------------------------
matFiles = dir(fullfile(dataDir,'**','*.mat'));   % recursive search
if isempty(matFiles)
    error('No MAT files found under %s', dataDir);
end

meanData = struct('t', [], 'u', [], ...
                  'v1', [], 'v1_filt', [], ...
                  'v2', [], 'v2_filt', []);
nFiles   = numel(matFiles);

for f = 1:nFiles
    S = load(fullfile(matFiles(f).folder, matFiles(f).name));

    fields = fieldnames(meanData);
    for k = 1:numel(fields)
        fn = fields{k};
        if isfield(S,fn)
            if isempty(meanData.(fn))
                meanData.(fn) = S.(fn);
            else
                meanData.(fn) = meanData.(fn) + S.(fn);
            end
        end
    end
end

% Convert sums to means
fields = fieldnames(meanData);
for k = 1:numel(fields)
    fn = fields{k};
    if ~isempty(meanData.(fn))
        meanData.(fn) = meanData.(fn)/nFiles;
    end
end

% Use filtered velocities if available
t      = meanData.t(:);
u_exp  = meanData.u(:);
if ~isempty(meanData.v1_filt), v1_exp = meanData.v1_filt(:); else, v1_exp = meanData.v1(:); end
if ~isempty(meanData.v2_filt), v2_exp = meanData.v2_filt(:); else, v2_exp = meanData.v2(:); end

%% ------------------------------------------------------------------------
% 2) Identified model parameters (note lower-case c)
% -------------------------------------------------------------------------

b1      = 0.0265;      
b2      = 0.0027;        

dcoul_1 = 0.5051;       
dcoul_2 = 0.0505;

m1      = 0.00015955;   
m2      = 0.00008515;

k       = 0.065;        
c       = 0.000;        


%% ------------------------------------------------------------------------
% 3) Input signal for Simulink – Inport/From-Workspace named 'usim'
%    Set that block's data variable to 'usim' and format to Timeseries.
% -------------------------------------------------------------------------
usim = timeseries(u_exp, t);
usim.TimeInfo.Units = 'seconds';
assignin('base','usim',usim);   % ensure base workspace has it

%% ------------------------------------------------------------------------
% 4) Run simulation
% -------------------------------------------------------------------------
stopTime = num2str(t(end));
load_system(mdl);

simOut = sim(mdl, ...
             'StopTime', stopTime, ...
             'ReturnWorkspaceOutputs','on');

% With your To Workspace settings (Save format: Array),
% v1sim and v2sim are value-only vectors; time is 'tsim'.
t_sim  = simOut.tsim(:);
v1_sim = simOut.v1sim(:);
v2_sim = simOut.v2sim(:);

% Interpolate onto experimental time grid
v1_sim_i = interp1(t_sim, v1_sim, t, 'linear', 'extrap');
v2_sim_i = interp1(t_sim, v2_sim, t, 'linear', 'extrap');

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
ylabel('v_1 [rad/s]', 'Color','w');
title('Section 2.5, U = \pm2.0 V: v_1', 'Color','w');
legend({'Experiment','Simulation'}, 'TextColor','w', 'Location','Best', 'Box','off');
%ylim([-150 150])
darken(ax1);

% ---- v2 ----
ax2 = subplot(2,1,2); hold(ax2,'on');
plot(t, v2_exp,   'w',  'LineWidth', 1.4);
plot(t, v2_sim_i, 'r','LineWidth', 2.0);
grid on;
xlabel('Time [s]', 'Color','w');
ylabel('v_2 [rad/s]', 'Color','w');
title('Section 2.5, U = \pm2.0 V: v_2', 'Color','w');
legend({'Experiment','Simulation'}, 'TextColor','w', 'Location','Best', 'Box','off');
%ylim([-100 100])
darken(ax2);

%% ------------------------------------------------------------------------
% 6) Save figure
% -------------------------------------------------------------------------
outDir = fullfile(baseDir,'Section_3_Plots');
if ~exist(outDir,'dir'); mkdir(outDir); end
outPng = fullfile(outDir,'MTE360_2-5_1_5V_sim_overlay.png');
exportgraphics(fig, outPng, 'Resolution',300,'BackgroundColor','black');
