%% Initialization and Parameters
clc; clear; close all;

% =========================================================================
% USER SELECTION: Choose Experimental Data Set
% =========================================================================
voltage_choice = '8v'; % Set this to '4v' or '8v'

if strcmpi(voltage_choice, '4v')
    exp_filename = '4v_sat.mat';
    disp(['Using ' exp_filename]);
elseif strcmpi(voltage_choice, '8v')
    exp_filename = '8v_sat.mat';
    disp(['Using ' exp_filename]);
else
    error('Invalid selection.');
end

% Create output directory
outDir = 'Section_3_Plots';
if ~exist(outDir, 'dir'), mkdir(outDir); end

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



Ts = 0.001; 

%% SECTION 3.2: Notch and Double Lead
disp('Processing Section 3.2...');
model_3_2 = fullfile('Part2_2', 'Simulation', 'flexibleDrive_NotchLead.slx');
exp_file_3_2 = fullfile('Part2_2', 'Experiment', exp_filename);

if isfile(model_3_2) && isfile(exp_file_3_2)
    % 1. Load Experiment FIRST to get the correct Reference
    expData = load(exp_file_3_2);
    
    % 2. Inject Experimental Reference into Sim Inputs (time, s)
    % This fixes the "1 peak vs 2 peak" issue by using the real experiment command
    [t_ref, x_ref] = extract_exp_data(expData, 'xrexp'); 
    
    if ~isempty(t_ref)
        time = t_ref; 
        s = x_ref;
    else
        % Fallback if xrexp is missing (try loading text file)
        warning('Reference not found in MAT file, loading default.');
        d = load('input_trajectory.txt'); time = d(:,1); s = d(:,2);
    end
    
    % 3. Run Simulation (Now uses the updated 'time' and 's')
    sim(model_3_2);
    
    % 4. Capture Results
    if exist('x1sim', 'var'), s_x1 = x1sim; else, s_x1 = []; end
    if exist('x2sim', 'var'), s_x2 = x2sim; else, s_x2 = []; end
    if exist('usim', 'var'),  s_u  = usim;  else, s_u  = []; end
    
    % 5. Plot
    plot_name = ['Section_3.2_' voltage_choice '_NotchLead'];
    save_path = fullfile(outDir, [plot_name '.png']);
    plot_results_dark(s_x1, s_x2, s_u, expData, time, s, ...
        ['Section 3.2 (' voltage_choice '): Notch & Double Lead'], save_path);
else
    warning('Skipping 3.2: Model or Data file not found.');
end

%% SECTION 3.3: Notch, Double Lead and Lag
disp('Processing Section 3.3...');
model_3_3 = fullfile('Part2_3', 'Simulation', 'flexibleDrive_NotchLeadLag.slx');
exp_file_3_3 = fullfile('Part2_3', 'Experiment', exp_filename);

if isfile(model_3_3) && isfile(exp_file_3_3)
    expData = load(exp_file_3_3);
    
    % Inject Reference
    [t_ref, x_ref] = extract_exp_data(expData, 'xrexp');
    if ~isempty(t_ref), time = t_ref; s = x_ref; end
    
    sim(model_3_3);
    
    if exist('x1sim', 'var'), s_x1 = x1sim; else, s_x1 = []; end
    if exist('x2sim', 'var'), s_x2 = x2sim; else, s_x2 = []; end
    if exist('usim', 'var'),  s_u  = usim;  else, s_u  = []; end
    
    plot_name = ['Section_3.3_' voltage_choice '_NotchLeadLag'];
    save_path = fullfile(outDir, [plot_name '.png']);
    plot_results_dark(s_x1, s_x2, s_u, expData, time, s, ...
        ['Section 3.3 (' voltage_choice '): Notch, Lead & Lag'], save_path);
end

%% SECTION 3.4: Feedforward
disp('Processing Section 3.4...');
model_3_4 = fullfile('Part2_4', 'Simulation', 'flexibleDrive_compensated_feedforward.slx');
exp_file_3_4 = fullfile('Part2_4', 'Experiment', exp_filename);

if isfile(model_3_4) && isfile(exp_file_3_4)
    expData = load(exp_file_3_4);
    
    % Inject Reference
    [t_ref, x_ref] = extract_exp_data(expData, 'xrexp');
    if ~isempty(t_ref), time = t_ref; s = x_ref; end
    
    sim(model_3_4);
    
    if exist('x1sim', 'var'), s_x1 = x1sim; else, s_x1 = []; end
    if exist('x2sim', 'var'), s_x2 = x2sim; else, s_x2 = []; end
    if exist('usim', 'var'),  s_u  = usim;  else, s_u  = []; end
    
    plot_name = ['Section_3.4_' voltage_choice '_Feedforward'];
    save_path = fullfile(outDir, [plot_name '.png']);
    plot_results_dark(s_x1, s_x2, s_u, expData, time, s, ...
        ['Section 3.4 (' voltage_choice '): Feedforward'], save_path);
end


%% --- DARK MODE PLOTTING FUNCTION ---

function plot_results_dark(simX1, simX2, simU, expStruct, t_ref, s_ref, plotTitle, savePath)
    
    % --- 1. DATA PREP ---
    % Sim
    [t_sim, x1_s] = extract_struct_data(simX1);
    [~,     x2_s] = extract_struct_data(simX2);
    [~,     u_s]  = extract_struct_data(simU);
    
    % Use the injected reference directly for plotting the Sim Reference
    % (Since we forced 'time' and 's' to be the input, s_ref is the input)
    if ~isempty(t_sim)
        % Interpolate s_ref to match t_sim exactly for error calc
        xr_s = interp1(t_ref, s_ref, t_sim, 'linear', 'extrap');
    else
        xr_s = [];
    end

    % Exp
    [t_exp, x1_e] = extract_exp_data(expStruct, 'x1exp');
    [~,     x2_e] = extract_exp_data(expStruct, 'x2exp');
    [~,     u_e]  = extract_exp_data(expStruct, 'uexp');
    [~,     xr_e] = extract_exp_data(expStruct, 'xrexp');
    
    % Exp Fallbacks
    if isempty(u_e) && isfield(expStruct, 'Vm'), [~, u_e] = extract_exp_data(expStruct, 'Vm'); end
    % If experiment ref is missing, assume it's the same as what we injected
    if isempty(xr_e), xr_e = interp1(t_ref, s_ref, t_exp, 'linear', 'extrap'); end

    % --- 2. FIGURE SETUP (DARK MODE) ---
    fig = figure('Color','k', 'Position',[100 100 1200 800], 'Name', plotTitle);
    
    % Anonymous function for styling axes
    darken = @(ax) set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', ...
                       'GridColor', [0.35 0.35 0.35], ...
                       'MinorGridColor', [0.25 0.25 0.25], ...
                       'FontSize', 10, 'LineWidth', 1.1);

    sgtitle(plotTitle, 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');
    
    % --- 3. SUBPLOTS ---
    
    % Subplot 1: Cart 1 Position
    ax1 = subplot(4,1,1); hold(ax1, 'on');
    if ~isempty(t_exp), plot(t_exp, x1_e, 'w', 'LineWidth', 1.4); end
    if ~isempty(t_sim), plot(t_sim, x1_s, 'r', 'LineWidth', 2.0); end
    ylabel('x_1 [m]', 'Color', 'w');
    title('1) Cart 1 Position', 'Color', 'w');
    grid on; axis tight; 
    darken(ax1);
    if ~isempty(t_exp), legend({'Exp', 'Sim'}, 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none', 'Location', 'best'); end

    % Subplot 2: Commanded vs Cart 2
    ax2 = subplot(4,1,2); hold(ax2, 'on');
    % Plot Reference (Cyan Dotted)
    if ~isempty(t_ref), plot(t_ref, s_ref, 'c:', 'LineWidth', 1.5); end 
    if ~isempty(t_exp), plot(t_exp, x2_e, 'w', 'LineWidth', 1.4); end
    if ~isempty(t_sim), plot(t_sim, x2_s, 'r', 'LineWidth', 2.0); end
    ylabel('Pos [m]', 'Color', 'w');
    title('2) Reference vs Cart 2', 'Color', 'w');
    grid on; axis tight;
    darken(ax2);
    if ~isempty(t_sim), legend({'Ref', 'Exp', 'Sim'}, 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none', 'Location', 'best'); end

    % Subplot 3: Tracking Error
    ax3 = subplot(4,1,3); hold(ax3, 'on');
    if ~isempty(t_exp) && ~isempty(xr_e)
        e_exp = xr_e - x2_e;
        plot(t_exp, e_exp, 'w', 'LineWidth', 1.4);
    end
    if ~isempty(t_sim) && ~isempty(xr_s)
        e_sim = xr_s - x2_s;
        plot(t_sim, e_sim, 'r', 'LineWidth', 2.0);
    end
    ylabel('Error [m]', 'Color', 'w');
    title('3) Tracking Error (e = x_r - x_2)', 'Color', 'w');
    grid on; axis tight;
    darken(ax3);

    % Subplot 4: Control Signal
    ax4 = subplot(4,1,4); hold(ax4, 'on');
    if ~isempty(t_exp), plot(t_exp, u_e, 'w', 'LineWidth', 1.4); end
    if ~isempty(t_sim), plot(t_sim, u_s, 'r', 'LineWidth', 2.0); end
    ylabel('Voltage [V]', 'Color', 'w');
    xlabel('Time [s]', 'Color', 'w');
    title('4) Control Signal', 'Color', 'w');
    grid on; axis tight;
    darken(ax4);
    
    % --- 4. EXPORT FIGURE ---
    drawnow;
    exportgraphics(fig, savePath, 'Resolution', 300, 'BackgroundColor', 'black');
    disp(['Saved figure to: ' savePath]);
end

%% --- DATA EXTRACTION HELPERS ---

function [t, v] = extract_struct_data(s_in)
    t = []; v = [];
    if isempty(s_in), return; end
    try
        if isfield(s_in, 'time'), t = s_in.time; end
        if isfield(s_in, 'signals') && isfield(s_in.signals, 'values')
            v = s_in.signals.values;
        end
    catch
    end
end

function [t, v] = extract_exp_data(expS, fieldName)
    t = []; v = [];
    if isfield(expS, fieldName)
        raw = expS.(fieldName);
        if size(raw, 2) >= 2
            t = raw(:,1);
            v = raw(:,2);
        end
    end
end