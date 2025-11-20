clc; clear; close all;

baseDir = fileparts(mfilename('fullpath'));

%% List all voltage subfolders (e.g., '1_0V', '1_5V')
voltageFolders = dir(fullfile(baseDir,'Lab3_Group6_Data/2-5/', '*V'));
voltageFolders = voltageFolders([voltageFolders.isdir]);  % keep only folders
 
% keep ONLY '1_75V' and '2_0V'
keep = ismember({voltageFolders.name}, {'1_75V','2_0V'});
voltageFolders = voltageFolders(keep);

u0 = zeros(1,2);
Vss = zeros(1,2);
p1 = zeros(1); 
m1 = zeros(1,2); % Solve p1 for both u0=1.75V,2V
w_d = zeros(1,2);

alpha = 0.54; % =m2/m1
beta = 0.1; % =b2/b1
gamma = 0.1; % =d2/d1

%% Filter Velocity Data to Find Vss
for i = 1:length(voltageFolders)
    vFolderName = voltageFolders(i).name;
    vFolderPath = fullfile(baseDir, 'Lab3_Group6_Data/2-5/', vFolderName);
    
    % Find all .mat files in this voltage folder0.5
    matFiles = dir(fullfile(vFolderPath, '*.mat'));
    
    meanData = struct('t', [], 'u', [], 'v1', [], 'v1_filt', [], 'v2', [], 'v2_filt', []);
    nFiles = length(matFiles);


     for f = 1:nFiles
        data = load(fullfile(vFolderPath, matFiles(f).name));

        fields = fieldnames(meanData);
        for k = 1:length(fields)
            field = fields{k};
            if isfield(data, field)
                if isempty(meanData.(field))
                    meanData.(field) = data.(field);  % first file
                else
                    meanData.(field) = meanData.(field) + data.(field); % accumulate
                end
            end
        end
    end

    %% Divide by number of files to get average
    for k = 1:length(fields)
        field = fields{k};
        if ~isempty(meanData.(field))
            meanData.(field) = meanData.(field) / nFiles;
        end
    end

    voltageStr = strrep(vFolderName, '_', '.');
    u0(i) = str2double(voltageStr(1:end-1));

    %% Process data 
    fs = 1000;
    total_time = 8.0;

    % mask for b/d work
    mask_bd = makeMask(meanData.v1_filt, fs, total_time,...
                        0.25, ... % intial fraction cut (25%)
                        1.0,  ... % original win_start
                        0.5,  ... % win_length
                        1.0);     % period
    
    % mask for m1 work
    mask_m1 = makeMask(meanData.v1_filt, fs, total_time,...
                        0.25, ... % intial fraction cut (25%)
                        1.5,  ... % original win_start
                        0.5,  ... % win_length
                        1.0);     % period

    % crop first 25% of data to remove transients
    initial_cut = 0.25 * 8; % =1.6s
    samples_cut = initial_cut * fs;
    v1_trim = meanData.v1_filt(samples_cut+1:end);
    v2_trim = meanData.v2_filt(samples_cut+1:end); 

    %% Compute d_t/b_t from Vss
    v1_clean = v1_trim(mask_bd);
    v2_clean = v2_trim(mask_bd);

    den = max([abs(v1_clean), abs(v2_clean)], [], 2);
    percent_diff = abs(v1_clean - v2_clean) ./ den;
    percent_diff(isnan(percent_diff)) = 0; %avoid divide by zero errors

    keep_mask = percent_diff <= 0.15;
    
    v1_band = v1_clean(keep_mask);
    v2_band = v2_clean(keep_mask);

    mean_v1 = mean(abs(v1_band));
    mean_v2 = mean(abs(v2_band));
    
    avg_v = (mean_v1+mean_v2)/2;
    avg_v(isnan(avg_v))=0; % set to zero if no data found

    Vss(i) = avg_v;

    %% Compute m1

    peak_period = 500; %peaks repeat every 500 samples in cleaned data

    m1_data = abs(v1_trim(mask_m1));
    m1_data = m1_data(1:end-1); % convert from 3001 -> 3000 entities
    m1_sum = reshape(m1_data, peak_period, []);
    m1_avg = mean(m1_sum,2);

    % trim till start of rise
    [~, minidx] = min(m1_avg);
    m1_avg = m1_avg(minidx:end);


    m1_ss = mean_v1; % switch to using actual steady state velocity
    rise_targ = [0.1, 0.9]; % 10-90% = rise time
    rise_to_tau = log(9); % relation between rise time and time constant
    m1_targ = rise_targ .* m1_ss;

    t_start_m1 = 0;
    t_time_c_m1 = 0;

    for j = 2:length(m1_avg)
    
        if abs(m1_avg(j)-m1_targ(1)) <= 0.1 * m1_targ(2) && t_start_m1 == 0
            t_start_m1 = j;
        end

        if abs(m1_avg(j) - m1_targ(2)) <= 0.1 * m1_targ(2)
            t_time_c_m1 = j;
            break;
        end
    end
    tau_m1 = (t_time_c_m1-t_start_m1)/rise_to_tau/1000; % in s

    K_v = mean_v1/u0(i); %in units of mm/(Vs)

    m1(i) = tau_m1/K_v;

    %% Compute find w_d

    %find all peaks
    [~, t_max] = findpeaks(v2_trim, meanData.t(samples_cut+1:end), 'MinPeakProminence', 10);
    [~, t_min] = findpeaks(-v2_trim, meanData.t(samples_cut+1:end), 'MinPeakProminence', 10);
    t_all = [t_max(:); t_min(:)];
    t_all = sort(t_all);

    dt = diff(t_all);

    % sort for f > 0.5Hz
    dt_min = min(dt);
    tol = max(dt)/2;
    valid_idx = (dt <= tol);

    dt_filt = dt(valid_idx);

    T_avg = mean(dt_filt);
    w_d(i) = 2*pi()/T_avg;
end

%% extract parameter Characteristics
A = [Vss', ones(2,1)];

x = A \ u0';

b_t = x(1);
d_t = x(2);

%b_t = (1+beta)b1 = (1/beta + 1)b2

b(1) = b_t/(1+beta);
b(2) = b_t/(1/beta + 1);

d(1) = d_t/(1+gamma);
d(2) = d_t/(1/gamma + 1);

m1_avg = mean(m1);


m_t = (1+alpha)*m1_avg;

m(1) = m1_avg;
m(2) = m_t - m1_avg;

% find real pole and mass constants
p1 = b_t/m_t;

% coverge to damping ratio
w_d_avg = mean(w_d);

% Iteration settings
zeta = 0.9; %initial guess
tol = 1e-5; % convergance tolerance
max_iter = 20; %max number of iterations

fprintf("Starting iteration...\n");
fprintf("Initial ζ guess = %.4f\n", zeta);

for iter = 1:max_iter
    w_n = w_d_avg / sqrt(1-zeta^2);

    k_approx = m(1) * m(2) * w_n^2 / m_t;
    
    s = tf('s');

    a_1= (m(1)*b(2) + m(2)*b(1)) / (m(1)*m(2));
    a_2= (b(1)*b(2) + m_t*k_approx) / (m(1)*m(2));
    a_3= (b_t * k_approx) / (m(1)*m(2));

    den = s*(s^3 + a_1 * s^2 + a_2 * s + a_3);
    sys = 1/den;

    [d_wn, d_z, d_p] = damp(sys); % each row: [wn  zeta  p]

    % find the complex pole pair (zeta < 1 and imaginary part ~= 0)
    idx = find(d_z < 1 & imag(d_p) ~= 0, 1, 'first');

    if isempty(idx)
        error('No complex poles found. Check parameters!');
    end

    new_zeta = d_z(idx);

    fprintf("Iter %d: ζ_old = %.5f → ζ_new = %.5f\n", iter, zeta, new_zeta);

    if abs(new_zeta - zeta) < tol
        fprintf("\nConverged after %d iterations!\n", iter);
        break;
    end

    zeta = new_zeta;
end

fprintf("\nFinal damping ratio ζ = %.6f\n", zeta);
fprintf("Final stiffness k = %.3f N/m\n", k);
fprintf("Final ω_n = %.3f rad/s\n", w_n);

%% Function that removes periodic mask from data:
function mask = makeMask(v1, fs, total_time, initial_cut_frac, ...
                         original_win_start, win_length, period)

    % --- Crop first X% of data (initial transient removal) ---
    initial_cut = initial_cut_frac * total_time;      % seconds
    samples_cut = round(initial_cut * fs);

    v1_trim = v1(samples_cut+1:end);   % only used for sizing mask
    N = length(v1_trim);

    % --- Convert timing values into trimmed reference frame ---
    trim_win_start  = original_win_start - initial_cut;
    trim_win_length = win_length;
    trim_period     = period;

    % --- Initialize mask (default = keep data) ---
    mask = true(N,1);

    % --- Apply periodic window removal ---
    t = trim_win_start;

    while t < (total_time - initial_cut)
        if t >= 0
            startSample = floor(t * fs) + 1;
            endSample   = floor((t + trim_win_length) * fs);
            % Bounds check
            startSample = max(startSample, 1);
            endSample   = min(endSample, N);

            mask(startSample:endSample) = false;
        end
        t = t + trim_period;
    end
end