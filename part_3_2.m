clc; clear; close all;

baseDir = fileparts(mfilename('fullpath'));

%% List all voltage subfolders (e.g., '1_0V', '1_5V')
voltageFolders = dir(fullfile(baseDir,'Lab3_Group6_Data/2-5/', '*V'));
voltageFolders = voltageFolders([voltageFolders.isdir]);  % keep only folders

u0 = zeros(1,2);
Vss = zeros(1,2);

alpha = 0.54; % =m2/m1
beta = 0.1; % =b2/b1
gamma = 0.1; % =d2/d1

%% Filter Velocity Data to Find Vss
for i = 1:length(voltageFolders)
    vFolderName = voltageFolders(i).name;
    vFolderPath = fullfile(baseDir, 'Lab3_Group6_Data/2-5/', vFolderName);
    
    % Find all .mat files in this voltage folder
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

    % mask for m2 work
    mask_m2 = makeMask(meanData.v1_filt, fs, total_time,...
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

    %% Compute m1 for give u0
    peak_period = 500; %peaks repeat every 500 samples in cleaned data

    m1 = abs(v1_trim(mask_m1));
    m1 = m1(1:end-1); % convert from 3001 -> 3000 entities
    m1_sum = reshape(m1, peak_period, []);
    m1_avg = mean(m1_sum,2);

    %m1_0 = m1_avg(1);
    m1_max = max(m1_avg);
    m1_targ = 0.63*m1_max;

    t_start_m1 = 0;
    t_time_c_m1 = 0;

    for j = 2:length(m1_avg)
        
        if m1_avg(j) > 0.01 && t_start_m1 == 0
            t_start_m1 = j;
        end

        if abs(m1_avg(j) - m1_targ) <= 0.1 * m1_targ
            t_time_c_m1 = j;
            break;
        end
    end
    
    tau_m1(i) = t_time_c_m1-t_start_m1;

    %% Compute m2 only for u_0=1.5V
    if u0(i) == 1.5
        start_trim = 100; % trim first 100ms to remove transients
        m2_pre = abs(v2_trim(mask_m2));
        m2_pre = m2_pre(start_trim:end-1);

        % make sample size a multiple of period:
        end_trim = mod(length(m2_pre),peak_period);
        m2 = m2_pre(1:end-end_trim);

        m2_sum = reshape(m2, peak_period, []);
        m2_avg = mean(m2_sum, 2);

        m2_max = max(m2_avg);
        m2_targ = 0.63 * m2_max;

        t_start_m2 = 0;
        t_time_c_m2 = 0;

        for j = 2:length(m2_avg)
            if m2_avg(j) > 0.01 && t_start_m2 == 0
                t_start_m2 = j;
            end

            if abs(m2_avg(j) - m2_targ) <= 0.1 * m2_targ
                t_time_c_m2 = j;
                break;
            end
        end
        tau_m2 = t_time_c_m2 - t_start_m2;
    end
end


A = [Vss', ones(2,1)];

x = A \ u0';

b_t = x(1);
d_t = x(2);


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