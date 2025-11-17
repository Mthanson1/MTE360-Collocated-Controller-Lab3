clc; clear; close all;

baseDir = fileparts(mfilename('fullpath'));

% List all voltage subfolders (e.g., '1_0V', '1_5V')
voltageFolders = dir(fullfile(baseDir,'Lab3_Group6_Data/2-5/', '*V'));
voltageFolders = voltageFolders([voltageFolders.isdir]);  % keep only folders

u0 = zeros(1,2);
Vss = zeros(1,2);

alpha = 0.54; % =m2/m1
beta = 0.1; % =b2/b1
gamma = 0.1; % =d2/d1

% Filter Velocity Data to Find Vss
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

    % Divide by number of files to get average
    for k = 1:length(fields)
        field = fields{k};
        if ~isempty(meanData.(field))
            meanData.(field) = meanData.(field) / nFiles;
        end
    end

    voltageStr = strrep(vFolderName, '_', '.');
    u0(i) = str2double(voltageStr(1:end-1));

    fs = 1000;

    % crop first 20% of data to remove transients
    initial_cut = 0.25 * 8; % =1.6s
    samples_cut = initial_cut * fs;
    v1_trim = meanData.v1_filt(samples_cut+1:end);
    v2_trim = meanData.v2_filt(samples_cut+1:end);

    % remove periodic 0.5s to isolate ss velocity
    original_win_start = 1.0;
    win_length = 0.5;
    period = 1.0;
    total_time = 8.0;

    % convert to trimmed-time frame
    trim_win_start = original_win_start - initial_cut;   
    trim_period = period;
    trim_win_length = win_length;

    mask = true(size(v1_trim));

    t = trim_win_start;
    while t < (total_time - initial_cut)
        if t >= 0
            startSample = floor(t * fs) + 1;
            endSample   = floor((t + trim_win_length) * fs);

            startSample = max(startSample, 1);
            endSample   = min(endSample, length(v1_trim));

            mask(startSample:endSample) = false;
        end
        t = t + trim_period;
    end

    % Compute d_t/b_t from Vss
    v1_clean = v1_trim(mask);
    v2_clean = v2_trim(mask);

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



end


A = [Vss', ones(2,1)];

x = A \ u0';

b_t = x(1);
d_t = x(2);


