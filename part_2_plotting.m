clc; clear; close all;
% --- CONFIG ----
baseDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(baseDir, 'Section_2_Plots');
cases   = {'2-1','2-2','2-3','2-4'};     % will run all and make five figures
% ----------------

for k = 1:numel(cases)
    caseName = cases{k};
    caseFolder = fullfile(baseDir, 'Lab3_Group6_Data', caseName);

    % Find .mat file inside folder
    files = dir(fullfile(caseFolder, '*.mat'));
    if isempty(files)
        error('No MAT file found in folder %s', caseFolder);
    end

    fname = files(1).name;
    filepath = fullfile(caseFolder, fname);

    % Extract last X-XXX pattern (e.g., 0-270)
    tokens = regexp(fname, '(\d+-\d+)(?=\.mat)', 'tokens');
    
    if ~isempty(tokens)
        raw = tokens{end}{1};         % e.g., '0-270'
        parts = split(raw, '-');      % {'0', '270'}
        major = parts{1};
        minor = parts{2};
        value = str2double([major '.' minor]);
    else
        % Default value
        value = 0.05;
    end

    % Load .mat
    data = load(fname);
    
    fig = plot_data_p_control(data,caseName,value);

    file_name = sprintf('MTE360_%s_overlay.png',caseName);
    outPng = fullfile(outputDir, file_name);

    exportgraphics(fig,outPng,'Resolution',300,'BackgroundColor','black');
end

%% List all voltage subfolders (e.g., '1_0V', '1_5V')
voltageFolders = dir(fullfile(baseDir,'Lab3_Group6_Data/2-5/', '*V'));
voltageFolders = voltageFolders([voltageFolders.isdir]);  % keep only folders

for v = 1:length(voltageFolders)
    vFolderName = voltageFolders(v).name;
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
    val = ['U=' voltageStr];

    fig = plot_velocity_response(meanData,'2-5',val);

    file_name = sprintf('MTE360_2-5_%s_overlay.png',vFolderName);
    outPng = fullfile(outputDir, file_name);

    exportgraphics(fig,outPng,'Resolution',300,'BackgroundColor','black');
end

%% Plot requested frequency responses
freq = {'0_5', '4_0', '8_0', '16_0'}; % Generate a figure for each frequency listed here in format XX.X = XX_X

for j = 1:length(freq)
    freqFileName = sprintf('2_6-%sHz.mat',freq{j});
    fFilePath = fullfile(baseDir, 'Lab3_Group6_Data/2-6/', freqFileName);

    freqStr = ['f=' strrep(freq{j}, '_', '.') 'Hz'];

    data = load(fFilePath);
    
    fig = plot_velocity_response(data,'2-6', freqStr);

    file_name = sprintf('MTE360_2-6_freq-%s.png',num2str(j));
    outPng = fullfile(outputDir, file_name);

    exportgraphics(fig,outPng,'Resolution',300,'BackgroundColor','black');
end