% Read dataset
dataset = readtable('/Users/noahschweitzer/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Microglia_project/experiments_avg20250724.csv');
%dataset = readtable('/Users/nes93/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Microglia_project/experiments_avg2.csv');



% change when you want to do young vs old

% finder = find(dataset.Age=="young");
% dataset = dataset(finder, :);
finder = find(dataset.Age=="old");
dataset = dataset(finder, :);


% Initialize new columns
dataset.AdjustedTime = NaN(height(dataset), 1);
dataset.PercentChange = NaN(height(dataset), 1);                % For mean_FullCellComplexity
dataset.PercentChange_t_index = NaN(height(dataset), 1);       % For mean_t_index
dataset.PercentChange_logBranch = NaN(height(dataset), 1);     % For mean_log_MaxBranchLength

% Get unique FileNames
fileNames = unique(dataset.FileName);




% Adjust time per FileName based on earliest PostInjection
for i = 1:length(fileNames)
    currentFile = fileNames{i};
    idx = strcmp(dataset.FileName, currentFile);
    fileData = dataset(idx, :);

    postTimes = fileData.Time(strcmp(fileData.Injection, 'PostInjection'));
    if isempty(postTimes) || all(isnan(postTimes))
        warning('No PostInjection time for file: %s', currentFile);
        continue;
    end

    minPostInjectionTime = min(postTimes);
    dataset.AdjustedTime(idx) = fileData.Time - minPostInjectionTime;
end

% Filter baseline data
baselineData = dataset(strcmp(dataset.Injection, 'PreInjection'), :);

% Compute percent change for each metric
for i = 1:length(fileNames)
    currentFile = fileNames{i};
    idx = strcmp(dataset.FileName, currentFile);
    fileData = dataset(idx, :);

    % Baselines
    baselineFCC = mean(baselineData.mean_FullCellComplexity(strcmp(baselineData.FileName, currentFile)), 'omitnan');
    baselineTI  = mean(baselineData.mean_t_index(strcmp(baselineData.FileName, currentFile)), 'omitnan');
    baselineLB  = mean(baselineData.mean_log_MaxBranchLength(strcmp(baselineData.FileName, currentFile)), 'omitnan');

    % Percent change
    dataset.PercentChange(idx) = ((fileData.mean_FullCellComplexity - baselineFCC) / baselineFCC) * 100;
    dataset.PercentChange_t_index(idx) = ((fileData.mean_t_index - baselineTI) / baselineTI) * 100;
    dataset.PercentChange_logBranch(idx) = ((fileData.mean_log_MaxBranchLength - baselineLB) / baselineLB) * 100;
end

% Grouping by APOE_Sex
groups = unique(dataset.APOE_Sex);
numGroups = length(groups);

% Common plot settings
xLimits = [0 80];

% ---- Plot 1: FullCellComplexity ----
figure('Name', 'Percent Change: FullCellComplexity');
for g = 1:numGroups
    currentGroup = groups{g};
    groupIdx = strcmp(dataset.APOE_Sex, currentGroup);
    groupData = dataset(groupIdx, :);
    groupFiles = unique(groupData.FileName);

    subplot(numGroups, 1, g);
    hold on;

    for i = 1:length(groupFiles)
        currentFile = groupFiles{i};
        fileIdx = strcmp(groupData.FileName, currentFile);
        fileData = groupData(fileIdx, :);

        if all(isnan(fileData.PercentChange))
            continue;
        end

        smoothedFCC = smooth(fileData.AdjustedTime, fileData.PercentChange, 0.2, 'lowess');
        plot(fileData.AdjustedTime, smoothedFCC, 'LineWidth', 1.5, 'DisplayName', currentFile);
    end

    ylabel('Percent Change (%)');
    title(['FullCellComplexity - APOE-Sex: ' currentGroup]);
    xlim(xLimits);
    if g == numGroups
        xlabel('Adjusted Time (Relative to Injection)');
    end
    legend('Interpreter', 'none');
    grid on;
    hold off;
end

% ---- Plot 2: t_index ----
figure('Name', 'Percent Change: t_index');
for g = 1:numGroups
    currentGroup = groups{g};
    groupIdx = strcmp(dataset.APOE_Sex, currentGroup);
    groupData = dataset(groupIdx, :);
    groupFiles = unique(groupData.FileName);

    subplot(numGroups, 1, g);
    hold on;

    for i = 1:length(groupFiles)
        currentFile = groupFiles{i};
        fileIdx = strcmp(groupData.FileName, currentFile);
        fileData = groupData(fileIdx, :);

        if all(isnan(fileData.PercentChange_t_index))
            continue;
        end

        smoothedTI = smooth(fileData.AdjustedTime, fileData.PercentChange_t_index, 0.2, 'lowess');
        plot(fileData.AdjustedTime, smoothedTI, 'LineWidth', 1.5, 'DisplayName', currentFile);
    end

    ylabel('Percent Change (%)');
    title(['t\_index - APOE-Sex: ' currentGroup]);
    xlim(xLimits);
    if g == numGroups
        xlabel('Adjusted Time (Relative to Injection)');
    end
    legend('Interpreter', 'none');
    grid on;
    hold off;
end

% ---- Plot 3: log_MaxBranchLength ----
figure('Name', 'Percent Change: log(MaxBranchLength)');
for g = 1:numGroups
    currentGroup = groups{g};
    groupIdx = strcmp(dataset.APOE_Sex, currentGroup);
    groupData = dataset(groupIdx, :);
    groupFiles = unique(groupData.FileName);

    subplot(numGroups, 1, g);
    hold on;

    for i = 1:length(groupFiles)
        currentFile = groupFiles{i};
        fileIdx = strcmp(groupData.FileName, currentFile);
        fileData = groupData(fileIdx, :);

        if all(isnan(fileData.PercentChange_logBranch))
            continue;
        end

        smoothedLB = smooth(fileData.AdjustedTime, fileData.PercentChange_logBranch, 0.2, 'lowess');
        plot(fileData.AdjustedTime, smoothedLB, 'LineWidth', 1.5, 'DisplayName', currentFile);
    end

    ylabel('Percent Change (%)');
    title(['log(MaxBranchLength) - APOE-Sex: ' currentGroup]);
    xlim(xLimits);
    if g == numGroups
        xlabel('Adjusted Time (Relative to Injection)');
    end
    legend('Interpreter', 'none');
    grid on;
    hold off;
end