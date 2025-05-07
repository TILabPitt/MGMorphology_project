
function [temps,temps2] = resultstocsv(filename)

load(filename)
% Initialize arrays to store data
numendpts = [];
fullcellcomplexity = [];
maxBranchLength = [];
minBranchLength = [];
avgBranchLength = [];
numbranchpts = [];
% MaxBranchLength_touchingvess = [];
% AvgBranchLength_touchingvess = [];
% TotalBranches_nottouchingvess = [];
% TotalBranches_touchingvess = [];

% Get the number of values per time point
values_per_timepoint = cellfun(@numel, Microglia.numendpts);


% Loop through all  time points to extract data from each field
for i = 1:length(values_per_timepoint)
    numendpts = [numendpts, cell2mat(Microglia.numendpts{i})];
    fullcellcomplexity = [fullcellcomplexity, cell2mat(Microglia.FullCellComplexity{i})];
    maxBranchLength = [maxBranchLength, (Microglia.MaxBranchLength{i})];
    minBranchLength = [minBranchLength, (Microglia.MinBranchLength{i})];
    avgBranchLength = [avgBranchLength, (Microglia.AvgBranchLength{i})];
    numbranchpts = [numbranchpts, (Microglia.numbranchpts{i})];
%     MaxBranchLength_touchingvess = [MaxBranchLength_touchingvess, (Microglia.MaxBranchLength_touchingvess{i})];
%     AvgBranchLength_touchingvess = [AvgBranchLength_touchingvess, (Microglia.AvgBranchLength_touchingvess{i})];
%     TotalBranches_nottouchingvess = [TotalBranches_nottouchingvess, (Microglia.TotalBranches_nottouchingvess{i})];
%     TotalBranches_touchingvess = [TotalBranches_touchingvess, (Microglia.TotalBranches_touchingvess{i})];
end

% Transpose the arrays to match table requirements
fullcellcomplexity = fullcellcomplexity';
numendpts = numendpts';
maxBranchLength = maxBranchLength';
minBranchLength = minBranchLength';
avgBranchLength = avgBranchLength';
numbranchpts = numbranchpts';
% MaxBranchLength_touchingvess = MaxBranchLength_touchingvess';
% AvgBranchLength_touchingvess = AvgBranchLength_touchingvess';
% TotalBranches_nottouchingvess = TotalBranches_nottouchingvess';
% TotalBranches_touchingvess = TotalBranches_touchingvess';


% Determine the total length of entries
numlength = length(numendpts);


% Define the variable names and types, adding the new fields
varNames = ["TimePoint", "Time", 'Numendpts', 'FullCellComplexity', ...
            'MaxBranchLength', 'MinBranchLength', 'AvgBranchLength', 'NumBranchPts',...
            'MaxBranchLength_touchingvess','AvgBranchLength_touchingvess',...
            'TotalBranches_nottouchingvess','TotalBranches_touchingvess'];

varTypes = ["string", "double", "double", "double", ...
            "double", "double", "double", "double"...
            "double", "double", "double", "double"];

% Initialize the table
temps = table('Size', [numlength, length(varNames)], ...
              'VariableTypes', varTypes, 'VariableNames', varNames);

% Set starting and final indices for assigning the data
starting_index = 1;



% Populate the table with data across time points
for i = 1:length(values_per_timepoint)
    if i == 1
        final_index = values_per_timepoint(i);
    else
        final_index = starting_index + values_per_timepoint(i) - 1;
    end

    % Create time point label
    timepoint = ['time', num2str(i)];
    
    % Fill the table for the current time point
    temps.TimePoint(starting_index:final_index) = timepoint;
    temps.Time(starting_index:final_index) = i;
    temps.Numendpts(starting_index:final_index) = numendpts(starting_index:final_index);
    temps.FullCellComplexity(starting_index:final_index) = fullcellcomplexity(starting_index:final_index);
    temps.MaxBranchLength(starting_index:final_index) = maxBranchLength(starting_index:final_index);
    temps.MinBranchLength(starting_index:final_index) = minBranchLength(starting_index:final_index);
    temps.AvgBranchLength(starting_index:final_index) = avgBranchLength(starting_index:final_index);
    temps.NumBranchPts(starting_index:final_index) = numbranchpts(starting_index:final_index);
%     temps.MaxBranchLength_touchingvess(starting_index:final_index) = MaxBranchLength_touchingvess(starting_index:final_index);
%     temps.AvgBranchLength_touchingvess(starting_index:final_index) = AvgBranchLength_touchingvess(starting_index:final_index);
%     temps.TotalBranches_nottouchingvess(starting_index:final_index) = TotalBranches_nottouchingvess(starting_index:final_index);
%     temps.TotalBranches_touchingvess(starting_index:final_index) = TotalBranches_touchingvess(starting_index:final_index);

    % Update the starting index for the next iteration
    starting_index = final_index + 1;
end



varNames = ["TimePoint",'AvgDist','PercentMgVol'];
varTypes = ["double", "double","double"];
temps2 = table('Size',[length(Microglia.AvgDist),length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);


x = [1:length(Microglia.AvgDist)]';

temps2.AvgDist = cell2table(Microglia.AvgDist);
temps2.PercentMgVol = cell2table(Microglia.PercentMgVol);
temps2.TimePoint = array2table(x);

