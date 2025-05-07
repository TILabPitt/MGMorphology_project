
%% If results were separate into several files, we need to combine them to make one file again

pre = readtable("20250320_TSer3133_microglia_Probabilities.xlsx");
post1 = readtable("20250320_TSer3136_microglia_Probabilities_1.xlsx");
post2 = readtable("20250320_TSer3136_microglia_Probabilities_2.xlsx");
post3 = readtable("20250320_TSer3137_microglia_Probabilities_1.xlsx");
post4 = readtable("20250320_TSer3137_microglia_Probabilities_2.xlsx");
% Initialize the offset for the starting time of the first table
startTime = 0;
% Adjust the Time column for each table
pre.Time = pre.Time + startTime;
startTime = pre.Time(end);  % Update startTime based on the last entry in pre
str3 = "PreInjection";
pre.Injection = repmat([str3],[size(pre,1),1]);
post1.Time = post1.Time + startTime;
startTime = post1.Time(end);  % Update startTime based on the last entry in post1
str3 = "PostInjection";
post1.Injection = repmat([str3],[size(post1,1),1]);
post2.Time = post2.Time + startTime;
startTime = post2.Time(end);  % Update startTime based on the last entry in post2
str3 = "PostInjection";
post2.Injection = repmat([str3],[size(post2,1),1]);
post3.Time = post3.Time + startTime;
startTime = post3.Time(end);  % Update startTime based on the last entry in post2
str3 = "PostInjection";
post3.Injection = repmat([str3],[size(post3,1),1]);
post4.Time = post4.Time + startTime;
startTime = post4.Time(end);  % Update startTime based on the last entry in post2
post4.Injection = repmat([str3],[size(post4,1),1]);
combinedTable = [pre; post1; post2;post3;post4];
combinedTable
writetable(combinedTable, "20250320_combined.xlsx")