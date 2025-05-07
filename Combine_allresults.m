

%%
data917 = readtable('20240917_Tser053microglia_rest-im_Probabilities.xlsx')
str = "PreInjection"
data917.Injection = repmat([str],[size(data917,1),1])
str1 = "APOE3"
data917.APOE = repmat([str1],[size(data917,1),1])
str2 = "Male"
data917.Sex = repmat([str2],[size(data917,1),1])
str3 = 'Mouse0917'
data917.FileName = repmat([str3],[size(data917,1),1]);


%%
data924 = readtable('20240924microglia_preinj-im_Probabilities.xlsx')
str = "PreInjection"
data924.Injection = repmat([str],[size(data924,1),1])
str1 = "APOE4"
data924.APOE = repmat([str1],[size(data924,1),1])
str2 = "Female"
data924.Sex = repmat([str2],[size(data924,1),1])
str3 = 'Mouse0924'
data924.FileName = repmat([str3],[size(data924,1),1]);

%% This file will append all excel results files to one excel sheet to be read into R.
% If results are first time, manually enter in APOE, sex and filename
% status for each xlsx file. 

data926 = readtable("20240926microglia_preinj-im_Probabilities.xlsx")
str = "PreInjection"
data926.Injection = repmat([str],[size(data926,1),1])
str1 = "APOE3"
data926.APOE = repmat([str1],[size(data926,1),1])
str2 = "Male"
data926.Sex = repmat([str2],[size(data926,1),1])
str3 = 'Mouse0926'
data926.FileName = repmat([str3],[size(data926,1),1]);
%%%
existing_data = readtable('existing_results.xlsx') % add to existing results

new_table = [data917;data924;data926;data1003; data1031;data1029;data1024;data1022;data1017;data1105;data1107;data1119;data1203;data1205;data1121]
writetable(new_table,'20241213_results_micelabel.xlsx')