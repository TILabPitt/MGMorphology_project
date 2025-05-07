

%% Enter file location of .mat results after processing through MGMorph. This code will convert to .xlsx file  
listdir = dir('*mat*');

for i = 1:length(listdir)
    filename = fullfile(listdir(i).folder,listdir(i).name);
    [temps,temps2] = resultstocsv(filename)
    newname = erase(listdir(i).name,'.mat');
    writetable(temps,[newname,'.xlsx'])
end