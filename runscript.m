
clear;clc;close all;
system('caffeinate &');


%% Enter your segmented volumes as a dir list here 
listdir_mg=dir('/Volumes/LaCie/MicrogliaProject_post2025/ilastik_microglia_processing/ilastik_processing03312025/*Prob*');



for i = 1:length(listdir_mg)
    
    %% Get files
    fullfilename = fullfile(listdir_mg(i).folder,listdir_mg(i).name);
    filename = listdir_mg(i).name;
%     vess_dir = dir(['Segmentations_20241101/pt2/Vessels/*',filename(1:13),'*']);
%     fullfilename_vess = fullfile(vess_dir(1).folder,vess_dir(1).name);
    disp(['doing file #',num2str(i),'/',num2str(length(listdir_mg))])
    
    
    %% process vessel first
    
%     vesselseg = h5read(fullfilename_vess,'/exported_data');
%     vessel_orig = squeeze(vesselseg(1,:,:,:));
%     vessel_orig(vessel_orig==2)=0;
%     noise2 = 50;
%     NoiseIm2=bwareaopen(vessel_orig,noise2);
%     se = strel('disk',3);
%     filled_ves = imfill(NoiseIm2,'holes');
%     dilated_ves = imdilate(filled_ves,se);
%     closed_ves = imclose(dilated_ves,se);
%     eroded_ves = imerode(closed_ves,se);
%     vessel_processed = imgaussfilt(uint8(eroded_ves),2);
    
    %% process microglia
    
    im = h5read(fullfilename,'/exported_data');
    im2 = squeeze(im(1,:,:,:));
    im2(im2>0.5)=1;
    im2(im2<1)=0;
    im2 = uint8(im2);
    maxSlices = 20;
    if size(im2,3)>=maxSlices
        numSlices = size(im2,3); 
        numSplits = ceil(numSlices/maxSlices);
        if numSplits == 2
            splitIdx = ceil(numSlices/2);
            im_new = im2(:, :, 1:splitIdx);
            Microglia  = MGMorph2D(im_new);
            save([erase(listdir_mg(i).name,'.h5'),'_1'],'Microglia','-v7.3')
            im_new2 = im2(:, :, splitIdx+1:end);
            Microglia  = MGMorph2D(im_new2);
            save([erase(listdir_mg(i).name,'.h5'),'_2'],'Microglia','-v7.3')
        elseif numSplits == 3
            splitIdx = ceil(numSlices/3);
            im_new = im2(:, :, 1:splitIdx);
            Microglia  = MGMorph2D(im_new);
            save([erase(listdir_mg(i).name,'.h5'),'_1'],'Microglia','-v7.3')
            im_new2 = im2(:, :, splitIdx+1:2*splitIdx);
            Microglia  = MGMorph2D(im_new2);
            save([erase(listdir_mg(i).name,'.h5'),'_2'],'Microglia','-v7.3')
            im_new3 = im2(:, :, 2*splitIdx+1:end);
            Microglia  = MGMorph2D(im_new3);
            save([erase(listdir_mg(i).name,'.h5'),'_3'],'Microglia','-v7.3')
        end
    else
        Microglia  = MGMorph2D(im2);
        save(erase(listdir_mg(i).name,'.h5'),'Microglia','-v7.3')
    end
end





system('killall caffeinate');
