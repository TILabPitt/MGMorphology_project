# MGMorphology_project

Code was partially adapted from https://github.com/ElisaYork/3DMorph. If used, please cite  "3DMorph automatic analysis of microglial morphology in 3 dimensions from ex vivo and in vivo imaging"

Requirements:

Ilastik for segmentation: https://www.ilastik.org/

Accurate Fast Marching Matlab add-on: https://www.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching



Processing steps:

Step 1:

Create a .h5 file of the maximum intensity projection of the GFP channel for each time frame for Ilastik segmentation. Make sure data is normalized 0-1 : 

load('TSer03141_res.mat')

vessel=squeeze(maxim_intc(:,:,2,:));gfp=squeeze(maxim_intc(:,:,3,:));

vessel=volwlevel(vessel,[],1);gfp=volwlevel(gcamp,[],1);

h5create("20250327_TSer3141_microglia.h5","/im",size(gfp))

h5write("20250327_TSer3141_microglia.h5","/im",gfp)


Step 2:

Open "MyProject2.ilp" in Ilastik which is the trained ilastik file. Go to bottom where it says export new batch and load in new data. Note- if you wish to use the trained Ilastik model, please request the data which was trained for the model or else it will not work.

Step 3:

Go back to Matlab and run "runscript.m". Make sure you update dir variable in line 7 to where new segmentations are located. Note- this function will take a while take a few hours to run if you are processing several time series. Step away and work on something else. 

Step 4:

Once finished, we need to convert results in .mat files to .xlsx. Run "resultstocsv_wrap.m" where you update dir file to the location of where your processed .mat files are located.

Step 5: 

If needed, append your excel tables corresponding the each time series collected in one imaging session. Run "Combine_tables.m" and specify which time series to append together.

Step 6:
Finally, combine all of your results together and specify sex, APOE statues, and file name in "Combine_allresults.m"
