function  Microglia  = MGMorph2D(im2)%vessel_processed)



noise = 250; %connected pixels smaller than this removed from image
%checkim_size = 2500; % check large cells to make sure you dont need to split into two or three

Microglia.vals = cell(size(im2,3),1);
Microglia.AvgDist = cell(size(im2,3),1);
%Microglia.SkelIdx = cell(size(im2,3),1);
Microglia.WholeSkel = cell(size(im2,3),1);
Microglia.numendpts = cell(size(im2,3),1);
Microglia.MaxBranchLength = cell(size(im2,3),1);
Microglia.MinBranchLength = cell(size(im2,3),1);
Microglia.AvgBranchLength = cell(size(im2,3),1);
Microglia.numbranchpts = cell(size(im2,3),1);
Microglia.PercentMgVol = cell(size(im2,3),1);
Microglia.FullCellComplexity = cell(size(im2,3),1);
%Microglia.vec_norm = cell(size(im2,3),1);
% Microglia.MaxBranchLength_touchingvess = cell(size(im2,3),1);
% Microglia.AvgBranchLength_touchingvess = cell(size(im2,3),1);
% Microglia.TotalBranches_touchingvess = cell(size(im2,3),1);
% Microglia.TotalBranches_nottouchingvess = cell(size(im2,3),1);

for volind = 1:size(im2,3)
    %% Get connected components for each mask
    
    img = im2(:,:,volind);
    sz = size(img);
    sz3 = size(im2,3);
    
    NoiseIm=bwareaopen(img,noise); %Removes objects smaller than set value (in pixels). For 3D inputs, uses automatic connectivity input of 26. Don't want small background dots left over from decreased threshold.
    
    %im_new = manual_touchups(NoiseIm);
    ConnectedComponents=bwconncomp(NoiseIm);
    NumberOfIdentifiedObjects = length(ConnectedComponents.PixelIdxList);
    connected_lengths = cellfun('length',ConnectedComponents.PixelIdxList);
    if max(connected_lengths)>10000
        ex=zeros(sz(1),sz(2));%Create blank image of correct size
        finder = find(connected_lengths==max(connected_lengths));
        NoiseIm(ConnectedComponents.PixelIdxList{finder(1)})=0;
       
        ex(ConnectedComponents.PixelIdxList{finder(1)})=1;
        se=strel('diamond',1);
        newmask = imopen(imerode(ex, se),se);
        newmask = bwareaopen(newmask,noise);
        NoiseIm = NoiseIm+newmask;
        ConnectedComponents=bwconncomp(NoiseIm);
        NumberOfIdentifiedObjects = length(ConnectedComponents.PixelIdxList);
    end

    for m = 1:NumberOfIdentifiedObjects
        Microglia.vals{volind}{m}=ConnectedComponents.PixelIdxList{1,m}; %Write this separated object to new cell array, Microglia in location 'col'.
    end
    
    %% Territorial volume 
    %Uses convhulln to create a 3D polygon around the object's external points.
    %For the total occupied vs unoccupied volume, don't want to exclude any
    %cells/processes. Use Microglia list here, not FullMg.
    
    ConvexVol = zeros(NumberOfIdentifiedObjects,1);
    %progbar = waitbar(0,'Finding territorial volume...');
    for i = 1:NumberOfIdentifiedObjects
        %waitbar (i/numObjSep, progbar);
        [x,y] = ind2sub(sz,[Microglia.vals{volind}{i}]); %input: size of array ind values come from, list of values to convert.
        obj = [y,x]; %concatenate x y z coordinates.
        [k,v] = convhull(obj);
        ConvexVol(i,:) = v;
    end
    
    
    TotMgVol = sum(ConvexVol); %Calculate total area of image covered by microglia. 
    CubeVol = (sz(1)*sz(2)); %area of image cube in um^3. 
    EmptyVol = CubeVol-TotMgVol;%And the remaining 'empty space'.
    Microglia.PercentMgVol{volind} = ((TotMgVol)/(CubeVol))*100;


    %% Volume of Full Cells
    % Determines the cell volume by the number of voxels multiplied by the
    % voxscale to convert into real world units. Also finds the convex
    % territorial volume of only full cells. Cell complexity or extent) is
    % calculated as territorial volume / cell volume and represents how
    % bushy/amoeboid or branched cells are within their territory.
    
    NumberOfPixelsPerCell = cellfun(@numel,Microglia.vals{volind});
    [biggest,idx] = max(NumberOfPixelsPerCell);
    CellVolume = (NumberOfPixelsPerCell);%*voxscale)'; %list of volume of each cell
    %MaxCellVol = biggest*voxscale; %Volume determined by microscope scale, and voxel number reported in cc.PixelIdxList. Should be in um^3
    
    
    % Find the convexvol of only full cells
    FullCellTerritoryVol = zeros(NumberOfIdentifiedObjects,1);
    for i = 1:NumberOfIdentifiedObjects
        [x,y] = ind2sub(sz,[Microglia.vals{volind}{i}]); %input: size of array ind values come from, list of values to convert.
        obj = [y,x]; %concatenate x y z coordinates.
        [k,v] = convhulln(obj);
        FullCellTerritoryVol(i,:) = v;%*voxscale;
    end
    
    % Find the complexity (or extent) of full cells.
    FullCellComplexity = zeros(NumberOfIdentifiedObjects,1);
    for i = 1:length(FullCellTerritoryVol)
    Microglia.FullCellComplexity{volind}{i} = FullCellTerritoryVol(i)/CellVolume(i);
    end
    


    %% Distance Between Centroids
    %Finds location of centroid of each cell, and measures distance between
    %them as a measure of cell density or dispersion.
    
    cent = (zeros(NumberOfIdentifiedObjects,2));

    for i=1:NumberOfIdentifiedObjects
        ex=zeros(sz(1),sz(2));%Create blank image of correct size
        ex(Microglia.vals{volind}{i})=1;%write in only one object to image. Cells are white on black background.
        CentroidCoord = SomaCentroid(ex);
        cent(i,:) = CentroidCoord;
    end

    centum = (zeros(NumberOfIdentifiedObjects,2));
    centum(:,1)=cent(:,1); %*scale; %Convert pixel location to microns so that distances are in correct scale
    centum(:,2)=cent(:,2); %*scale;
    centdist = pdist2(centum,centum); %Calculate distance from each centroid to all other centroids
    centdist=nonzeros(centdist); %Remove all 0s (distance from one centroid to itself)
    Microglia.AvgDist{volind} = mean(centdist);



    kernel = [1 1 1; 1 0 1;1 1 1];
    numendpts = zeros(NumberOfIdentifiedObjects,1);
    numbranchpts = zeros(NumberOfIdentifiedObjects,1);
    MaxBranchLength = zeros(NumberOfIdentifiedObjects,1);
    MinBranchLength = zeros(NumberOfIdentifiedObjects,1);
    AvgBranchLength = zeros(NumberOfIdentifiedObjects,1);

    %% Skeletonize
    for i=1:NumberOfIdentifiedObjects
        ex=zeros(sz(1),sz(2));%Create blank image of correct size
        ex(Microglia.vals{volind}{i})=1;%write in only one object at a time to image. 
        SmoothEx = imgaussfilt(ex); %Smooth the cell so skeleton doesn't pick up many fine hairs
        %skel = bwskel(logical(SmoothEx)); 
        disp(['doing skeletonization of microglia #',num2str(i),'/',num2str(NumberOfIdentifiedObjects),' and time point #' ,num2str(volind),'/', num2str(sz3)])
        
        % seems to have issues at boundary
        maskSize = size(SmoothEx); % Get size of the image
        boundaryTouched = any(SmoothEx(1, :)) || any(SmoothEx(end, :)) || ...
                  any(SmoothEx(:, 1)) || any(SmoothEx(:, end));
        
        % Pad the image if it touches the boundary and apply a bit bigger
        % gauss smooth
        if boundaryTouched
            % Pad the image by 10 pixels on each side
            paddedMask = padarray(SmoothEx, [10 10], 0); % Pad with zeros (background)
            
            fastskel = bwskel(logical(paddedMask));
            fastbranch = bwmorph(fastskel,'branchpoints');
            if length(find(fastbranch(:))) > 60 % struggles with very complex boundaries
                paddedMask = imgaussfilt(paddedMask,2);
                paddedMask = imbinarize(paddedMask);
            else
                paddedMask = imgaussfilt(paddedMask,1);
                paddedMask = imbinarize(paddedMask);
            end

            skeletonizedMask = skeleton(paddedMask);%Find the skeleton! This uses msfm3d and rk4 files, which have been compiled and the .mexw64 versions included. If errors, re-run compilation of these files (in FastMarching_version3b folder), and add the folder and subfolders to path. 
            %Convert cell output of branches into one image for further processing.
            WholeSkel=zeros(size(paddedMask));
            WholeList = round(vertcat(skeletonizedMask{:}));
            SkelIdx = sub2ind(size(paddedMask),WholeList(:,1),WholeList(:,2));
            WholeSkel(SkelIdx)=1;
            % Remove the 10-pixel padding from all sides
            WholeSkel = WholeSkel(11:end-10, 11:end-10);

        else
            % If the mask doesn't touch the boundary, no need to pad
            paddedMask = SmoothEx;

            skeletonizedMask = skeleton(paddedMask);%Find the skeleton! This uses msfm3d and rk4 files, which have been compiled and the .mexw64 versions included. If errors, re-run compilation of these files (in FastMarching_version3b folder), and add the folder and subfolders to path. 

            WholeSkel=zeros(size(ex));
            WholeList = round(vertcat(skeletonizedMask{:}));
            SkelIdx = sub2ind(size(ex),WholeList(:,1),WholeList(:,2));
            WholeSkel(SkelIdx)=1;
        end
        
        %check if skeleton got split up into multiple connected components
        [labelMatrix, numObjects] = bwlabel(WholeSkel);
        if numObjects>1
            WholeSkel = connectpointsskeleton(WholeSkel);
        end

        %Microglia.SkelIdx{volind}{i}=SkelIdx;
        Microglia.WholeSkel{volind}{i} = WholeSkel;
        

        [BoundedSkel, right, left, top, bottom]  = BoundingBoxOfCell(WholeSkel); %Create a bounding box around the skeleton and only analyze this area to significantly increase processing speed. 
        si = size(BoundedSkel);

        % Find endpoints, and trace branches from endpoints to centroid    
        i2 = floor(cent(i,:)); %From the calculated centroid, find the nearest positive pixel on the skeleton, so we know we're starting from a pixel with value 1.
        closestPt = NearestPixel(WholeSkel,i2,1); %scale=1 for now
        i2 = closestPt; %Coordinates of centroid (endpoint of line).
        i2(:,1)=(i2(:,1))-left+1;
        i2(:,2) = (i2(:,2))-bottom+1;
        
        endpts = (conv2(BoundedSkel,kernel,'same')==1)& BoundedSkel; %convolution, overlaying the kernel cube to see the sum of connected pixels.      
        %endpts = bwmorph(BoundedImg, 'endpoints');
        EndptList = find(endpts==1);
        [r,c]=ind2sub(si,EndptList);%Output of ind2sub is row column plane
        EndptList = [r c];
        Microglia.numendpts{volind}{i} = length(EndptList);
        %Microglia.vec_norm{volind}{i} = vecnorm_calc(EndptList,i2);
            
        
        masklist =zeros(si(1),si(2),length(EndptList));
        ArclenOfEachBranch = zeros(length(EndptList),1);

        for j=1:size(EndptList,1)%Loop through coordinates of endpoint.
            i1 = EndptList(j,:); 
            mask = ConnectPointsAlongPath(BoundedSkel,i1,i2);
            masklist(:,:,j)=mask;
            if any(size(mask)==1) & length(find(WholeSkel))<4
                ArclenOfEachBranch(j,1) = NaN;
            else
                % Find the mask length in microns
                pxlist = find(masklist(:,:,j)==1);%Find pixels that are 1s (branch)
                distpoint = reorderpixellist(pxlist,si,i1,i2); %Reorder pixel lists so they're ordered by connectivity
                %Convert the pixel coordinates by the scale to calculate arc length in microns.
                distpoint(:,1) = distpoint(:,1); %If 1024 and downsampled, these scales have been adjusted
                distpoint(:,2) = distpoint(:,2); %If 1024 and downsampled, these scales have been adjusted
                if size(distpoint,1) < 2
                    ArclenOfEachBranch(j,1) = NaN;
                else
                    [arclen,seglen] = arclength(distpoint(:,1),distpoint(:,2));%Use arc length function to calculate length of branch from coordinates
                    ArclenOfEachBranch(j,1)=arclen; %Write the length in microns to a matrix where each row is the length of each branch, and each column is a different cell.
                end
            end
        end


        %Find average min, max, and avg branch lengths
        if ~isempty(ArclenOfEachBranch)
            Microglia.MaxBranchLength{volind}(i) = max(ArclenOfEachBranch);
            Microglia.MinBranchLength{volind}(i) = min(ArclenOfEachBranch);
            Microglia.AvgBranchLength{volind}(i) = mean(ArclenOfEachBranch,"omitnan");

%             overlap_skel = zeros(length(ArclenOfEachBranch),1);
%             for k = 1:size(masklist,3)
%                 mask = masklist(:,:,k);
%                 mask_imgspace = RestoreToOriginalSpace(mask, right, left, top, bottom, sz); % restore branch of skeleton to image space from bonding box space
%                 se = strel('disk',3);
%                 dilated_branch = imdilate(mask_imgspace,se); % dilate branch
%                 mask_skel = ex.*dilated_branch; 
%                 overlap_skel(k) = any(ex(:) & vessel_processed(:)); %see if dilated skeleton is touching vessel
%                 %now I want to get the max arclength of a branch touching a vessel
%                 %and the minimum arclength of a branch not touching a vessel
%                 %and number of branches touching a vessel
%                 %and mean length of branches touching a vessel
%             end
%             if any(overlap_skel)
%                 finder = find(overlap_skel);
%                 finder_none = find(~overlap_skel);
%                 Microglia.MaxBranchLength_touchingvess{volind}(i) = max(ArclenOfEachBranch(finder));
%                 Microglia.AvgBranchLength_touchingvess{volind}(i) = mean(ArclenOfEachBranch(finder));  
%                 Microglia.TotalBranches_touchingvess{volind}(i) = length(ArclenOfEachBranch(finder));
%                 Microglia.TotalBranches_nottouchingvess{volind}(i) = length(ArclenOfEachBranch(finder_none));
%             else
%                 Microglia.MaxBranchLength_touchingvess{volind}(i) = NaN;
%                 Microglia.AvgBranchLength_touchingvess{volind}(i) = NaN;
%                 Microglia.TotalBranches_touchingvess{volind}(i) = NaN;
%                 Microglia.TotalBranches_nottouchingvess{volind}(i) = length(ArclenOfEachBranch); % all branches would not be touching a vessel
%             end

        else
            Microglia.MaxBranchLength{volind}(i) = NaN;
            Microglia.MinBranchLength{volind}(i) = NaN;
            Microglia.AvgBranchLength{volind}(i) = NaN;  
%             Microglia.MaxBranchLength_touchingvess{volind}(i) = NaN;
%             Microglia.AvgBranchLength_touchingvess{volind}(i) = NaN;
%             Microglia.TotalBranches_touchingvess{volind}(i) = NaN;
%             Microglia.TotalBranches_nottouchingvess{volind}(i) = NaN;
        end
        %Microglia.ArclenOfEachBranch{volind}(i) = ArclenOfEachBranch;



        fullmask = sum(masklist,3);%Add all masks to eachother, so have one image of all branches.
        fullmask(fullmask(:,:,:)>3)=4;%So next for loop can work, replace all values higher than 3 with 4. Would need to change if want more than quaternary connectivity.
        

        % Define branch level and display all on one colour-coded image.
        pri = (fullmask(:,:,:))==4;
        sec = (fullmask(:,:,:))==3;
        tert = (fullmask(:,:,:))==2;
        quat = (fullmask(:,:,:))==1;

         % Find branchpoints
        brpts =zeros(si(1),si(2),3);
        for kk=1:2 %For branchpoints not connected to end branches (ie. not distal branches). In fullmask, 1 is branch connected to end point, so anything greater than that is included. 
        temp = (fullmask(:,:,:))>kk;
        tempendpts = (conv2(temp,kernel,'same')==1)& temp; %Get all of the 'distal' endpoints of kk level branches
        brpts(:,:,kk+1)=tempendpts;
        end

        % Find any branchpoints of 1s onto 4s (ie. final branch coming off of main trunk). 
        quatendpts = (conv2(quat,kernel,'same')==1)& quat; %convolution, overlaying the kernel cube onto final branches only.
        quatbrpts = quatendpts - endpts; %Have points at both ends of final branches. Want to exclude any distal points (true endpoints)
        %Only want to keep these quant branchpoints if they're connected to a 4(primary branch). Otherwise, the branch point will have been picked up in the previous for loop. 
        fullrep= fullmask;
        fullrep(fullrep(:,:)<4)=0;%Keep only the 4s, as 4s (don't convert to 1)
        qbpts = fullrep+quatbrpts;%Add the two vectors, so should have 4s and 1s.
        qbpts1 = convn(qbpts,ones([3 3]),'same'); %convolve with cube of ones to get 'connectivity'. All 1s 
        brpts(:,:,1) = (quatbrpts.*qbpts1)>= 5;
        allbranch = sum(brpts,3); %combine all levels of branches
        BranchptList = find(allbranch==1);%Find how many pixels are 1s (branchpoints)
        [r,c]=ind2sub(si,BranchptList);%Output of ind2sub is row column plane
        BranchptList = [r c];
        numbranchpts(i,:) = length(BranchptList);
        Microglia.numbranchpts{volind}(i) = length(BranchptList);

    end
end


    


