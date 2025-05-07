
listdir_mg=[dir('Results11042024/*1022_TS*'); dir('Results11042024/*1031*');dir('Results11042024/*1105*');dir('Results11042024/*1107*')];
sz=[512,512];
kernel = [1 1 1; 1 0 1;1 1 1];

for k = 1:length(listdir_mg)
    load(fullfile(listdir_mg(k).folder,listdir_mg(k).name))
    Microglia.vec_norm = cell(length(Microglia.vals),1);
    Microglia.mean_theta = cell(length(Microglia.vals),1);
    Microglia.std_vec = cell(length(Microglia.vals),1);
    disp(['doing file #',num2str(k),'/',num2str(length(listdir_mg))])
    for volind = 1:length(Microglia.vals)
        for i = 1:length(Microglia.vals{volind})
            ex=zeros(sz(1),sz(2));%Create blank image of correct size
            ex(Microglia.vals{volind}{i})=1;%write in only one object to image. Cells are white on black background.
            CentroidCoord = SomaCentroid(ex);
            WholeSkel = Microglia.WholeSkel{volind}{i};
            
            [BoundedSkel, right, left, top, bottom]  = BoundingBoxOfCell(WholeSkel); %Create a bounding box around the skeleton and only analyze this area to significantly increase processing speed. 
            si = size(BoundedSkel);
    
            % Find endpoints, and trace branches from endpoints to centroid    
            i2 = floor(CentroidCoord); %From the calculated centroid, find the nearest positive pixel on the skeleton, so we know we're starting from a pixel with value 1.
            closestPt = NearestPixel(WholeSkel,i2,1); %scale=1 for now
            i2 = closestPt; %Coordinates of centroid (endpoint of line).
            i2(:,1)=(i2(:,1))-left+1;
            i2(:,2) = (i2(:,2))-bottom+1;
            
            endpts = (conv2(BoundedSkel,kernel,'same')==1)& BoundedSkel; %convolution, overlaying the kernel cube to see the sum of connected pixels.      
            %endpts = bwmorph(BoundedImg, 'endpoints');
            EndptList = find(endpts==1);
            [r,c]=ind2sub(si,EndptList);%Output of ind2sub is row column plane
            EndptList = [r c];
            EndptList2 = EndptList;
            findermember = ismember(EndptList2,i2,'rows');
            EndptList2(findermember) = []; %remove any endpoints that are the same

            if ~isempty(EndptList2)
                if size(EndptList2,1)>=2
                    [vec_norm,mean_theta,std_theta] = vecnorm_calc(EndptList2,i2);
                    Microglia.vec_norm{volind}{i} = vec_norm;
                    Microglia.mean_theta{volind}{i} = mean_theta;
                    Microglia.std_vec{volind}{i} = std_theta;
                else
                    Microglia.vec_norm{volind}{i} = NaN;
                    Microglia.mean_theta{volind}{i} = NaN;
                    Microglia.std_vec{volind}{i} = NaN;
                end
            else
                Microglia.vec_norm{volind}{i} = NaN;
                Microglia.mean_theta{volind}{i} = NaN;
                Microglia.std_vec{volind}{i} = NaN;
            end
        end
    end
    save([erase(listdir_mg(k).name,'.mat'),'_vecnorms'],'Microglia')
end
