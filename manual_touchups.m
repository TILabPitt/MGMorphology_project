function [ im_new ] = manual_touchups( NoiseIm, finder1 )

% manually inspect images to make sure you don't need to split cells touching each other
    sz = size(NoiseIm);
    ConnectedComponents=bwconncomp(NoiseIm); %returns structure with 4 fields. PixelIdxList contains a 1-by-NumObjects cell array where the k-th element in the cell array is a vector containing the linear indices of the pixels in the k-th object. 26 defines connectivity. This looks at cube of connectivity around pixel.
    numObj = numel(ConnectedComponents.PixelIdxList); %PixelIdxList is field with list of pixels in each connected component. Find how many connected components there are.
    objlength = length(ConnectedComponents.PixelIdxList);    

    inpt = [];
    for m = 1:length(finder1) % visually check each segmentation, report if bad
        ex=zeros(sz(1),sz(2));%Create blank image of correct size
        inpt_find = finder1(m);
        ex(ConnectedComponents.PixelIdxList{1,m})=1;%write in only one object at a time to image.
        figure,show(ex)
        shg
        channel=input('Enter 1 for good, 2 for bad: ');
        if channel==2 % just get index for now
            inpt = [inpt,m]; 
        end
    end
    close all
    

    isdone = false;
    while isdone == false
        if isempty(inpt)
            isdone = true;
        else
            for i = 1:length(inpt)
                inpt_new = finder1(inpt(i));
                ex=zeros(sz(1),sz(2));%Create blank image of correct size
                ex(ConnectedComponents.PixelIdxList{1,inpt_new})=1;%write in only one object at a time to image.
                figure,show(ex)
                shg
                method = input('Enter 1 for erode, 2 for manual: ');

                if method == 1
                    methodgood=2;
                    sph = 1;
                    while methodgood == 2
                       
                        se = strel('sphere',sph);
                        Mask = imerode(ex, se); 
                        Mask = imdilate(Mask, se);
                        figure,show(Mask)
                        shg
                        opts = input('Enter 2 for if still bad, 1 for if good, 3 for switch to manual: ');
                        if opts == 2
                            methodgood = 2;
                            sph = sph+1;
                        elseif opts == 1
                            sub_mask = ex - Mask; % remove from image
                            finder = find(sub_mask(:));
                            NoiseIm(finder) = 0;
                            methodgood = 1;
                        elseif opts == 3
                            mask_reg = bwlabel(selectMask(ex));
                            Mask = ex.*mask_reg;
                            sub_mask = ex - Mask; % remove from image
                            finder = find(sub_mask(:));
                            NoiseIm(finder) = 0;
                            methodgood = 1;
                        end
                    end
                    isdone = true;

                elseif method == 2
                    mask_reg = bwlabel(selectMask(ex));
                    Mask = ex.*mask_reg;
                    sub_mask = ex - Mask; % remove from image
                    finder = find(sub_mask(:));
                    NoiseIm(finder) = 0;
                    isdone = true;
                end
            end
        end
    end
    im_new = NoiseIm;
end