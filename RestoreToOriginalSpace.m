function OriginalSpaceImg = RestoreToOriginalSpace(BoundedImg, right, left, top, bottom, originalSize)
    % RestoreToOriginalSpace
    % This function takes a bounded (cropped) binary mask and positions it
    % back into the original image space.
    % INPUT:
    %   - BoundedImg: The cropped binary mask image.
    %   - right, left, top, bottom: Bounding box coordinates from original
    %     image used to crop the BoundedImg.
    %   - originalSize: Size of the original image [rows, cols].
    % OUTPUT:
    %   - OriginalSpaceImg: The restored mask in the original image space.
    % Initialize an empty image of the original size
    OriginalSpaceImg = zeros(originalSize);
    % Insert the bounded image back to its original position
    OriginalSpaceImg(left:right, bottom:top) = BoundedImg;
end