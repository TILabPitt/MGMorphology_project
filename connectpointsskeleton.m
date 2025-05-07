function skeletonMask = connectpointsskeleton(skeletonMask)


[labelMatrix, numObjects] = bwlabel(skeletonMask);

% Check if there are exactly two connected components
if numObjects ~= 2
    error('The skeleton does not contain exactly two components.');
end

% Step 2: Extract the two separate skeleton components
component1 = (labelMatrix == 1);
component2 = (labelMatrix == 2);

% Step 3: Compute the Euclidean Distance Transform of the second component
[distToComp2, idxNearestInComp2] = bwdist(component2); % Distance to component 2

% Step 4: Find the coordinates of the nearest point in component 1
[y1, x1] = find(component1); % Coordinates of component 1 pixels
[minDist, nearestIdxInComp1] = min(distToComp2(component1)); % Nearest point in component 1
point1 = [y1(nearestIdxInComp1), x1(nearestIdxInComp1)]; % Coordinates of the nearest point in component 1

% Step 5: Get the corresponding nearest point in component 2 using the index from bwdist
[y2, x2] = ind2sub(size(component2), idxNearestInComp2(y1(nearestIdxInComp1), x1(nearestIdxInComp1))); 
point2 = [y2, x2]; % Coordinates of the nearest point in component 2

% Step 6: Use Bresenham's line algorithm to draw the shortest path between the points
lineCoords = drawLine(point1(2), point1(1), point2(2), point2(1));

% Step 7: Add the line to the skeleton to connect the components
for k = 1:size(lineCoords, 1)
    skeletonMask(lineCoords(k, 2), lineCoords(k, 1)) = 1; % Set the path pixels to 1 (skeletonized)
end

end