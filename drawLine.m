function lineCoords = drawLine(x1, y1, x2, y2)
    % Bresenham's line algorithm to get the coordinates of the line between (x1, y1) and (x2, y2)
    lineCoords = [];
    steep = abs(y2 - y1) > abs(x2 - x1);
    if steep
        [x1, y1] = deal(y1, x1);
        [x2, y2] = deal(y2, x2);
    end
    if x1 > x2
        [x1, x2] = deal(x2, x1);
        [y1, y2] = deal(y2, y1);
    end
    dx = x2 - x1;
    dy = abs(y2 - y1);
    error = dx / 2;
    ystep = -1;
    if y1 < y2
        ystep = 1;
    end
    y = y1;
    for x = x1:x2
        if steep
            lineCoords = [lineCoords; [y, x]]; %#ok<AGROW>
        else
            lineCoords = [lineCoords; [x, y]]; %#ok<AGROW>
        end
        error = error - dy;
        if error < 0
            y = y + ystep;
            error = error + dx;
        end
    end
end