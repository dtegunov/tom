function markerset = tom_mark_initMarker(image, nPoints, minDistance, border)


isize = size(image);
[eimage, threshold, bx, by] = edge(image, 'sobel');
eimage(:) = 1;

idx = find(eimage);
nidx = length(idx);
idxx = zeros(1, nidx);
idxy = zeros(1, nidx);

for (i=1:nidx)
    [idxy(i), idxx(i)] = ind2sub(isize,idx(i));
end;

eimage_val = zeros(nidx, 1);
for (i=1:nidx)
    x = idxx(i);
    y = idxy(i);
    eimage_val(i) = abs(bx(y,x)) * abs(by(y,x));
end;

[eimage_val, eimage_validx] = sort(eimage_val, 'descend');

idxx = idxx(eimage_validx);
idxy = idxy(eimage_validx);

markerset = nan(2, 1, nPoints);
nmarkerset = 0;



%figure(1);colormap gray;axis equal off;
%imagesc(image);

%hold on;
%plot(idxx(:), idxy(:), 'g.');



i = 1;
minDistance = minDistance^2;
while (nmarkerset < nPoints && i<=nidx)
    x = idxx(i);
    y = idxy(i);

    if (x <= border  ||  x > isize(2)-border || ...
        y <= border  ||  y > isize(1)-border)
        %1;
    else
        dist = ((markerset(1, 1, 1:nmarkerset) - idxx(i)) .^ 2) + ((markerset(2, 1, 1:nmarkerset) - idxy(i)) .^ 2);
        if (any(dist < minDistance))
            %plot(x, y, 'ro', 'MarkerSize', 5, 'LineWidth', 1);
        else
            nmarkerset = nmarkerset + 1;
            markerset([1 2], 1, nmarkerset) = [ idxx(i) idxy(i) ];
            %plot(x, y, 'bd', 'MarkerSize', 15, 'LineWidth', 3);

        end;
    end;
    i = i+1;
end;


markerset = markerset([1 2], 1, 1:nmarkerset);





