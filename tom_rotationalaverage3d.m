function [ rotaverage ] = tom_rotationalaverage3d( map )

coordmin = -size(map, 1) / 2;
coordmax = size(map, 1) / 2 - 1;
maxr = size(map, 1) / 2;

[x, y, z] = ndgrid(coordmin:coordmax, coordmin:coordmax, coordmin:coordmax);
r = round(sqrt(x.^2+y.^2+z.^2));
r = min(r, maxr);

rotaverage = zeros(maxr + 1, 1);
for i = 0:maxr
    rotaverage(i + 1) = sum(map(r == i));
end;

rotaverage = rotaverage(1:end-1);

end

