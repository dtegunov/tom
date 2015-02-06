function [pointidx, pointcoords, distance] = tom_nearestpoint(refpoint, points)
%TOM_NEARESTPOINT determines the nearest point in a point stack
%
%   [pointidx, pointcoords, distance] = tom_nearestpoint(refpoint, points)
%
%nearestpoint determines the nearest point in a point
%stack relative to a reference point
%
%PARAMETERS
%
%  INPUT
%   refpoint            [x,y,z] 3D or 2D coordinates of reference point
%   points              array of 3D or 2D coordinates of points
%  
%  OUTPUT
%   pointidx            index of the nearest point
%   pointcoords         [x,y,z] 3D or 2D coordinates of the nearest point
%   distance            distance from reference point to nearest point
%
%EXAMPLE
%   [pointidx, pointcoords, distance] = tom_nearestpoint([1 1 1], [5 4 3;2 5 7;-4 5 0])
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 1/12/05
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


if nargin ~=2
    error('Two arguments required.');
end

if (max(size(refpoint))==1)

    distancematrix = zeros(1,size(points,1));
    for i=1:max(size(points))
        d = 0;
        d = d + (points(i)-refpoint).^2;
        distancematrix(i) = sqrt(d);
    end

    [distance,pointidx] = min(distancematrix);
    pointcoords = points(pointidx);

    
    
else

    distancematrix = zeros(1,size(points,1));
    for i=1:size(points,1)
        d = 0;
        for j=1:size(points,2)
            d = d + (points(i,j)-refpoint(j)).^2;
        end
        distancematrix(i) = sqrt(d);
    end

    [distance,pointidx] = min(distancematrix);
    pointcoords = points(pointidx,:);

end;

