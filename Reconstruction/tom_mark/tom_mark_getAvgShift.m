function [shift_mean, shift_std] = tom_mark_getAvgShift(markerset)
%TOM_MARK_GETAVGSHIFT returns a matrix of shifts between images.
% 
%   [shift_mean, shift_std] = TOM_MARK_GETAVGSHIFT(markerset)
%
%   Given a markerset (i.e. an array of size (2 x nImage x nMarkers)
%   the average displacement between every two images is computed.
%   
%   Can be used together with the function TOM_MARK_FINDMATCHES
%
%PARAMETERS
%  INPUT
%     markerset: The coordinates of the marker in the corresponding image.
%         If a position is unknown, its pixelcoordinate must be set to <0
%         or to NaN. The mean is computed only over the known
%         markerpositions.
%         To get a result with every displacement set to unknown, markerset
%         must be a single integer: the number of images.
%  OUTPUT:
%     shift_mean: Is a matrix of ( nImage, nImage, 2 ) where the image #i
%         and #j are displaced by [dx, dy] = shift_mean(#i, #j, [1 2])
%     shift_std: A matrix with the deviations from the known displacements.
%
%   See also TOM_MARK_FINDMATCHES
%
%REFERENCES
%
%   20/03/07
%
%   created by Thomas Haller
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msize = size(markerset);
if (length(msize)==2)
    if (msize(1)==1 && msize(2)==1) 
        msize = [2, markerset, 0];
    else
        msize(3) = 1;
    end;
end;

shift_mean = nan(msize(2), msize(2), 2);
shift_dev = nan(msize(2), msize(2), 2);


if (msize(3)<1)
    return;
end;

markerset(markerset < 0) = nan;

for (i=1:msize(2))
    shift_mean(i,i,1:2) = 0;
    shift_std(i,i,1:2) = 0;
    for (j=(i+1):msize(2))
        xdiff = markerset(1, j, :) - markerset(1, i, :);
        ydiff = markerset(2, j, :) - markerset(2, i, :);

        idx = ~isnan(xdiff) & ~isnan(ydiff);
        xdiff = xdiff(idx);
        ydiff = ydiff(idx);

        if (~isempty(xdiff))
            shift_mean(i, j, 1) = mean(xdiff);
            shift_mean(i, j, 2) = mean(ydiff);
            shift_mean(j, i, 1) = -shift_mean(i, j, 1);
            shift_mean(j, i, 2) = -shift_mean(i, j, 2);

            shift_std(i, j, 1) = std(xdiff);
            shift_std(i, j, 2) = std(ydiff);
            shift_std(j, i, 1) = shift_std(i, j, 1);
            shift_std(j, i, 2) = shift_std(i, j, 2);
%        else
%            shift_mean(i, j, [1 2]) = nan;
%            shift_mean(j, i, [1 2]) = nan;
%            shift_std(i, j, [1 2]) = nan;
%            shift_std(j, i, [1 2]) = nan;
        end;
    end;
end;






