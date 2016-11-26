function vol=tom_spheremask(vol, radius,sigma,center)
%TOM_SPHEREMASK masks volume with sphere of radius r around center
%
%   vol=tom_spheremask(vol, radius,sigma,center);
%
%PARAMETERS
%
%  INPUT
%   vol                 volume
%   radius              radius determining radius of sphere
%   sigma               smoothing of mask; if entered mask will be smoothened;
%                       every voxel outside radius gets smoothened by a gaussian
%   center              vector determining center of sphere
%  
%  OUTPUT
%   vol                 masked volume
%
%EXAMPLE
%   xxx= ones(64,64);
%   yyy = tom_spheremask(xxx,4,10,[16 16 1]);
%   imagesc(yyy);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FF 08/14/02
%   updated by FF25/03/04
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

error(nargchk(1,4,nargin))
if (nargin < 4)
    center=[floor(size(vol,1)/2)+1, floor(size(vol,2)/2)+1, floor(size(vol,3)/2)+1];
end;

[x,y,z]=ndgrid(gpuArray(0:size(vol,1)-1), gpuArray(0:size(vol,2)-1), gpuArray(0:size(vol,3)-1));

if (nargin < 2)
    radius = floor((min(min(size(vol,1),size(vol,2)),size(vol,3))-1)/2) ;
end;

x = sqrt((x+1-center(1)).^2+(y+1-center(2)).^2+(z+1-center(3)).^2);
clear y z;

mask = gpuArray.ones(size(vol,1), size(vol,2), size(vol,3));
mask(x >= radius) = 0;

if (nargin > 2) 
    if (sigma > 0)
%         mask(ind) = exp(-((x(ind) -radius)/sigma).^2);
%         ind = find(mask < exp(-2));
%         mask(ind) = 0;
        x = max(0, min((x - radius) / sigma, 1));
        x = 0.5.*(1 + cos(x.*pi));
        mask = x;
    end;
end;

vol = vol.*mask;
