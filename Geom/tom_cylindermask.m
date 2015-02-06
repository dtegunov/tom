
function vol = tom_cylindermask(vol, radius,sigma,center)
%TOM_CYLINDERMASK masks volume with a cylinder of radius r around center
%
%   vol = tom_cylindermask(vol, radius,sigma,center)
%
%PARAMETERS
%
%  INPUT
%   vol                 volume
%   radius              radius determining radius of cylinder (in 2D!)
%   sigma               smoothing of mask; if entered mask will be smoothened;
%                       every voxel outside radius gets smoothened by a gaussian
%                       function exp(-((r-radius)/simga)^2)
%   center              vector determining center of sphere (in x-y plane!)
%  
%  OUTPUT
%   vol                 masked volume
%
%EXAMPLE
%   xxx= ones(64,64,64);
%   yyy = tom_cylindermask(xxx,4,10,[16 16]);
%   tom_dspcub(yyy);
%
%REFERENCES
%
%SEE ALSO
%   tom_spheremask, tom_sphere, tom_ellipsemask
%
%   created by FF 12/09/03
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
    center=[floor(size(vol,1)/2)+1, floor(size(vol,2)/2)+1];
end;
if (nargin < 2)
    radius = floor(min(size(vol,1),size(vol,2))/2) ;
end;
[x,y]=ndgrid(0:size(vol,1)-1,0:size(vol,2)-1);
x = sqrt((x+1-center(1)).^2+(y+1-center(2)).^2);
ind = find(x>=radius);
clear y;
mask = ones(size(vol,1), size(vol,2));
mask(ind) = 0;
if (nargin > 2) 
    if (sigma > 0)
        mask(ind) = exp(-((x(ind) -radius)/sigma).^2);
        ind = find(mask < exp(-2));
        mask(ind) = 0;
    end;
end;
for iz=1:size(vol,3)
    vol(:,:,iz) = vol(:,:,iz).*mask;
end;