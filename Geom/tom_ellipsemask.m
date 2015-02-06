function vol=tom_ellipsemask(vol,a,b,c,sigma,center);
%TOM_ELLIPSEMASK masks volume with ellipse of radius r around center
%
%   vol=tom_ellipsemask(vol,a,b,c,sigma,center);
%
%PARAMETERS
%
%  INPUT
%   vol                 input volume
%   a                   radius 1
%   b                   radius 2
%   c                   radius 3
%   sigma               decay
%   center              coord of center
%  
%  OUTPUT
%   vol                 mask volume
%
%EXAMPLE
%   xxx= ones(64,64);
%   yyy = tom_ellipsemask(xxx,4,6,6,10,[16 16 1]);
%   imagesc(yyy);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FF 07/07/03
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

error(nargchk(4,6,nargin))
if (nargin < 6)
    center=[floor(size(vol,1)/2)+1, floor(size(vol,2)/2)+1, floor(size(vol,3)/2)+1];
end;
mask = ones(size(vol,1), size(vol,2), size(vol,3));
limmax = max(max(size(vol,1), size(vol,2)), size(vol,3));
[x,y,z]=ndgrid(-center(1)+1:-center(1)+size(vol,1),-center(2)+1:-center(2)+size(vol,2), ...
    -center(3)+1:-center(3)+size(vol,3));
x = (x./a).^2;
y = (y./b).^2;
z = (z./c).^2;
ind = find( sqrt(x+y+z) > 1);

mask(ind) = 0;
if (nargin > 4) 
    if (sigma > 0)
        mask(ind) = exp(-((sqrt( x(ind) + y(ind) + z(ind)) -1)/sigma).^2);
        ind = find(mask < exp(-2));
        mask(ind) = 0;
    end;
end;
vol = vol.*mask;
