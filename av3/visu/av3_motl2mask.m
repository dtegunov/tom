function mask = av3_motl2mask(motl,radius,dimx,dimy,dimz,rsmooth);
% AV3_MOTL2MASK creates a mask from a MOTL
%
%   mask = av3_motl2mask(motl,radius,dimx,dimy,dimz,rsmooth);
%
% PARAMETERS
%  INPUT
%   MOTL        motive list
%   RADIUS      radius of spheres to be used as a mask
%   DIMX        dimension in x of mask
%   DIMY        dimension in y of mask
%   DIMZ        dimension in z of mask
%   RSMOOTH     smoothing radius to prevent sharp edges
%
%  OUTPUT
%   MASK        generated mask with spheres at coordinates
%   
%   FF 12/17/03

error(nargchk(4,7,nargin))
if (nargin < 4)
    dimy=dimx;
    dimz=dimx;
end;
if (nargin < 7)
    rsmooth = 0;
end;
mask = zeros(dimx,dimy,dimz);
for ind=1:size(motl,2)
    x = motl(8,ind);y= motl(9,ind);z= motl(10,ind);
    mask = mask + tom_spheremask(ones(dimx,dimy,dimz),radius,rsmooth,[x y z]);
    disp(['added particle no ' num2str(ind)]);
end;
mask = tom_limit(mask,0,1);