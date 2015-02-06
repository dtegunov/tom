function mask = av3_motl2mask(motl,radius,dimx,dimy,dimz);
% AV3_MOTL2MASK creates a mask from a MOTL
%   mask = av3_motl2mask(motl,radius,dimx,dimy,dimz);
error(nargchk(3,5,nargin))
if (nargin < 4)
    dimy=dimx;
    dimz=dimx;
end;
for ind=1:size(motl,2)
    x = motl(8,ind);y= motl(9,ind);z= motl(10,ind);
    mask = mask + tom_spheremask(ones(dimx,dimy,dimz),radius);
    
end;