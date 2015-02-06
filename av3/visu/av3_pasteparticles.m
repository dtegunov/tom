function mask = av3_pasteparticles(motl,ref,dimx,dimy,dimz,iobj);
% AV3_MOTL2MASK creates a mask from a MOTL
%   mask = av3_pasteparticles(motl,ref,dimx,dimy,dimz,iobj);
%
%   MOTL        motive list
%   RADIUS      radius of spheres to be used as a mask
%   DIMX        dimension in x of mask
%   DIMY        dimension in y of mask
%   DIMZ        dimension in z of mask
%   IOBJ        index of object (optional). E.g. assume particles are
%               located at different viruses labelled by IOBJ (column 6 in
%               MOTL).
%
%   MASK        generated mask with spheres at coordinates
%
%   The idea is to visualize the locations of particles in the context of
%   your cell or whatever. Current implementation is relatively slow...
%   
%   FF 12/17/03

error(nargchk(4,7,nargin))
if (nargin < 4)
    dimy=dimx;
    dimz=dimx;
end;
if (nargin < 6)
    indx = 1:size(motl,2);
else
    indx = find(motl(6,:) == iobj);
end;
mask = zeros(dimx,dimy,dimz);
tmpvol1 = mask;tmpvol2=mask;
dimhalf = floor(size(ref,1)/2);
for ind=1:size(indx,2)
    class = motl(20,indx(ind));
    if (class == 1)
        x = motl(8,indx(ind));y= motl(9,indx(ind));z= motl(10,indx(ind));
        xshift = motl(11,indx(ind));yshift= motl(12,indx(ind));zshift= motl(13,indx(ind));
        phi = motl(17,indx(ind));psi= motl(18,indx(ind));the= motl(19,indx(ind));
        ref4paste = tom_spheremask(double(tom_rotate(ref,[phi,psi,the])),dimhalf-2,1);
        %ref4paste = tom_shift(ref4paste,[xshift yshift zshift]);
        tmpvol2 = tom_paste(tmpvol1,ref4paste,[x+round(xshift)-dimhalf y+round(yshift)-dimhalf z+round(zshift)-dimhalf ]);
        mask=tmpvol2+mask;
        disp(['added particle no ' num2str(ind)]);
    else
        disp(['particle no ' num2str(ind) ' discarded']);
    end;
end;
