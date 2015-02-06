function binned = tom_downscale(im,nscale)
%TOM_DOWNSCALE performs downscaling of 1D, 2D or 3D images. In contrast to
%tom_bin downscale factor is not necessarily power of 2.
%
%   binned = tom_downscale(im,nscale)
%
%PARAMETERS
%
%  INPUT
%   im                  image
%   nbin                downscale factor
%  
%  OUTPUT
%   binned              downscaled image
%
%EXAMPLE
%   tom_downscale(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_BIN
%
%   created by FF 02/18/13
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2013
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom
 

im=single(im);

imsz = size(im);
if size(imsz,2) == 2
    imsz(3) = 1;
elseif size(imsz,2) == 1
    imsz(2) = 1;
    imsz(3) = 1;
end

if imsz(3) == 1 && imsz(2) > 1 && imsz(1) > 1
    binned = zeros(ceil(imsz(1)/nscale),ceil(imsz(2)/nscale));
    for ix=1:nbin
        for iy=1:nbin
            binned(1:ceil((imsz(1)-ix+1)/nscale),1:ceil((imsz(2)-iy+1)/nscale)) = ...
                binned(1:ceil((imsz(1)-ix+1)/nscale)) + im(ix:nscale:imsz(1),iy:nscale:imsz(2));
        end;
    end;
    % normalize
    binned = binned/(nscale^2);
    xmod = mod(imsz(1),nscale);
    ymod = mod(imsz(2),nscale);
    if xmod > 0
        binned(ceil(imsz(1)/nscale),:) = binned(ceil(imsz(1)/nscale),:)*(nscale/(xmod));
    end;
    if ymod > 0
        binned(:,ceil(imsz(2)/nscale)) = binned(:,ceil(imsz(2)/nscale))*(nscale/(ymod));
    end;
elseif imsz(3) > 1 && imsz(2) > 1 && imsz(1) > 1
    binned = zeros(ceil(imsz(1)/nscale),ceil(imsz(2)/nscale),ceil(imsz(2)/nscale));
    for ix=1:nscale
        for iy=1:nscale
            for iz=1:nscale
                binned(1:ceil((imsz(1)-ix+1)/nscale),1:ceil((imsz(2)-iy+1)/nscale),1:ceil((imsz(3)-iz+1)/nscale)) = ...
                    binned(1:ceil((imsz(1)-ix+1)/nscale),1:ceil((imsz(2)-iy+1)/nscale),1:ceil((imsz(3)-iz+1)/nscale))...
                    + im(ix:nscale:imsz(1),iy:nscale:imsz(2),iz:nscale:imsz(3));
            end;
        end;
    end;
    % normalize 
    binned = binned/(nscale^3);
    xmod = mod(imsz(1),nscale);
    ymod = mod(imsz(2),nscale);
    zmod = mod(imsz(3),nscale);
    if xmod > 0
        binned(ceil(imsz(1)/nscale),:,:) = binned(ceil(imsz(1)/nscale),:,:)*(nscale/(xmod));
    end;
    if ymod > 0
        binned(:,ceil(imsz(2)/nscale),:) = binned(:,ceil(imsz(2)/nscale),:)*(nscale/(ymod));
    end;
    if zmod > 0
        binned(:,:,ceil(imsz(3)/nscale)) = binned(:,:,ceil(imsz(3)/nscale))*(nscale/(zmod));
    end;
else
    binned = zeros(ceil(imsz(1)/nscale),1);
    for ix=1:nscale
        binned(1:ceil((imsz(1)-ix+1)/nscale)) = binned(1:ceil((imsz(1)-ix+1)/nscale))...
            + im(ix:nscale:imsz(1));
    end;
    % normalize
    binned = binned/nscale;
    xmod = mod(imsz(1),nscale);
    if xmod > 0
        binned(ceil(imsz(1)/nscale)) = binned(ceil(imsz(1)/nscale))*(nscale/(xmod));
    end;
end;
