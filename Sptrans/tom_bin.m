function binned = tom_bin(im,nbin)
%TOM_BIN performs binning of 1D, 2D or 3D images
%
%   binned = tom_bin(im,nbin)
%
%PARAMETERS
%
%  INPUT
%   im                  image
%   nbin                number of binnings (default = 1)
%  
%  OUTPUT
%   binned              binned image
%
%EXAMPLE
%   tom_bin(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FF 03/14/03
%   updated by VL 01/13/05 added nbin = 0 case
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
 

error(nargchk(1,2,nargin))
if (nargin < 2)
    nbin = 1;
end;
im=single(im);
for ibin =1:nbin
    imsz = size(im);
    if size(imsz,2) == 2
        imsz(3) = 1;
    elseif size(imsz,2) == 1
        imsz(2) = 1;        
        imsz(3) = 1;
    end
    
    if imsz(3) == 1 && imsz(2) > 1 && imsz(1) > 1
        binned = (im(1:2:imsz(1)-1,1:2:imsz(2)-1)+im(1:2:imsz(1)-1,2:2:imsz(2))+im(2:2:imsz(1),1:2:imsz(2)-1)+im(2:2:imsz(1),2:2:imsz(2)))/4;
    elseif imsz(3) > 1 && imsz(2) > 1 && imsz(1) > 1
        binned = (im(1:2:imsz(1)-1,1:2:imsz(2)-1,1:2:imsz(3)-1)+im(1:2:imsz(1)-1,2:2:imsz(2),1:2:imsz(3)-1)+...
            im(2:2:imsz(1),1:2:imsz(2)-1,1:2:imsz(3)-1)+im(2:2:imsz(1),2:2:imsz(2),1:2:imsz(3)-1) +...
            im(1:2:imsz(1)-1,1:2:imsz(2)-1,2:2:imsz(3))+im(1:2:imsz(1)-1,2:2:imsz(2),2:2:imsz(3))+...
            im(2:2:imsz(1),1:2:imsz(2)-1,2:2:imsz(3))+im(2:2:imsz(1),2:2:imsz(2),2:2:imsz(3)))/8;
    else
        binned = (im(1:2:imsz(1)-1)+im(2:2:imsz(1)))/2;
    end;
    im = binned;
end;

% if no binning
if nbin == 0
    binned = im;
end;

