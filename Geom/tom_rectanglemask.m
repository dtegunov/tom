function im=tom_rectanglemask(im,size_m,sigma,center)
%TOM_RECTANGLEMASK creates a rectengular 2d mask
%
%   im=tom_rectanglemask(im,size_m,sigma,center)
%
%PARAMETERS
%
%  INPUT
%   im                 input image
%   size_m             size vector of the rectangle
%   sigma              sigma for smoothing
%   center             (optional) center of the rectangle
%  
%  OUTPUT
%   im                 outputmask
%
%EXAMPLE
%   im=tom_rectanglemask(zeros(56,56),[20 30],2);
%   figure; tom_imagesc(im);
%
%REFERENCES
%
%SEE ALSO
%   tom_filtergui, tom_apply_filter 
%
%   created by fb dont know
%   updated by fb 12/21/06
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
    center=[floor(size(im,1)/2)+1, floor(size(im,2)/2)+1, floor(size(im,3)/2)+1];
end;

startx=center(1)-floor(size_m(1)./2);
stopx=center(1)+floor(size_m(1)./2);
starty=center(2)-floor(size_m(2)./2);
stopy=center(2)+floor(size_m(2)./2);

if (size(im,3)>2)
    startz=center(3)-floor(size_m(3)./2);
    stopz=center(3)+floor(size_m(3)./2);
    im(startx:stopx,starty:stopy,startz:stopz)=1;
else
    im(startx:stopx,starty:stopy)=1;
end;

if (nargin > 2)
    if (sigma~=0)
        im=tom_filter(im,sigma);
    end;
end;