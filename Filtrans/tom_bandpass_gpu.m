function image = tom_bandpass(image,low,hi,smooth)
%TOM_BANDPASS performs bandpass filtering of image or volume
%
%   image = tom_bandpass(image,low,hi,smooth)
%
%PARAMETERS
%
%  INPUT
%   image               iamge or volume to be filtered
%   low                 lowest frequence (in pixels)
%   hi                  highest frequence
%   smooth              smoothing (optional) - if low = 0 no smoothing around
%                           zero frequency
%  
%  OUTPUT
%   image               filtered image or volume
%
%EXAMPLE
%   tom_bandpass(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_FILTER
%
%   created by FF 02/14/03
%   updated by FF 03/31/05
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


error(nargchk(3,4,nargin));

if nargin<4 % fast version
    [x,y,z]=ndgrid(gpuArray(-floor(size(image,1)/2):-floor(size(image,1)/2)+(size(image,1)-1)),...
                   gpuArray(-floor(size(image,2)/2):-floor(size(image,2)/2)+size(image,2)-1), ...
                   gpuArray(-floor(size(image,3)/2):-floor(size(image,3)/2)+size(image,3)-1));
    len = sqrt(x.^2 +y.^2+z.^2);
    clear y; clear z; clear x;
    lowp = ifftshift(len <= hi);
    
    image = fftn(image);
    image = real(ifftn(lowp.*image));
else % smoothened version
    image = fftshift(fftn(image));
    if low > 0
        image= real(ifftn(ifftshift(tom_spheremask_gpu(image,hi,smooth) - tom_spheremask_gpu(image,low,smooth))));
    else
        image= real(ifftn(ifftshift(tom_spheremask_gpu(image,hi,smooth))));
    end;
end;
    
