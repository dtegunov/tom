function image = av3_bandpass(image,low,hi,smooth)
%AV3_BANDPASS performs bandpass filtering of image or volume
%
%   image = av3_bandpass(image,low,hi,sigma)
%   
%   An IMAGE (1, 2, or 3D) is bandpass filtered using FFT. LOW and HI
%   specifies the low- and highpass. Additionally a smoothing parameter
%   SMOOTH can be specified to avoid hard cropping (as in TOM_SPHEREMASK).
%
%   In contrast to TOM_BANDPASS AV3_BANDPASS can also handle volumes with
%   different dimensions in x, y, and z.
%
%PARAMETERS 
%  IN
%   image         image or volume to be filtered
%   low           lowest frequ (in Ny)
%   hi            highest frequ
%   sigma         smoothing parameter (in Ny)
%
%  OUT
%   image         filtered image or volume
%
%SEE ALSO
%   TOM_BANDPASS, TOM_FILTER
%
%    Copyright (c) 2004
%    AV3 toolbox for Electron Tomography
%    Max-Planck-Institute for Biochemistry
%    Dept. Molecular Structural Biology
%    82152 Martinsried, Germany
%    http://www.biochem.mpg.de/tom
%
%   11/17/04 FF
%   last revision 03/31/05 FF - updated docu
error(nargchk(3,4,nargin));
%scf = 1/(size(image,1)*size(image,2)*size(image,3));
scf = 1;
if nargin<4 % fast version
    [x,y,z]=ndgrid( -1:1/(floor(size(image,1)/2)):(size(image,1)-1)*1/(floor(size(image,1)/2))-1,...
        -1:1/(floor(size(image,2)/2)):(size(image,2)-1)*1/(floor(size(image,2)/2))-1,...
        -1:1/(floor(size(image,3)/2)):(size(image,3)-1)*1/(floor(size(image,3)/2))-1);
    len = sqrt(x.^2 +y.^2+z.^2);clear y; clear z;clear x;
    lowp = len <= hi;
    highp=len >= low;
    image = fftshift(tom_fourier(image));
    image = scf*real(tom_ifourier(ifftshift(highp.*lowp.*image)));
else % smoothened version
    image = fftshift(tom_fourier(image));
    if low > 0
        a = hi*(floor(size(image,1)/2));
        b = hi*(floor(size(image,2)/2));
        c = hi*(floor(size(image,3)/2));
        d = low*(floor(size(image,1)/2));
        e = low*(floor(size(image,2)/2));
        f = low*(floor(size(image,3)/2));
        image= scf*real(tom_ifourier(ifftshift(tom_ellipsemask(image,a,b,c,smooth) - tom_ellipsemask(image,d,e,f,smooth))));
    else
        a = hi*(floor(size(image,1)/2));
        b = hi*(floor(size(image,2)/2));
        c = hi*(floor(size(image,3)/2));
        image= scf*real(tom_ifourier(ifftshift(tom_ellipsemask(image,a,b,c,smooth))));
    end;
end;

