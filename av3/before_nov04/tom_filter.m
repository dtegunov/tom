function im = tom_filter(im,radius,flag)
% TOM_FILTER convolutes with spherical or quadratic kernel
%
%   filtered = tom_filter(im,radius,flag)
%
%   This procedure filters using FFT.
%
% INPUT
%   IM          one, two or 3D volume
%   RADIUS      radius of kernel for spherical kernel
%                   for quadratic kernel: side length of kernel
%   FLAG        'circ' for spherical kernel, 'quadr' for quadratic kernel
%
% OUTPUT
%   FILTERED    filtered image
%
% FF 10/22/03
error(nargchk(2,3,nargin));
if (nargin < 3)
    flag = 'circ';
end;
switch lower(flag)
    case 'circ'
        mask=ones(size(im,1),size(im,2),size(im,3));
        mask=tom_spheremask(mask,radius);
    case 'quadr'
        mask=zeros(size(im,1),size(im,2),size(im,3));
        if size(im,3) > 1
            cent=[floor(size(im,1)/2)+1 floor(size(im,2)/2)+1 floor(size(im,3)/2)+1];
            mask(cent(1)-floor(radius/2):cent(1)-floor(radius/2)+radius-1,...
                cent(2)-floor(radius/2):cent(2)-floor(radius/2)+radius-1,cent(3)-floor(radius/2):cent(3)-floor(radius/2)+radius-1)=1;
        elseif ((size(im,3) == 1) &  (size(im,2) > 1))
            cent=[floor(size(im,1)/2)+1 floor(size(im,2)/2)+1];
            mask(cent(1)-floor(radius/2):cent(1)-floor(radius/2)+radius-1,...
                cent(2)-floor(radius/2):cent(2)-floor(radius/2)+radius-1)=1;
        else
            cent=[floor(size(im,1)/2)+1];
            mask(cent(1)-floor(radius/2):cent(1)+radius-1)=1;
        end;
end;
npix=sum(sum(sum(mask)));
im=real(fftshift(tom_ifourier(tom_fourier(mask).*tom_fourier(im))/npix));