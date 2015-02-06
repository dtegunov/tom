function psd = tom_calc_periodogram_parallel(image,sz,flag,smooth)
%TOM_CALC_PERIODOGRAM_PARALLEL calculates a periodogram to be used in the xmipp ctf
%determination functions parallelized.
%
% psd = tom_calc_periodogram(image,sz)
%
%PARAMETERS
%
%  INPUT
%   image           input 2D image
%   sz              (optional) output size of the averaged periodogram (default 256)
%   flag            (optional) currently not supported
%   smooth          smooth the border of the subregion by 'smooth' pixels
% 
%  OUTPUT
%   psd             periodogram
%
%EXAMPLE
%   psd = tom_calc_periodogram_parallel(image,256,0,256./16);
%
%
%REFERENCES
%
%
%
%SEE ALSO
%   
%   TOM_CALC_PERIODOGRAM, TOM_XMIPP_ADJUST_CTF, TOM_XMIPP_ENHANCE_PSD
%
%   created by AK & COS 27/08/07
%   parallel version by SN 01/10/09
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

if nargin < 4
    smooth = 0;
end

if nargin < 3
    flag = false;
end


if nargin < 2
    sz = 256;
end

image = single(image);

imsz = size(image);
shift = sz/2;

psd = zeros(sz,'single');

counter = 0;
idx=0;
area=0;
for i=1:shift:imsz(1)-sz;
    for j=1:shift:imsz(2)-sz;
        idx=idx+1;
        area(idx,1)=i;
        area(idx,2)=j;
    end
end

parfor idx=1:size(area,1)
    
    im_tmp = image(area(idx,1):area(idx,1)+sz-1,area(idx,2):area(idx,2)+sz-1);
    if smooth > 0
        im_tmp = tom_smooth(im_tmp,smooth);
    end
    psd = psd + abs(fft2(im_tmp)).^2;
    counter = counter + 1.0;
    
end;

psd = (psd./counter)./sz.^4;
