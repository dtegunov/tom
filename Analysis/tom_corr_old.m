function ccf = tom_corr(a,b,flag,mask)
% TOM_CORR computes cross correlation function 
%
%   ccf = tom_corr(a,b,flag)
%
%PARAMETERS
%
%  INPUT
%   a       input array 1 - 1d, 2d or 3d data
%   b       input array 2 - same dimension as array 1
%   flag    can be set to 'norm' to obtain normalized CCF
%   mask    to apply a bandpassfilter
%
%  OUTPUT
%   ccf     cross coreelation function
%
%   TOM_CORR computes cross correlation function according to
%             /
%   ccf(t) = | a(t').*b(t+t') dt'
%            /
%   Therefore, the function is computed in reciprocal space. CCF is real,
%   thus A and B are exepected to be real. If normalization is chosen, A is
%   replaced by (A-mean(A))/std(A).
%
%EXAMPLE
%   im = tom_emread('proteasome.em');
%   ccf = tom_corr(im.Value,im.Value,'norm');
%   % example computes normalized autocorrelation fuction of im.Value. 
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FF 03/30/04
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

error(nargchk(0, 4, nargin, 'struct'))

if (nargin > 2)
    n = (size(a,1)*size(a,2)*size(a,3));
    mn = sum(sum(sum(a)))/n;
    stdv = sqrt(sum(sum(sum(a.*a)))/n - mn^2);
    if (stdv~=0)
        a = (a -mn)/stdv;
    end;
    mn = sum(sum(sum(b)))/n;
    stdv = sqrt(sum(sum(sum(b.*b)))/n - mn^2);
    if  (stdv~=0)
        b = (b - mn)/stdv;
    end;
    end;

if (nargin > 3)
    mask=fftshift(mask);
    ccf = real(ifftshift(tom_ifourier((tom_fourier(b).*mask).*conj((tom_fourier(a).*mask) ))))/(size(a,1)*size(a,2)*size(a,3));%perform correlation
else
    ccf = real(ifftshift(tom_ifourier(tom_fourier(b).*conj(tom_fourier(a)))))/(size(a,1)*size(a,2)*size(a,3));%perform correlation
end;




