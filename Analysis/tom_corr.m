function ccf = tom_corr(a,b,flag,mask,psf_wedge_1,psf_wedge_2)
% TOM_CORR computes cross correlation function 
%
%   ccf = tom_corr(a,b,flag,mask)
%
% PARAMETERS
%   INPUT
%    a                  input array 1 - 1d, 2d or 3d data
%    b                  input array 2 - same dimension as array 1
%    flag               can be set to 'norm' to obtain normalized CCF or
%                       'norm_mask' to calculate CCF applying a mask
%                       (see also Roseman, Ultramicroscopy 94 (2003), 225-236)
%   mask                calculate correlation with mask
%   psf_wedge_1         ...
%   psf_wedge_2         ...
%
%   OUTPUT
%   ccf                 cross correlation function
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
%
%SEE ALSO
%   Roseman, Ultramicroscopy 94 (2003), 225-236 (SN)
%
%   created by FF 03/30/04
%   updated by ...
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

error(nargchk(0, 6, nargin, 'struct'))

if (nargin==2)
    ccf = real(fftshift(tom_ifourier(tom_fourier(b).*conj(tom_fourier(a)))))/(size(a,1)*size(a,2)*size(a,3));%perform correlationend;
end;


if (nargin == 3)
    n = (size(a,1)*size(a,2)*size(a,3)); % number of voxel
    mn = sum(sum(sum(a)))/n; % mean of a 
    %stdv = sqrt(sum(sum(sum(a.*a)))/n - mn^2); % standard deviation of a
    stdv = std2(a);
    if (stdv~=0)
        a = (a -mn)/stdv;
    end;
    mn = sum(sum(sum(b)))/n;
    %stdv = sqrt(sum(sum(sum(b.*b)))/n - mn^2);
    stdv = std2(b);
    if  (stdv~=0)
        b = (b - mn)/stdv;
    end;
    ccf = real(fftshift(tom_ifourier(tom_fourier(b).*conj(tom_fourier(a)))))/(size(a,1)*size(a,2)*size(a,3));%perform correlationend;
end;

if (nargin == 4)
    [a mea st n]=norm_inside_mask(a,mask);
    [b mea st n]=norm_inside_mask(b,mask);
    ccf = real(fftshift(tom_ifourier(tom_fourier(b).*conj(tom_fourier(a)))))./n;%perform correlation
end;
    


if (nargin > 4)

    % wedge particles with their original wedge and crosswise
    %a=tom_ifourier(tom_fourier(a).*tom_fourier(psf_wedge_1));
    %b=tom_ifourier(tom_fourier(b).*tom_fourier(psf_wedge_2));
    a=tom_ifourier(tom_fourier(a).*tom_fourier(psf_wedge_2));
    b=tom_ifourier(tom_fourier(b).*tom_fourier(psf_wedge_1));
    
    
    [a mea st n]=norm_inside_mask(a,mask);
    [b mea st n]=norm_inside_mask(b,mask);
    
    ccf = real(fftshift(tom_ifourier(tom_fourier(b).*conj(tom_fourier(a)))))./n;%perform correlation
    
end;


function [normed_vol mea st n]=norm_inside_mask(vol,mask)

% quick and dirty
% sum(sum(sum(mona_mask.*(mask~=0))))./sum(sum(sum(mask~=0)))
% n=0;
% s=0;
% for x=1:size(vol,1)
%     for y=1:size(vol,2)
%         for z=1:size(vol,3)
%             if mask(x,y,z)~=0
%                 n=n+1;
%                 s=s+vol(x,y,z);
%             end;
%         end;
%     end;
% end;
% mea=s./n;

mask=mask~=0;

n=sum(sum(sum(mask~=0)));
mea=sum(sum(sum((vol.*mask).*(vol~=0))))./n;
st=sqrt(sum(sum(sum((((mask==0).*mea)+(vol.*mask) -mea).^2)))./(n - 1));
normed_vol=((vol-mea)./st).*mask;

% n=0;
% s=0;
% for x=1:size(vol,1)
%     for y=1:size(vol,2)
%         for z=1:size(vol,3)
%             if mask(x,y,z)~=0
%                 n=n+1;
%                 s=s+(vol(x,y,z)-mea).^2;
%             end;
%         end;
%     end;
% end;
% st=sqrt(s./(n));
% 
% normed_vol=((vol-mea)./st).*mask;




