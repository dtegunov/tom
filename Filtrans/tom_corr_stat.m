function [vol_corr vol_corr_bin]=tom_corr_stat(vol,cer,mask,thr4corr,dust_size)
%tom_corr_stat corrects statitics of a volume
%
%   vol_corr=tom_corr_stat(vol,cer,mask,dust_size)
%
%PARAMETERS
%
%  INPUT
%   vol                 volume
%   cer                 cernel size
%   mask                opt.
%   thr4corr            iso thr 4 corr volume
%   dust_size           dust size
%                           
%
%  OUTPUT
%   vol_corr            corrected volume 
%
%EXAMPLE
%   
% corr=tom_corr_stat(vv,50,mask);
%
%REFERENCES
%
%SEE ALSO
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

if (nargin < 4)
    thr4corr='';
end;

if (nargin < 5)
    dust_size='';
end;



disp(' ');

cer_mask=tom_spheremask(ones(size(vol)),cer);

m_vol=tom_os3_mean(vol,cer_mask);
m_std=tom_os3_std(vol,m_vol,cer_mask);

vol_corr=(vol-m_vol)./m_std;

vol_corr=vol_corr.*mask;

if (isempty(thr4corr)==0)
    vol_corr=vol_corr>thr4corr;
end;


if (isempty(dust_size)==0)
    vol_corr=bwareaopen(vol_corr,dust_size,6);
end;



