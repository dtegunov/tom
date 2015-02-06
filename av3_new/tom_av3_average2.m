function avg=tom_av3_average2(vols,psfs,avg_out)
%TOM_AV3_AVERAGE2 adds volumes in vol_list and weights every volume with
%the corresponding in psf_list
%
%   avg=tom_av3_average2(vol_list,psf_list,avg_out)
%
%PARAMETERS
%
%  INPUT
%   vol_list            can be a wildcard or a cell
%   psf_list            can be a wildcard or a cell
%   avg_out             filename of the output average
%
%  
%  OUTPUT
%   avg           avergae
%
%
%EXAMPLE
%   read image
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN/FB 01/24/06
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

if (nargin < 3)
    avg_out='';
end;


%parse inputs

if (ischar(vols))
    dd=dir(vols);
    root_p=fileparts(vols);
    if (isempty(root_p))
        root_p='.';
    end;
    for i=1:length(dd)
        vol_list{i}=[root_p '/' dd(i).name];
    end;
else
    vol_list=vols;
end;

if (ischar(psfs))
    dd=dir(psfs);
    root_p=fileparts(psfs);
    if (isempty(root_p))
        root_p='.';
    end;
    for i=1:length(dd)
        psf_list{i}=[root_p '/' dd(i).name];
    end;
else
    psf_list=psfs;
end;

em_flag=0;
spider_flag=0;

if (length(psf_list)~=length(vol_list))
   disp([num2str(length(vol_list)) ' vols found'] );
   disp([num2str(length(psf_list)) ' psfs found'] );
  % error(['number of vols and psfs has 2 be the same']);
end;

if (tom_isemfile(vol_list{1}))
    v=tom_emreadc(vol_list{1});
end;
if (tom_isspiderfile(vol_list{1}))
    v=tom_spiderread(vol_list{1});
end;

if (tom_isemfile(psf_list{1}))
    psf=tom_emreadc(psf_list{1});
end;
if (tom_isspiderfile(psf_list{1}))
    psf=tom_spiderread(psf_list{1});
    spider_flag=1;
end;


avg=tom_apply_weight_function(v.Value,psf.Value);

disp([vol_list{1} ' * ' psf_list{1} ' + ']);

for i=2:length(vol_list)
    if (tom_isemfile(vol_list{i}))
        v=tom_emreadc(vol_list{i});
    end;
    if (tom_isspiderfile(vol_list{i}))
        v=tom_spiderread(vol_list{i});
    end;
    if (tom_isemfile(psf_list{i}))
        psf=tom_emreadc(psf_list{i});
    end;
    if (tom_isspiderfile(psf_list{i}))
        psf=tom_spiderread(psf_list{i});
    end;
    
    avg=avg+tom_apply_weight_function(v.Value,psf.Value);
    disp([vol_list{i} ' * ' psf_list{i} ' + ']);
end;


if (isempty(avg_out)==0)
    disp([' = ' avg_out]);
    if (em_flag==1)
        tom_emwrite(avg_out,avg);
    end;
    if (spider_flag==1)
        tom_spiderwrite(avg_out,avg);
    end;
end;   











