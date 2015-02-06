function new_mask=tom_grow_mask(mask,factor,max_error,filter_cer)
%TOM_GROW_MASK grows a give binary mask by a certain grow factor
%
%   new_mask=tom_grow_mask(mask,factor)
%
%PARAMETERS
%
%  INPUT
%   mask                binary mask
%   factor              factor ...determines the new number of voxels
%   max_error           (2) allowed error in percent of input mask 
%   filter_cer          (2% of mask size) in pixels
%
%  OUTPUT
%   vol                 masked volume
%
%EXAMPLE
%   xxx= ones(64,64,64);
%   yyy = tom_spheremask(xxx,10);
%   tom_dspcub(yyy);
%    
%   %to make the mask 10 percent bigger use factor 1.1
%   new_mask=tom_grow_mask(yyy,1.3)  
%   figure; tom_dspcub(new_mask);
% 
%   length(find(new_mask>0))./length(find(yyy>0)) 
%
%
%NOTE:
%
% new numb of voxels = old num of voxels .* factor
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   fb april 2012
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
    max_error=2;
end;

if (nargin < 4)
    fact4filt_cer=size(mask,1)./128;
    filter_cer=round(3.*fact4filt_cer);
end;


thr_inc=0.1;
thr_tmp=0.4;

org_num_of_vox=length(find(mask > 0));
dust_size=round(org_num_of_vox-(0.3.*org_num_of_vox));

max_itr=1000;


mask_filt=tom_filter(mask,filter_cer);

for ii=1:30
    thr_inc=thr_inc./ii;
    thr_start=thr_tmp;
    for i=1:max_itr
        thr_tmp=thr_start-(thr_inc.*i);
        new_num=length(find(bwareaopen(mask_filt > thr_tmp,dust_size,6) ));
        if (new_num > (org_num_of_vox * factor))
            break;
        end;
    end;
    act_error= ((new_num-(org_num_of_vox * factor))./org_num_of_vox) .* 100;
    thr_tmp=thr_start-(thr_inc.*(i-1));
    %disp(num2str(act_error));
    if (act_error < max_error)
        break;
    end;
end;
new_mask=bwareaopen(mask_filt>(thr_start-(thr_inc.*i)),dust_size,6);
    





