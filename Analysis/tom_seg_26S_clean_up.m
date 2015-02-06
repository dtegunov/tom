function [d_6clean,d_3clean,final_mask]=tom_seg_26S_clean_up(d_6fold,d_3fold,pixs,dust_size,dust_size_4_cc,cernel,expand_thr,mass_with_cc,mass_without_cc,mask_diff)
%tom_seg_26S_clean_up 
%   [d_6clean,d_3clean]=tom_seg_26S_clean_up(d_6fold,d_3fold,pixs,dust_size,dust_size_4_cc,cernel,expand_thr)
%PARAMETERS
%
%  INPUT
%   d_6fold             6 fold sym map
%   d_3fold             3 fold sym map
%   pixs                pixelsize
%   dust_size           size of dust for aaa seg
%   dust_size_4_cc      size of dust for coil seg 
%   cernel              cernel for matlab bwopen  
%   expand_thr          thres for expanding coils
%   mass_with_cc        mass 6 fold  
%   mass_without_cc     mass 3 fold   
%
%  OUTPUT 
%   final_mask         mask
%
%
%
%EXAMPLE
%   
% [d_6clean,d_3clean]=tom_seg_26S_clean_up(d_6fol,d_3fold,4.42,2000,80,18,2.3,280,225);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB late june 10
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



cc_mass=mass_with_cc-mass_without_cc;

[mask_6 dust]=tom_remove_dust(d_6fold,mass_without_cc,pixs,dust_size,cernel,20);
 
 
 
 %grow mask a little
 mask_6_grow=tom_filter(mask_6,2)>0.08;
 
 if (exist('mask_diff','var'))
     mask_diff_grow=tom_filter(mask_diff,2)>0.04;
 else
     mask_diff_grow=ones(size(d_3fold));
 end;
 
 [mask_3_rings dust]=tom_remove_dust(d_3fold.*mask_6_grow,mass_without_cc,pixs,dust_size,cernel,20);
 
 
 no_cc_mask=(mask_3_rings==0);
 
 for i=1:3
    final_mask=tom_mass_grow((d_3fold.*no_cc_mask).*mask_diff_grow,mask_3_rings,cc_mass,expand_thr,pixs);
    mask_clean=bwareaopen(final_mask,dust_size_4_cc,cernel);
    dust=final_mask-mask_clean;
    num_of_dust=length(find(dust>0));
    disp([num2str(num_of_dust) ' pixels removed!']);
    if (num_of_dust==0)
        break;
    end;
    no_cc_mask=(((no_cc_mask==0)+(dust))>0)==0;
end;
 
d_3clean=(mask_3_rings+final_mask)>0;
 
d_6clean=mask_6;
 
 