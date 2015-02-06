function [vol_clean dust thresh num_dust_vox]=tom_remove_dust(in,mass,pixs,dust_size,bw_cer,max_iter,verbose)
%TOM_REMOVE_DUST  performs an iterative dust removel and mass estimation 2
%                 obtain the final mass without dust  
%   
%   [vol_clean dust thresh num_dust_vox]=tom_remove_dust(in,mass,pixs,dust_size,bw_cer,max_iter,verbose)
%
%PARAMETERS
%
%  INPUT
%   in                 input volume 
%   mass               mass array or single value
%   pixs               pixelsize in Angstroem
%   dust_size          (200) dust size for matlab function bwareaopen
%   bw_cer             (6) kernel for matlab function bwareaopen
%   max_iter           (30) max number of cleaning iterations 
%   verbose            (1) verbose 0/1   
% 
%  OUTPUT 
%   vol_clean          cleaned volume
%   dust               dust mask 
%   thresh             threshold
%   num_dust_vox       number of dust voxels
%
%
%  NOTE
%   mass is white !!! 
%   the dust removel is performed by bwareaopen
%     
%
%EXAMPLE
%   
%  sp=tom_spheremask(ones(64,64,64),5,6) + tom_spheremask(ones(64,64,64),2,4,[18 18 18]);
%  
%  % example without mass growth 
%  [vol_clean dust]=tom_remove_dust(sp,2500,8.84,400,6); 
%  figure; tom_dspcub(vol_clean); set(gcf,'Name','without');
%
%  % example with mass growth 
%  [vol_clean dust]=tom_remove_dust(sp,1000:100:2500,8.84,400,6);
%  figure; tom_dspcub(vol_clean); set(gcf,'Name','with mass growth 1000:100:2500');
% 
%
%REFERENCES
%
%SEE ALSO
%   
%   bwareaopen,tom_calc_isosurface
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

if (nargin<4)
    dust_size=200;
end;

if (nargin<5)
    bw_cer=6;
end;

if (nargin<6)
    max_iter=30;
end;

if (nargin<7)
    verbose=1;
end;

in=tom_norm(in,'mean0+1std')+100; %...dust is set to 0 so shift mass above

dust=ones(size(in));

for ii=1:length(mass)
    num_of_dust=length(find(dust==0));
    if (verbose==1)
        disp(['Starting Iterative Removal with: ' num2str(mass(ii)) ' and '  num2str(num_of_dust) ' dust voxels']);
    end; 
    
    mask_start=in.*dust;
   
    for i=1:max_iter
        [thresh mass_out vol_out mask_out]=tom_calc_isosurface(mask_start,mass(ii),pixs,0.02);
        mask_clean=bwareaopen(mask_out,dust_size,bw_cer);
        dust=dust.*((mask_out-mask_clean)==0);
        if (verbose==1)
            disp(['   Iter ' num2str(i) ': num of dust pixel: '  num2str(length(find(mask_out-mask_clean))) ' thresh: ' num2str(thresh)]) ;
        end;
        mask_start=in.*dust;
        if (length(find(mask_out-mask_clean))==0)
            break;
        end;
    end;
    mask=mask_clean;
end;


vol_clean=mask;


num_dust_vox=length(find(dust==0));











