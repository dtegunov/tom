function tom_filter2mass2(input_vol,output_vol,mass,pixs,dust_size,bw_cer,max_iter,verbose)
%tom_filter2mass filters a volume with a binary mask according 2 the mass of the structure 
%   
%    tom_filter2mass(input_vol,output_vol,mass,pixs,dust_size,bw_cer,max_iter,verbose)
%
%PARAMETERS
%
%  INPUT
%   input_vol          input volume  
%   output_vol         output volume
%   mass               mass array or single value
%   pixs               pixelsize in Angstroem
%   dust_size          (200) dust size for matlab function bwareaopen
%   bw_cer             (6) kernel for matlab function bwareaopen
%   max_iter           (30) max number of cleaning iterations 
%   verbose            (1) verbose 0/1   
% 
%  OUTPUT 
%
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
%  tom_spiderwrite('sp.vol',sp);
%  tom_filter2mass('sp.vol','sp_f.vol',2500,8.84)
%
%
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



if (nargin<5)
    dust_size=200;
end;

if (nargin<6)
    bw_cer=6;
end;

if (nargin<7)
    max_iter=30;
end;

if (nargin<8)
    verbose=1;
end;

%mask=tom_spheremask(ones(256,256,256),105); % not working!!!

mask=tom_spiderread('/fs/gpfs01/lv01/gpfs/gpfs-baumei/nickell/26S/spwt/rec/classify_11__42_chi128/mask_bin256_big_thr7_upright.spi'); mask=mask.Value;

% mask=tom_rotate(tom_cylindermask(ones(256,256,256),79),[270 90 90])>0;
% mask(1:12,:,:)=0;
% mask(end-12:end,:,:)=0;
% idx_annu=find(mask==0);
in=tom_spiderread(input_vol);


%[vol_clean dust]=tom_remove_dust(in.Value,mass,pixs,dust_size,bw_cer,max_iter,verbose);

%mask_big=(tom_filter(vol_clean,25)>0.003); % working !!!
%mask_big=mask;



%idx_annu=find(mask_big==0);

%m_anu=mean(in.Value(idx_annu));
%m_std=std(in.Value(idx_annu));
% 
% 
% 
% idx_no_mask=find((vol_clean==0));
% 
%my_noise=(tom_norm(rand(length(idx_annu),1),'mean0+1std').*m_std)+m_anu;
% 
out=tom_norm(tom_filter(mask,50).*in.Value,'mean0+1std');
% 

%out(idx_annu)=my_noise;

%out=tom_norm(in.Value,'mean0+1std',vol_clean).*vol_clean;

%out=tom_permute_bg(in.Value,tom_filter(vol_clean,16)>0.01);

%out=in.Value.*mask;

tom_spiderwrite([input_vol '.org'],in.Value);

tom_spiderwrite(output_vol,out);






