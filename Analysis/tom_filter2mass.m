function tom_filter2mass(input_vol,output_vol,mass,pixs,dust_size,bw_cer,max_iter,noise_fill,verbose)
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
%   noise_fill         (1) replace zeros by noise with the same stat as the background   
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
    dust_size=2000;
end;

if (nargin<6)
    bw_cer=6;
end;

if (nargin<7)
    max_iter=5;
end;

if (nargin<8)
    noise_fill=0;
end;

if (nargin<9)
    verbose=1;
end;


in=tom_spiderread(input_vol);

[vol_clean dust thr]=tom_remove_dust(in.Value,mass,pixs,dust_size,bw_cer,max_iter,verbose);


if (noise_fill)
    mask=tom_spheremask(ones(size(in.Value)),round(size(in.Value,1)./2));
    idx_annu=find(mask==0);
    m_anu=mean(in.Value(idx_annu));
    m_std=std(in.Value(idx_annu));
    idx_vol_clean=find(vol_clean==0);
    my_noise=(tom_norm(rand(length(idx_vol_clean),1),'mean0+1std').*m_std)+m_anu;
    out=in.Value;
    out(idx_vol_clean)=my_noise;
else
    out=vol_clean.*in.Value;
end;


tom_spiderwrite(output_vol,out);



