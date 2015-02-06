function tom_analyse_3d_res(res_map,samp,vol,flag,param,output_folder)
% TOM_ANALYSE_3D_RES analyses 3d resolution map 
%  
%     [res mat_out]=tom_compare_deluxe(map,flag,param)
%  
%  TOM_ANALYSE_3D_RES by building res interval volumes
%  resolution map yes!!
%
%
%  PARAMETERS
%  
%    INPUT
%     res_map                resolution map 
%     samp                   sampling of res map
%     vol                    binarized volume desity map (use tom_calc_isosurface)  
%     flag                   flag for buildin res groups (num_of_groups,increment or mygroups)
%                              
%     param                  number of groups or size of interval e.g. 5 in
%                            Ang
%     output_folder          folder for output 
%
%                          
%    
%    OUTPUT
%
%    NOTE:
%    param for mygroup is a 2 dim matrix with start and stop
%      
%     v(1,:)=[5 10];
%     v(2,:)=[11 15];
%     v(3,:)=[16 20]; 
%     v(4,:)=[25 100];
%
%  EXAMPLE
%  
%  v(1,:)=[0 5]; 
%  v(2,:)=[5 10];
%  v(3,:)=[10 15];
%  v(4,:)=[15 20]; 
%  v(5,:)=[20 25]; 
%  v(6,:)=[25 10000];
% 
%  vol=tom_spiderread('proj_match_split_1.vol');
%  [thresh mass vol vol_bin]=tom_calc_isosurface(vol.Value,2500,4.42,.01);
%  res_map=tom_emread('vol_out_kern10.em'); 
%  tom_analyse_3d_res(res_map.Value,2,vol_bin,'mygroups',v,'out_k10');
%   
%
%  REFERENCES
%  
%  SEE ALSO
%     tom_compare_deluxe,tom_calc_isosurface
%  
%     created by FB 01/24/06
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom



if (samp~=1)
    res_map_small=zeros(size(tom_bin(res_map)));
    res_map_small=res_map(1:samp:end,1:samp:end,1:samp:end);
    res_map=tom_rescale3d(res_map_small,size(res_map),'bilinear');
end;
 tom_emwrite([output_folder '/vol_rescale.em'],res_map.*vol);

[mean, max, min, std, variance] = tom_dev(res_map,'noinfo');

if (strcmp(flag,'mygroups'))
    for i=1:size(param,1)
        vol_bin=((res_map >= param(i,1)) .*  (res_map < param(i,2)));
        vol_bin=vol_bin.*vol;
        tom_emwrite([output_folder '/vol_' num2str(param(i,1)) '_' num2str(param(i,2)) '.em'],vol_bin);
    end;
end;

figure; isosurface(vol_bin,0.5,res_map.*vol); axis image;
