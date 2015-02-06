function tom_av2_remove_feat_picklist(picklist,find_what,replace_with,mask,align_flag,rpl_flag)
%TOM_AV2_REMOVE_FEAT_PICKLIST 
%   
%
%  tom_av2_remove_feat_picklist(picklist,find_what,replace_with,mask,flag)
%
%  tom_av2_remove_feat_picklist remove features at position of the picklist
%  
%
%PARAMETERS
%
%  INPUT
%  picklist         align2d struct with pick coord
%  find_what        rplace param for output files  
%  replace_with     replace param for output files
%  mask             mask for replacement 
%  align_flag       (not implemented) align mask according param in alg list
%  rpl_flag         rpl_flag only rand_stream implemented
%
%
% OUTPUT
%
%EXAMPLE
%   matlabpool open local;  
%   tom_av2_remove_feat_picklist('pl_tops.mat','/high/','/high_rmf/',mask,0,'rand_stream')
%
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_em_classify3d
%
%   created by FB 08/09/09
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

warning off;


disp(['Reading ' picklist]);
load(picklist)


all_files=unique({align2d(1,:).filename});
disp([num2str(length(all_files)) ' images found ']);


parfor i=1:length(all_files)
    do_remove(all_files{i},align2d,find_what,replace_with,mask,align_flag,rpl_flag);
end;


function do_remove(f_name,align2d,find_what,replace_with,mask,align_flag,rpl_flag)

tic;

img=tom_emreadc(f_name);
img=img.Value;
sz_m=size(mask);
idx_mask=find(mask>0);
r_stat_mask=round(sz_m(1)*1.5);


p_idx=find(ismember({align2d(1,:).filename},f_name));


figure;

img_clean=img;
for i=1:length(p_idx)
    mask_stat=tom_spheremask(ones(size(img)),r_stat_mask,0,[align2d(1,p_idx(i)).position.x align2d(1,p_idx(i)).position.y 1]);
    idx=find(mask_stat>0);
    ime=mean2(img(idx));
    istd=std2(img(idx));
    pos=[align2d(1,p_idx(i)).position.x align2d(1,p_idx(i)).position.y]-round(sz_m./2)-[1 1];
    try
        img_cut=tom_cut_out(img,pos,sz_m);
    catch ME
        disp(['cannot cut pos: ' num2str(pos)]);
    end;
    
    if (strcmp(rpl_flag,'rand_stream'))
        rnd_str=tom_norm(rand(length(idx_mask),1),'mean0+1std');
        rnd_str=(rnd_str.*istd)+ime;
    end;
    
    img_cut(idx_mask)=rnd_str;
    img_clean=tom_paste(img_clean,img_cut,pos);
    %img_perm=tom_permute_bg(img.Value,m_tmp==0);
   
end;


new_name=strrep(f_name,find_what,replace_with);

if (strcmp(f_name,new_name)==1)
    error('new and old filename are the sampe check replacement');
end;


tom_emwrite(new_name,img_clean);

disp([f_name ' ==> ' new_name ' ' num2str(length(p_idx)) ' features removed ']);
toc;









