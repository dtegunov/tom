function tom_av2_permute_bg(f_part_stack,f_part_struct,f_mask3d)
%TOM_AV2_PERMUTE_BG premutes background pixels
%   
%
%  tom_av2_permute_bg(f_part_stack,f_part_struct,f_mask3d)
%
%  TOM_AV2_PERMUTE_BG premutes background pixels
%  
%
%PARAMETERS
%
%  INPUT
%   f_part_stack           *.doc filename use abs filename (... for further processing)
%   f_part_struct          name of output em stack
%   f_mask3d               name of the output struct name 
%   
%
%  OUTPUT
%
%EXAMPLE
%     
%  tom_av2_permute_bg('st_out.em','st_out.mat','mask3d.em');
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



load(f_part_struct);

stack=tom_emread(f_part_stack);
stack=stack.Value;

mask3d=tom_emread(f_mask3d);
mask3d=mask3d.Value;
mask_rest_idx=find(tom_spheremask(ones(size(stack,1),size(stack,2)),round((size(stack,1)./2)-1))==0);

doc=tom_xmippdocread(st.doc_name);

all_cl=unique(st.ref_nr);

cl_db=zeros(64,64);

for i=1:length(all_cl)
    cl_idx=find([st.ref_nr==all_cl(i)]);
    mask=sum(tom_rotate(mask3d,st.euler_ang_zxz_proj(all_cl(i),:)),3)>0;
    try
        for ii=1:length(cl_idx)
            
            part_idx=cl_idx(ii);
            tmp_im=stack(:,:,part_idx);
            
            indd=find(((tmp_im==0)==0).*(mask==0));
            ind_rand=randperm(length(indd));
            
            cl_db=cl_db+tmp_im;
            
            tmp_im(indd)=tmp_im(indd(ind_rand));
            
            idx_z=find(tmp_im==0);
            fill=tom_norm(rand(length(idx_z),1),'mean0+1std').*std(tmp_im(mask_rest_idx)) + mean2(tmp_im(mask_rest_idx));
            tmp_im(idx_z)=fill;
            
            %rotate Back!
            tmp_sh=[-doc(part_idx).xoff -doc(part_idx).yoff];
            
            if (doc(part_idx).flip==0)
                % im_tmp_alg=tom_rotate(tom_shift(tmp_im,tmp_sh), -st.euler_ang_zyz(part_idx,3));
                im_tmp_alg=tom_shift(tom_rotate(tmp_im+100, -st.euler_ang_zyz(part_idx,3)),tmp_sh);
                
            else
                %im_tmp_alg=tom_rotate(tom_shift(tmp_im,[-tmp_sh(1) tmp_sh(2)]), -st.euler_ang_zyz(part_idx,3));
                im_tmp_alg=(tom_shift(tom_rotate(tmp_im+100,-st.euler_ang_zyz(part_idx,3)),[-tmp_sh(1) tmp_sh(2)]));
                im_tmp_alg=tom_mirror(im_tmp_alg,'x');
            end;
            [a b c]=fileparts(st.part_names{part_idx});
            new_name=[a '_perm/' b c];
            warning off;
            mkdir([a '_perm/']);
            warning on;
            idx_z=find(im_tmp_alg<0.5);
            im_tmp_alg=im_tmp_alg-100;
            fill=tom_norm(rand(length(idx_z),1),'mean0+1std').*std(tmp_im(mask_rest_idx)) + mean2(tmp_im(mask_rest_idx));
            im_tmp_alg(idx_z)=fill;
            tom_spiderwrite(new_name,im_tmp_alg);
        end;
    catch ME
        disp(['Error: ' st.part_names{i}])
        disp(ME.message);
    end;
    %debug it baby
    cl_db=zeros(64,64);
    disp([num2str(i) ' of ' num2str(length(all_cl)) ' done!' ]);
end;

