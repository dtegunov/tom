function tom_av2_xmipp_align_stack(doc_current_ang_filename,output_stack_name,output_struct_name,binning,filter_k,post_alg_flag)
%TOM_AV2_XMIPP_ALIGN_STACk aligns particle stack according 2 a xmipp doc
%   file
%
%  tom_av2_xmipp_align_stack(doc_current_ang_filename,output_stack_name,output_struct_name,binning)
%
%  TOM_AV2_XMIPP_ALIGN_STACK aligns a particle stack according 2 a xmipp doc
%                  in proj-match: Iter_xx_current_angles.doc
%                  in ml3d      : ml3d_it0000xx.doc
%  
%
%PARAMETERS
%
%  INPUT
%   doc_current_ang_filename         *.doc filename use abs filename (... for further processing)
%   output_stack_name                name of output em stack
%   output_struct_name               name of the output struct name 
%   binning                          (0) binning of the images 
%   filter_k                         (0) filter kernel (for xxl datasets 
%                                                   ...its faster to pre filter the stack and switch off the filter in tom_av2_em_classify3d)   
%   post_alg_flag                    (0)flag for post centering (usful for C2)                     
%                                    
%
%  OUTPUT
%
%EXAMPLE
%     
%  tom_av2_xmipp_align_stack('/fs/scratch/fbeck/test_var/dros_classic/Iter_10_current_angles.doc','st_out.em','st_out.mat',0);
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

if (nargin < 4)
    binning=0;
end;

if (nargin < 5)
    filter_k=0;
end;

if (nargin < 6)
    post_alg_flag=0;
end;


max_wrirte_chunk=250000;
%max_wrirte_chunk=2000;

[a b c]=fileparts(doc_current_ang_filename);

if (isempty(a) || strcmp(doc_current_ang_filename(1),'/')==0 )
    answ=questdlg('Warning rel doc path used (use abs path!) Continue anyway ?');
    if (strcmp(answ,'Cancel') || strcmp(answ,'No'))
        return;
    else
        disp(['Warning rel doc path: ' doc_current_ang_filename]);
    end;
end;


fprintf('%s ', 'reading currant ang doc file...');
st=tom_xmippdocread(doc_current_ang_filename);
fprintf('%s \n', ['...done! ' ]); 
doc=st;

[stack_alg st]=do_align(st,binning,doc_current_ang_filename,filter_k);


%tom_emwritec(output_stack_name,stack_alg);

if (size(stack_alg,3) < max_wrirte_chunk)
    st.split_part_stack{1}=output_stack_name;
    tom_emwritec(output_stack_name,stack_alg,'standard','single');
else
    disp(['lenght stack: ' num2str(size(stack_alg,3)) ' > ' num2str(max_wrirte_chunk) ' ==> splitting stack'  ])
    [f_base f_name f_ext]=fileparts(output_stack_name);
    packages=tom_calc_packages(ceil(size(stack_alg,3)./max_wrirte_chunk),size(stack_alg,3));
    for ip=1:size(packages,1)
        if (isempty(f_base))
            f_base=pwd();
        end;
        tmp_name=[f_base '/' f_name '_' num2str(ip) f_ext];
        st.split_part_stack{ip}=tmp_name;
        disp(['   stack alg ==> ' num2str(packages(ip,1)) ' to ' num2str(packages(ip,2)) ' ==> ' tmp_name] );
        tom_emwritec(tmp_name,stack_alg(:,:,packages(ip,1): packages(ip,2)),'standard','single');
    end;
end;

save(output_struct_name,'st');


if (post_alg_flag)
    vol=build_models(st,st.euler_ang_zxz_proj,stack_alg);
    shifts=classify(stack_alg,st,vol);
    tom_dev(shifts);
    for i=1:length(doc)
        doc(i).xoff=doc(i).xoff-shifts(i,1);
        doc(i).yoff=doc(i).yoff+shifts(i,2);
    end;
    [stack_alg st]=do_align(doc,binning,doc_current_ang_filename,filter_k);

    save(output_struct_name,'st');
   
    tom_emwritec(output_stack_name,stack_alg,'standard','single');
    
end 





function [stack_alg st]=do_align(st,binning,doc_current_ang_filename,filter_k)


fprintf('%s ', ['Converting Angles: ' ]);

num_of_entries=size(st,1);

ref_nr=zeros(num_of_entries,1);
euler_ang_zxz=zeros(num_of_entries,3);
euler_ang_zyz=zeros(num_of_entries,3);
euler_ang_zxz_proj=zeros(max([st(:).ref]),3);
euler_ang_zyz_proj=zeros(max([st(:).ref]),3);


for i=1:num_of_entries
    ref_nr(i)=st(i).ref;
    [xx,angles] = tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, st(i).psi);
    euler_ang_zxz(i,:)=angles;
    euler_ang_zyz(i,:)=[st(i).rot st(i).tilt st(i).psi];
    [aa tttemp]=tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, 0);
    euler_ang_zxz_proj(st(i).ref,:)=tttemp;
    euler_ang_zyz_proj(st(i).ref,:)=[st(i).rot st(i).tilt 0];
    if (mod(i,5000)==0)
        fprintf('%s','.');
    end;
end;
fprintf('%s \n','done!');
%tmp=euler_ang_zxz_proj(1:max(ref_nr));
%euler_ang_zxz_pro=tmp;




fprintf('%s ', ['Allocating Memory: ' ]);
try
    im_tmp=tom_spiderread(['./' st(1).name]);
    path_flag='rel_man';
catch ME
    im_tmp=tom_spiderread([st(1).name]);
    path_flag='abs_man';
end
sz_org=size(im_tmp.Value);

if (max(size(binning))==1 )
    im_tmp.Value=tom_bin(im_tmp.Value,binning);
else
    im_tmp.Value=imresize(im_tmp.Value,binning);
end;

stack_alg=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),num_of_entries,'single');
fprintf('%s \n','done!');


error_idx=zeros(num_of_entries,1);

fprintf('%s ', ['Aligning Stack: ' ]); 


parfor i=1:num_of_entries
    if (strcmp(path_flag,'rel_man'))
        tmp_name=['./' st(i).name];
    else
        tmp_name=[st(i).name];
    end;
    part_names{i}=tmp_name;
    
    try
        im=tom_spiderread(tmp_name);
    catch ME
        disp(['Error: cannot read file: ' tmp_name ' ' ME.message ]);
        im.Value=rand(sz_org(1),sz_org(2));
        error_idx(i)=1;
    end;
    
    if (max(size(binning))==1 )
        im.Value=tom_bin(im.Value,binning);
    else
        im.Value=imresize(im.Value,binning);
    end;
    
    im.Value=tom_norm(im.Value,'mean0+1std');
    
    if (max(size(binning))==1 )
        tmp_sh=[st(i).xoff st(i).yoff]./(2^binning);
    else
        tmp_sh=[st(i).xoff st(i).yoff]./(sz_org(1)./binning(1));
    end;
    
    
    if (st(i).flip==0)
        im_tmp_alg=tom_rotate(tom_shift(im.Value,tmp_sh), euler_ang_zyz(i,3));
    else
        im.Value=tom_mirror(im.Value,'x');
        im.Value=im.Value;
        im_tmp_alg=tom_rotate(tom_shift(im.Value,[-tmp_sh(1) tmp_sh(2)]), euler_ang_zyz(i,3));
    end;
    
    if (filter_k > 0)
        stack_alg(:,:,i)=single(tom_filter(im_tmp_alg,filter_k));
    else
        stack_alg(:,:,i)=single(im_tmp_alg);
    end;
    
    if (mod(i,5000)==0)
        fprintf('%s','.');
    end;
end;

fprintf('%s \n', ['...done! ' ]); 

disp([num2str(length(find(error_idx==1)) ) ' errors in alignment!']);

clear('st');

st.ref_nr=ref_nr;
st.euler_ang_zxz=euler_ang_zxz;
st.euler_ang_zyz=euler_ang_zyz;
st.euler_ang_zxz_proj=euler_ang_zxz_proj;
st.error_idx=error_idx;
st.part_names=part_names;
st.doc_name=doc_current_ang_filename;
st.euler_ang_zyz_proj=euler_ang_zyz_proj;
st.pre_filter_k=filter_k;

function vol=build_models(part_st,proj_angles,part_stack)



classes=build_classes(part_st,part_stack,length(proj_angles));


vol=backproj(classes,proj_angles,'C1');
    
   


function classes=build_classes(part_st,stack_alg,l_ang)

ref_nr=part_st.ref_nr;
classes=zeros(size(stack_alg,1),size(stack_alg,2),l_ang);
part_per_class=zeros(l_ang,1);

for i=1:size(stack_alg,3)
    classes(:,:,ref_nr(i))=classes(:,:,ref_nr(i))+ stack_alg(:,:,i);
    part_per_class(ref_nr(i))=part_per_class(ref_nr(i))+1;
end;

for i=1:size(classes,3)
   if part_per_class(i) > 0
       classes(:,:,i)=tom_norm(classes(:,:,i)./part_per_class(i),'mean0+1std');
   end;
end;




function vol=backproj(classes,proj_angles,sym)



numc=1;
sym_angs(1)=0;
if (strcmp(sym,'C1')==0)
    numc=str2double(sym(2));
    sym_angs(1)=0;
    for i=2:numc;
        sym_angs(i)=360./numc;
    end;
end;



zz=1;
tmpp=zeros(size(classes,3),3);
for i=1:size(classes,3)
    if (std2(classes(:,:,i) )~=0 )
        for ii=1:numc
            tmpp(zz,:)=[(proj_angles(i,1)+sym_angs(ii)) proj_angles(i,2) proj_angles(i,3)];
            zz=zz+1;
        end;
    end;
end;
proj_angles_clean=tmpp(1:(zz-1),:);
clear('tmpp');



THICK=size(classes,1);

%backproject image stack
vol=zeros(size(classes,1),size(classes,1),size(classes,1),'single');
for i=1:size(proj_angles,1)
    if (std2(classes(:,:,i) )~=0 )
        
        for ii=1:numc
            w_func=tom_calc_weight_functionc([size(vol,1) size(vol,2)],proj_angles_clean,THICK,[(proj_angles(i,1) +sym_angs(ii)) proj_angles(i,2) proj_angles(i,3)]);
            w_proj=tom_apply_weight_function(classes(:,:,i),w_func);
            w_proj=single(w_proj);
            tom_backproj3d_euler(vol,w_proj,proj_angles(i,1)+sym_angs(ii),proj_angles(i,2),proj_angles(i,3),[0 0 0]);
            if isempty(find(isnan(vol)))==0
                disp([num2str(length(isnan(vol))) 'NAN found']);
                vol=tom_remove_nan(vol);
                vol=single(vol);
                tom_dev(vol);
            end
        end;
        
    end;
    
end;


function [shifts]=classify(stack_alg_all,part_st,vol)


ref=zeros(size(vol,1),size(vol,2),size(part_st.euler_ang_zxz_proj,1));
disp('Projecting ')

for i=1:size(part_st.euler_ang_zxz_proj,1)
    tmp_vol=tom_rotate(vol,part_st.euler_ang_zxz_proj(i,:));
    ref(:,:,i)=sum(tmp_vol,3);
end;


shifts=zeros(size(stack_alg_all,3),2);
fprintf('%s ', 'Classifying: ' );

num_of_nodes=floor(size(stack_alg_all,3)./12000);

%num_of_nodes=floor(size(stack_alg_all,3)./150);
if (num_of_nodes==0)
    num_of_nodes=1;
end;

packages=tom_calc_packages(num_of_nodes,size(stack_alg_all,3));

mask_part=tom_spheremask(ones(size(vol,1),size(vol,2)),size(vol,2)-3,1);

mask_ccf=tom_spheremask(ones(size(vol,1),size(vol,2)),2,3);

mid_img=[floor(size(vol,1)./2)+1 floor(size(vol,2)./2)+1];

for ii=1:size(packages,1)
    
    offs=packages(ii,1)-1;
    stack_alg=stack_alg_all(:,:,(packages(ii,1):packages(ii,2)));
    tic;
    parfor i=1:size(stack_alg,3)  %...parfor it
        try
            ref_nr=part_st.ref_nr(i+offs);
            ref_tmp=squeeze(ref(:,:,ref_nr,1:end));
            stack_tmp=stack_alg(:,:,i);
            ccf=tom_corr(stack_tmp.*mask_part,ref_tmp.*mask_part,'norm',mask_part);
            [pos val]=tom_peak((ccf.*mask_ccf)+(-1+mask_ccf),'spline');
            sh_tmp=pos(1,:)-mid_img;
            shifts(i+offs,:)=sh_tmp;
         catch
            disp(['Error in cc-module: ']);
            disp(['Part Nr: ' num2str(i) ]);
            disp(['Size ref tmp ' num2str(size(ref_tmp)) 'Size stack_alg: ' num2str(size(stack_alg)) ]);
            disp(lasterr);
        end;
    end;
    
  end;




function [ref_nr cc shift_max]=stack_ccc(ref,stack,data,mask,norm_flag,filter,calc_shifts,mask_ccf)


cc=zeros(data.l_ref,1);
shift_vals=zeros(data.l_ref,2);

for ii=1:data.l_ref
    
    
    if (filter > 0)
        ref_tmpo=tom_filter(ref(:,:,ii),filter);
        part_tmpo=tom_filter(stack,filter);
    else
        ref_tmpo=ref(:,:,ii);
        part_tmpo=stack;
    end;
    
    if (strcmp(norm_flag,'norm') && calc_shifts > 0)
        ccf=tom_corr(ref_tmpo.*mask,part_tmpo.*mask,norm_flag,mask);
        [pos val]=tom_peak((ccf.*mask_ccf)+(-1+mask_ccf));
        sh_tmp=pos(1,:)-data.mid_img;
    end;
    
    if (strcmp(norm_flag,'norm') && calc_shifts==0)
        idx=find(mask==1);
        val=sum((tom_norm(ref_tmpo(idx),'mean0+1std').*tom_norm(part_tmpo(idx),'mean0+1std')))./ (length(idx));
        pos=data.mid_img;
        sh_tmp=pos(1,:)-data.mid_img;
    end;
    
    if (strcmp(norm_flag,'particle') && calc_shifts==1)
        ccf=tom_corr(ref(:,:,ii).*mask,stack.*mask);
        [pos val]=tom_peak((ccf.*mask_ccf)+(-1+mask_ccf));
        sh_tmp=pos(1,:)-data.mid_img;
    end;
    
    if (strcmp(norm_flag,'variance'))
        [a1 b c d e1]=tom_dev(ref(:,:,ii).*mask,'noinfo',mask);
        [a2 b c d e2]=tom_dev(stack.*mask,'noinfo',mask);
        ccf=(a1-a2).*(a1-a2);
        val=ccf;
        sh_tmp=[0 0];
    end;
    
    cc(ii)=val;
    shift_vals(ii,:)=sh_tmp;
    
end;


[max_val pos_max]=max(cc);
shift_max=shift_vals(pos_max,:);
ref_nr=pos_max;








