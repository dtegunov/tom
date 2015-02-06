function tom_av2_em_refine3d(part_stack_f,part_st_f,vol_mask,num_of_classes,num_of_iter,outputdir,filter,sym,restart_flag,norm_flag,max_sh,max_rot)
%tom_av2_em_refine3d classifys a particle stack into n classes in 3d
%  
%
%PARAMETERS
%
%  INPUT
%   part_stack_f                     filename of the aligned part stack  (...from tom_av2_xmipp_align_stack)
%   part_st_f                        filename of the output structure    (...from tom_av2_xmipp_align_stack)
%   vol_mask                         3d mask for focusing classification
%   num_of_classes                   number of classes
%   num_of_iter                      number of iterations  
%   outputdir                        output directory (use abs path for further processing !!)  
%   filter                           (0) filterkernel
%   sym                              (C1) symmetry flag ...only C2 implemented so far!
%                                         for using C2 the dataset must have been processed C2      
%   restart_flag                     (0) restart flag to restart the classification   
%   norm_flag                        (norm) flag for norming 
%   max_sh                           (0) max_shift  
%
%  OUTPUT
%
%EXAMPLE
%     
% tom_av2_em_classify3d('st_out.em','st_out.mat','mask.em',3,50,'/fs/scratch/fbeck/test_var/out/');
%
% %shifts
% tom_av2_em_classify3d('st_out.em','st_out.mat','mask.em',3,50,'/fs/scratch/fbeck/test_var/out/',0,'C1',0,'norm',5); 
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_xmipp_align_stack, tom_av2_em_classify3d_2xmipp
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


part_st=load(part_st_f);
part_st=part_st.st;

if (nargin < 7)
    filter=0;
end;

if (nargin < 8)
    sym='C1';
end;

if (nargin < 9)
    restart_flag=0;
end;

if (nargin < 10)
    norm_flag='norm';
end;

if (nargin < 11)
    max_sh=0;
end;

if (nargin < 11)
    max_rot=0;
end;


xmipp_output_flag=1;

if (ischar(vol_mask))
    vol_mask=tom_emread(vol_mask);
    vol_mask=vol_mask.Value;
end;

doc_file=part_st.doc_name;
part_stack=tom_emread(part_stack_f);
part_stack=part_stack.Value;

if (isempty(vol_mask))
    vol_mask=ones(size(part_stack,1),size(part_stack,1),size(part_stack,1));
end;


[a b c]=fileparts(outputdir);

if (isempty(a) || strcmp(outputdir(1),'/')==0 )
    answ=questdlg('Warning rel output path used (use abs path!) Continue anyway ?');
    if (strcmp(answ,'Cancel') || strcmp(answ,'No'))
        return;
    else
        disp(['Warning rel doc path: ' outputdir]);
    end;
end;


%load spider file and transform angles from zxz 2 zyz 
%proj_angles=transform_angles(porjangles_spifile);

proj_angles=part_st.euler_ang_zxz_proj;

if (restart_flag==0)
     %split in random subsets ...ini
    part_st=build_subsets(part_st,num_of_classes,'by_order',1);
    iter_start=1;
else
    iter_start=last_iteration(num_of_classes,outputdir);
    disp(['Restarting with iteration Nr: ' num2str(iter_start) ]);
end;


disp('Building mask projections:');
build_mask_projections(vol_mask,proj_angles,outputdir,'beckster');



%main loop rock Expectation Max
for iter_count=iter_start:num_of_iter
     
     if (restart_flag==0)
        %build models ..backporj
        build_models(part_st,proj_angles,num_of_classes,part_stack,outputdir,iter_count,sym);
     end;   
     
     %project
     build_projections(proj_angles,'beckster',num_of_classes,vol_mask,outputdir,iter_count);    
    
     %classify
     [part_st cc_mean class_count movin_part shift_stat]=classify(part_stack,part_st,num_of_classes,outputdir,iter_count,norm_flag,filter,max_sh,max_rot);
    
     disp(['Iteration Nr: ' num2str(iter_count) '  Mean cc is: ' num2str(cc_mean)  ' Classes Nr. ' num2str(class_count') ' Movin Part: ' num2str(movin_part)  ' Sum Shifts: '  num2str(shift_stat.sum) ' Mean Shifts: '  num2str(shift_stat.mean)  ' Max shifts: '  num2str(shift_stat.max)]);
     restart_flag=0;
     stat(iter_count).cc_mean=cc_mean;
     stat(iter_count).class_count=class_count;
     stat(iter_count).movin_part=movin_part;
     stat(iter_count).mean_shift=shift_stat.mean;
    
     part_st.doc_name=doc_file;
     part_st.outputdir=outputdir;
     save([outputdir '/stat.mat'],'stat');
     save([outputdir '/part_st.mat'],'part_st');
end;


if (xmipp_output_flag==1)
    disp('Creating new sel and doc files !');
    warning off;
    mkdir([outputdir '/xmipp']);
    warning on;
    for i=1:num_of_classes
        tom_av2_em_classify3d_2xmipp([outputdir '/part_st.mat'],i,[outputdir '/xmipp/cl_' num2str(i) '.sel'],[outputdir '/xmipp/cl_' num2str(i) '.doc']);
    end;
end;



function part_st=build_subsets(part_st,num_of_classes,split_flag,iter_num)

num_of_part=length(part_st.ref_nr);

if (strcmp(split_flag,'by_order'))
    for i=1:num_of_part
        part_st.class(iter_num,i)=mod(i,num_of_classes)+1;
    end;
end;


function build_models(part_st,proj_angles,num_of_classes,part_stack,outputdir,iter_num,sym)

warning off;

%if (exist([outputdir '/models' ],'dir')==0)
    mkdir([outputdir '/models' ]);
    unix(['chmod ugo+rwx ' outputdir '/models']);
%end;

%if (exist([outputdir '/classes' ],'dir')==0)
    mkdir([outputdir '/classes' ]);
    unix(['chmod ugo+rwx ' outputdir '/classes']);
%end;

warning on;

%allocate memory
classes=zeros(size(part_stack,1),size(part_stack,2),length(proj_angles),num_of_classes);
%vol=zeros(size(part_stack,1),size(part_stack,1),size(part_stack,1));

num_of_nodes=floor(size(part_stack,3)./25000);

num_of_nodes=1;

if (num_of_nodes==0)
    num_of_nodes=1;
end;

packages=tom_calc_packages(num_of_nodes,size(part_stack,3));

fprintf('%s ', 'Building classes: ' );

for ii=1:size(packages,1)
    tmp_st.ref_nr=part_st.ref_nr(packages(ii,1):packages(ii,2));
    tmp_st.class=part_st.class(:,packages(ii,1):packages(ii,2));
    tmp_st_part=part_stack(:,:,packages(ii,1):packages(ii,2));
    if (isfield(part_st,'shifts')==1)
        tmp_st.shifts=part_st.shifts;
    end;
    
%     for model_nr=1:num_of_classes
%         classes_tmp(:,:,:,model_nr)=build_classes(tmp_st,tmp_st_part,model_nr,iter_num,length(proj_angles));
%     end;
    
    for model_nr=1:num_of_classes
        classes(:,:,:,model_nr)=build_classes(tmp_st,tmp_st_part,model_nr,iter_num,length(proj_angles));
    end;
    fprintf('%s','.');
end;
clear('classes_tmp');

for model_nr=1:num_of_classes
    %write out model with iter and mode post-fixluja
    tom_emwrite([outputdir '/classes/classes_iter_' num2str(iter_num) '_model_' num2str(model_nr) '.em'],classes(:,:,:,model_nr));
    unix(['chmod ugo+rwx ' outputdir '/classes/classes_iter_' num2str(iter_num) '_model_' num2str(model_nr) '.em']);
end;
 
fprintf('%s \n', '...done ' );

parfor model_nr=1:num_of_classes
% for model_nr=1:num_of_classes
    vol=backproj(classes(:,:,:,model_nr),proj_angles,sym);
    %write out model with iter and mode post-fixluja
    tom_emwrite([outputdir '/models/model_iter_' num2str(iter_num) '_model_' num2str(model_nr) '.em'],tom_norm(vol,'mean0+1std'));
    unix(['chmod ugo+rwx ' outputdir '/models/model_iter_' num2str(iter_num) '_model_' num2str(model_nr) '.em']);
end;



        
%old paraellization works only for small stacks

% parfor model_nr=1:num_of_classes
% % for model_nr=1:num_of_classes
%     classes=build_classes(part_st,part_stack,model_nr,iter_num,length(proj_angles));
%     %write out model with iter and mode post-fixluja
%     tom_emwrite([outputdir '/classes/classes_iter_' num2str(iter_num) '_model_' num2str(model_nr) '.em'],classes);
%     unix(['chmod ugo+rwx ' outputdir '/classes/classes_iter_' num2str(iter_num) '_model_' num2str(model_nr) '.em']);
%     vol=backproj(classes,proj_angles,sym);
%     %write out model with iter and mode post-fixluja
%     tom_emwrite([outputdir '/models/model_iter_' num2str(iter_num) '_model_' num2str(model_nr) '.em'],vol);
%     unix(['chmod ugo+rwx ' outputdir '/models/model_iter_' num2str(iter_num) '_model_' num2str(model_nr) '.em']);
% end;




function proj=build_projections(proj_angles,proj_flag,num_of_classes,vol_mask,outputdir,iter_num)

if (exist([outputdir '/projections' ],'dir')==0)
    mkdir([outputdir '/projections' ]);
    unix(['chmod ugo+rwx ' outputdir '/projections' ]);
end;

tom_emwrite([outputdir '/mask.em'],vol_mask);

if (strcmp(proj_flag,'beckster'))
    for model_nr=1:num_of_classes
        vol=tom_emread([outputdir '/models/model_iter_' num2str(iter_num) '_model_' num2str(model_nr) '.em']);
        %vol=vol.Value.*vol_mask;
        vol=vol.Value;
        proj=project(vol,proj_angles);
        tom_emwrite([outputdir '/projections/proj_iter_' num2str(iter_num) '_model_' num2str(model_nr) '.em'],proj);
    end;
end;

if (strcmp(proj_flag,'spider'))
    for i=1:num_of_classes
        tom_spider_proj('.','vol.spi','../../../241108/spider/logs/refangles1.spi','projstack_test.spi');
    end;
end;


function [part_st mean_cc class_count movin_part shift_stat]=classify(stack_alg_all,part_st,num_of_classes,outputdir,iter_num,norm_flag,filter,max_sh,max_rot)





%load mask_stack in memory ...to be changed
mask_stack=tom_emread([outputdir '/mask_proj/mask_proj.em']);
mask_stack=mask_stack.Value;


try 
    h=tom_reademheader([outputdir '/projections/proj_iter_' num2str(iter_num) '_model_' num2str(1) '.em']);
    ref=zeros(h.Header.Size(1),h.Header.Size(2),h.Header.Size(3),num_of_classes);
catch ME
    disp('Allocation Ref Stack Error in Classification Module!');
    error(ME.message);
end;

%load references in memory ...to be changed
for model_nr=1:num_of_classes
    tmp=tom_emread([outputdir '/projections/proj_iter_' num2str(iter_num) '_model_' num2str(model_nr) '.em']);
    tmp=tmp.Value;
    ref(:,:,:,model_nr)=tmp;
end;
clear('tmp');


new_ref_nr=zeros(size(stack_alg_all,3),1);
class_count=zeros(num_of_classes,1);
cc_max=zeros(size(stack_alg_all,3),1);
shifts=zeros(size(stack_alg_all,3),2);
fprintf('%s ', 'Classifying: ' );

num_of_nodes=floor(size(stack_alg_all,3)./6000);

%num_of_nodes=floor(size(stack_alg_all,3)./150);
if (num_of_nodes==0)
    num_of_nodes=1;
end;

packages=tom_calc_packages(num_of_nodes,size(stack_alg_all,3));


data.l_ref=num_of_classes;
data.mid_img=[floor(h.Header.Size(1)./2)+1 floor(h.Header.Size(2)./2)+1];

if (max_sh > 0)
    mask_ccf=tom_spheremask(ones(h.Header.Size(1),h.Header.Size(2)),max_sh,1);
else
    mask_ccf=ones(h.Header.Size(1),h.Header.Size(2));
end;

polar_size=size(tom_cart2polar(ones(h.Header.Size(1),h.Header.Size(2))));
if (max_rot > 0)
    sz_y=(max(polar_size)./360).*max_rot;
    mask_rot=tom_rectanglemask(zeros(polar_size),[polar_size(1)-1 sz_y]);
    mask_part=tom_spheremask(ones(h.Header.Size(1),h.Header.Size(2)),round(h.Header.Size(1)./2)-1);
else
    mask_rot=ones(polar_size);
    mask_part=ones(h.Header.Size(1),h.Header.Size(2));
end;

for ii=1:size(packages,1)
    
    offs=packages(ii,1)-1;
    stack_alg=stack_alg_all(:,:,(packages(ii,1):packages(ii,2)));
    tic;
    %for i=1:size(stack_alg,3)
    parfor i=1:size(stack_alg,3)  %...parfor it
        try
            ref_nr=part_st.ref_nr(i+offs);
            ref_tmp=squeeze(ref(:,:,ref_nr,1:end));
            stack_tmp=stack_alg(:,:,i);
            mask_tmp=mask_stack(:,:,ref_nr);
            [ref_nr cc_vals shift_max]=stack_ccc(ref_tmp,stack_tmp,data,mask_tmp,norm_flag,filter,max_sh,mask_ccf,mask_rot,mask_part);
            new_ref_nr(i+offs)=ref_nr;
            cc_max(i+offs)=cc_vals(ref_nr);
            shifts(i+offs,:)=shift_max;
        catch ME
            disp(['Error in cc-module: ']);
            disp(['Part Nr: ' num2str(i) ]);
            disp(['Size ref tmp ' num2str(size(ref_tmp)) 'Size stack_alg: ' num2str(size(stack_alg)) ]);
            disp(lasterr);
        end;
    end;
    fprintf('%s','.');
    toc;
    
end;

fprintf('%s \n', ['done! ' ]); 
mean_cc=mean(cc_max);
shift_stat.mean=mean(shifts,1);
shift_stat.max=max(abs(shifts));
shift_stat.std=std(shifts,1);
shift_stat.sum=sum(abs(shifts),1);

part_st.class(iter_num+1,:)=new_ref_nr;
part_st.shifts(iter_num+1,:,:)=shifts;


%build part stat
for i=1:num_of_classes
    idx=find(part_st.class(iter_num+1,:)==i);
    class_count(i)=length(idx);
end;

try
  movin_part=length(find( (part_st.class(iter_num,:)-part_st.class(iter_num+1,:))>0));  
catch ME
    disp('warning; n-1 mod nr missing!');
end;




function classes=build_classes(part_st,stack_alg,model_nr,iter_num,l_ang)

ref_nr=part_st.ref_nr;
classes=zeros(size(stack_alg,1),size(stack_alg,2),l_ang);
part_per_class=zeros(l_ang,1);

for i=1:size(stack_alg,3)
    if (part_st.class(iter_num,i)==model_nr)
        %if (part_st.all_dists(i) < 4)
            if (isfield(part_st,'shifts'))
                classes(:,:,ref_nr(i))=classes(:,:,ref_nr(i))+ tom_shift(stack_alg(:,:,i),squeeze(-[part_st.shifts(iter_num,i,:)])');
               % disp('all shifts');
            else
                classes(:,:,ref_nr(i))=classes(:,:,ref_nr(i))+ stack_alg(:,:,i);
            end;
            part_per_class(ref_nr(i))=part_per_class(ref_nr(i))+1;

        %end;
    end;
end;

%fprintf('%s \n', '...done! '); 

%fprintf('%s ', 'Norming classes: ' ); 
for i=1:size(classes,3)
   if part_per_class(i) > 0
       classes(:,:,i)=tom_norm(classes(:,:,i)./part_per_class(i),'mean0+1std');
   end;
    
%     if (mod(i,20)==0)    
%         fprintf('%s','.');
%     end;
end;
%fprintf('%s \n', '...done! '); 

function vol=backproj(classes,proj_angles,sym)

fprintf('%s ', ['Reconstructing volume: ' ]); 

numc=1;
sym_angs(1)=0;
if (strcmp(sym,'C1')==0)
    numc=str2double(sym(2));
    sym_angs(1)=0;
    for i=2:numc;
        sym_angs(i)=360./numc;
    end;
end;




%remove empty classes
% zz=1;
% tmpp=zeros(size(classes,3),3);
% for i=1:size(classes,3)
%     if (std2(classes(:,:,i) )~=0 )
%         tmpp(zz,:)=proj_angles(i,:);
%         zz=zz+1;
%     end;
% end;
% proj_angles_clean=tmpp(1:(zz-1),:);
% clear('tmpp');

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
        
        
        if (mod(i,50)==0)
            fprintf('%s','.');
        end;
    end;
    
end;
fprintf('%s \n', ['...done! ' ]); 




function proj_angles_zxz=transform_angles(spider_angles_file)

%load anglefile
proj_ang_zyz=importdata(spider_angles_file);
proj_ang_zyz=proj_ang_zyz.data;
proj_ang_zyz=proj_ang_zyz(:,4:5);

parfor i=1:size(proj_ang_zyz,1)
    [xx,angles] = tom_eulerconvert_xmipp(proj_ang_zyz(i,2),proj_ang_zyz(i,1),0);
    proj_angles_zxz(i,:) = angles';
end;


function proj=project(vol,proj_angles)

proj=zeros(size(vol,1),size(vol,2),size(proj_angles,1));
disp('Projecting ')
parfor i=1:size(proj_angles,1)
    tmp_vol=tom_rotate(vol,proj_angles(i,:));
    proj(:,:,i)=sum(tmp_vol,3);
end;
fprintf('%s \n', ['...done! ' ]); 


function [ref_nr cc shift_max ang_max]=stack_ccc(ref,stack,data,mask,norm_flag,filter,calc_shifts,mask_ccf,mask_rot,mask_part)


cc=zeros(data.l_ref,1);
ang=zeros(data.l_ref,1);
shift_vals=zeros(data.l_ref,2);

rot_flag=std2(mask_rot);

for ii=1:data.l_ref
    
    
    if (filter > 0)
        ref_tmpo=tom_filter(ref(:,:,ii),filter);
        part_tmpo=tom_filter(stack,filter);
    else
        ref_tmpo=ref(:,:,ii);
        part_tmpo=stack;
    end;
    
    if (strcmp(norm_flag,'norm') && calc_shifts >= 0 && rot_flag>0)
        filter_st.Apply=0;
        [ang_out sh_tmp val]=tom_av2_align(ref_tmpo,part_tmpo,mask_part,mask_rot,mask_ccf,filter_st,1);
        sh_tmp=-sh_tmp(1:2);
        ang_out=-ang_out(1);
    end;
    
    if (strcmp(norm_flag,'norm') && calc_shifts > 0 && rot_flag==0)
        ccf=tom_corr(ref_tmpo.*mask,part_tmpo.*mask,norm_flag,mask);
        [pos val]=tom_peak((ccf.*mask_ccf)+(-1+mask_ccf));
        sh_tmp=pos(1,:)-data.mid_img;
    end;
    
    if (strcmp(norm_flag,'norm') && calc_shifts==0 && rot_flag==0)
        idx=find(mask==1);
        val=sum((tom_norm(ref_tmpo(idx),'mean0+1std').*tom_norm(part_tmpo(idx),'mean0+1std')))./ (length(idx));
        pos=data.mid_img;
        sh_tmp=pos(1,:)-data.mid_img;
    end;
    
    if (strcmp(norm_flag,'particle') && calc_shifts==1 && rot_flag==0)
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
    ang(ii)=ang_out;
end;


[max_val pos_max]=max(cc);
shift_max=shift_vals(pos_max,:);
ang_max=ang(pos_max);
ref_nr=pos_max;



    


function build_mask_projections(mask,proj_angles,outputdir,proj_flag)

if (exist([outputdir '/mask_proj' ],'dir')==0)
    mkdir([outputdir '/mask_proj' ]);
    unix(['chmod ugo+rwx ' outputdir '/mask_proj']);
end;

if (strcmp(proj_flag,'beckster'))
    
    if (std(std(std(mask)))==0)
        proj=ones(size(mask,1),size(mask,2),size(proj_angles,1));
    else
        proj=project(mask,proj_angles);
        for i=1:size(proj,3)
            proj(:,:,i)=proj(:,:,i) > 0; %make it binary again!
        end;
    end;
    tom_emwrite([outputdir '/mask_proj/mask_proj.em'],proj);
end;



function iter_num=last_iteration(num_of_classes,outputdir)

dd=dir([outputdir '/models/model_iter_*_model_' num2str(num_of_classes) '.em']);

all_iter=zeros(length(dd),1);
for i=1:length(dd)
  all_iter(i)=str2double(strtok(dd(i).name,'model_iter_'));
end;

iter_num=max(all_iter);



