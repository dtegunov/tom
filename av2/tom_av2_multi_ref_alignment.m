function class_st=tom_av2_multi_ref_alignment(path_stack,path_ref_stack,filter_param,correction_flag,alg_iter,package,task_number,disp_flag,disp_info)
%TOM_AV2_MULTI_REF_ALILGNMENT performs iterative rotational, translational alignment of an image stack (im) relative to reference stack (ref)
%
%   class_st=tom_av2_multi_ref_alignment(path_stack,path_ref_stack,filter_param,correction_flag,alg_iter,package,task_number,disp_flag,disp_info)
%
%PARAMETERS
%
%  INPUT
%   path_stack          ...
%   path_ref_stack      ...
%   filter_param        contains filter and mask structures:
%                                filter_param.mask.classify1 
%                                filter_param.mask.classify2
%                                filter_param.mask.align 
%                                filter_param.mask.ccf_rot
%                                filter_param.mask.ccf_trans
%                                filter_param.filter.classify
%                                filter_param.filter.align
%                                                
%                                filter_st -structure
%                                filter_st.Apply:  1 apply filter, 2 use default values, 0 not
%                                filter_st.Value:  vector with parameters [low, high, smooth ...]
%                                filter_st.Method: i.e. 'circ', 'quadr', 'bandpass'
%                                filter_st.Space: 'real' or 'fourier'
%                                filter_st.Times: 'apply filter n-times                                  
%
%                                mask_st -structure
%                                mask_st.Apply:  1 built mask according to the struct values, 2 use default values, 0 create ones
%                                mask_st.Value:  vector with parameters [size,radius,smooth,center ...]
%                                mask_st.Method: 'sphere' 'sphere3d' 'cylinder3d' 'rectangle'...
%   correction_flag     'no_correction' no alignment done just classification
%                       'post_correction' alignment done after classification 
%                       'pre_correction' alignment before
%                       classification
%   alg_iter            number of iterations for the alignment step
%   package             [start stop] of the particle in the image stack ...use '' for whole stack 
%   task_number         number of the task ...needed for more than one CU ...use '' for one                                     disp_flag_alg
%   disp_flag           flag (1: show multiref, 2: show alignment  3: show multiref+alignment  0: off) for demo mode via graphical interface
%                                 works only with on CU 
%   disp_info           ...
%  
%  OUTPUT
%   class_st            class_st(1,part_num):  class according to the ref stack;
%                       class_st(2,part_num):  value of the ccf maximum
%                       class_st(3,part_num):  shift x 
%                       class_st(4,part_num):  shift y
%                       class_st(5,part_num):  rotation angle
%
%EXAMPLE
%   ... = tom_av2_multi_ref_alignment(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB 01/24/06
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


root_path_stack=fileparts(path_stack);
if (isempty(root_path_stack)==1)
    if (task_number > 1)
        error('stack path must be absolute!');
    end;
end;

root_path_ref=fileparts(path_ref_stack);
if (isempty(root_path_ref)==1)
    if (task_number > 1)
        error('ref path must be absolute!');
    end;
end;


ref_head=tom_reademheader(path_ref_stack); ref_head=ref_head.Header; 
stack_head=tom_reademheader(path_stack); stack_head=stack_head.Header;

 
if (nargin==2)
   package=[1 stack_head.Size(3)];
   task_number=0; disp_flag=0; 
   filter_param=tom_av2_build_filter_param('',[Header.Size(1) Header.Size(2)],'default','tom_av2_multi_ref_alignment');
   rotation_flag='no_correction'; 
   alg_iter=2;
end;

if (nargin==3)
    task_number=0; disp_flag=0; 
    package=[1 stack_head.Size(3)];
    rotation_flag='no_correction'; 
    alg_iter=2;
end;

if (nargin==4)
     disp_flag=0; task_number=0;
     package=[1 stack_head.Size(3)];
     alg_iter=2;
end;

if (nargin==5)
    disp_flag=0; disp_info=zeros(ref_head.Size(3)+1); 
    package=[1 stack_head.Size(3)];
    task_number=0;
end;
 
if (nargin==6)
    disp_info=zeros(ref_head.Size(3)+1); 
    task_number=0;
    disp_flag=0;
end;

if (nargin==7)
    disp_info=zeros(ref_head.Size(3)+1); 
    disp_flag=0;
end;

if (nargin==8)
    disp_info=zeros(ref_head.Size(3)+1); 
end;


if (isempty(filter_param)==1)
    filter_param=tom_av2_build_filter_param('',[stack_head.Size(1) stack_head.Size(2)],'default','tom_av2_multi_ref_alignment');
end;

%check for empty structure
filter_param=tom_av2_build_filter_param(filter_param,[stack_head.Size(1) stack_head.Size(2)],'default','tom_av2_multi_ref_alignment');

if (isempty(correction_flag))
    correction_flag='no_correction';
end;

if (isempty(alg_iter))
    alg_iter=2;
end;

if (isempty(package))
   package=[1 stack_head.Size(3)];
end;

if (isempty(task_number))
    task_number=0;
end;

if (isempty(disp_flag))
    disp_flag=0;
end;

if (isempty(disp_info))
    disp_info=zeros(ref_head.Size(3)+1);
end;


% transfer often used data
rot_angle=0;
size_stack=stack_head.Size;
size_ref_stack=ref_head.Size;
middle_part(1)=(size_stack(1)./2)+1;
middle_part(2)=(size_stack(2)./2)+1;


% create masks
mask_classify1=tom_create_mask(filter_param.mask.classify1);
mask_classify2=tom_create_mask(filter_param.mask.classify2);
mask_align=tom_create_mask(filter_param.mask.align);
mask_ccf_rot=tom_create_mask(filter_param.mask.cc_rot);
mask_ccf_trans=tom_create_mask(filter_param.mask.cc_trans);


%command window print
store.end=size_stack(3); store.num_of_mesure=1;  store.mesured=0;
[store]=tom_disp_estimated_time(store,'start','classify');

%display   filter_param.filter
disp_classify(1,1,1,1,1,'new',1,disp_flag); 

all_ref=tom_emread(path_ref_stack);
all_part=tom_emread(path_stack);

zz=1;
for i=package(1):package(2) % loop over all particles

    %read the particle of stack
    %part=tom_emread([path_stack],'subregion',[1 1 i],[size_stack(1)-1 size_stack(2)-1 0]);
    warning('off');
    part.Value=all_part.Value(:,:,i);
    warning('on');
    
    part_org=part.Value;
    
     if (strcmp(correction_flag,'pre Alignment')==0)
        part.Value=tom_apply_filter(part.Value,filter_param.filter.classify);
     end;
        
    part=double(part.Value).*mask_classify1.*mask_classify2;
    
     
    for ii=1:size_ref_stack(3) %loop over all references
        
        %ref=tom_emread(path_ref_stack,'subregion',[1 1 ii],[size_stack(1)-1 size_stack(2)-1 0]);
         warning('off');
        ref.Value=all_ref.Value(:,:,ii);
         warning('on');
        
        ref_org=ref.Value;
        if (strcmp(correction_flag,'pre Alignment')==0)
            ref.Value=tom_apply_filter(ref.Value,filter_param.filter.classify);
        end;
        
        ref=double(ref.Value).*mask_classify1.*mask_classify2;
        
        if (size_ref_stack(3)==1)
            ccc=1; ccc_pos=[1 1];
            break;
        end;
        
        
        if (strcmp(correction_flag,'pre Alignment'))
            %alg_iter=1;
            [tmp_a tmp_sh_st(ii,:) ccc(ii)]=tom_av2_align(ref,part,mask_align,mask_ccf_rot,mask_ccf_trans,filter_param.filter.align,alg_iter,disp_flag.show_alignment);
            rot_angle_st(ii)=tmp_a(1);
            ccc_pos(ii,:)=[1 1];
            ccf=ones(size(mask_align));
        else
            %calculate correlation in cart coordinates ... no Alignment
            ccf=tom_corr(part,ref,'norm');
            ccf=ccf.*mask_ccf_trans;
            [ccc_pos(ii,:) ccc(ii)]=tom_peak(ccf,'spline');
        end;
        
         % display 
        if (disp_flag.developer==1)
            title{1}=['Particle Nr.: ' num2str(i)]; 
            title{2}=['Projection Nr.: ' num2str(ii) ' angle: ' num2str(disp_info(2,ii)) ' ' num2str(disp_info(1,ii)) ]; 
            title{3}=['Pos.: ' num2str(ccc_pos(ii,1)) ' '  num2str(ccc_pos(ii,2)) ' CCC:' num2str(ccc(ii))];
            if (strcmp(correction_flag,'pre Alignment'))
                disp_classify(tom_shift(tom_rotate(part,rot_angle_st(ii)),[tmp_sh_st(ii,:)]),ref,ccf,title,ccc_pos(ii,:),'disp','up',disp_flag);   
            else
                disp_classify(part,ref,ccf,title,ccc_pos(ii,:),'disp','up',disp_flag);   
            end;
        end;
       
    end; % end of all references
    
    
    
    [max_value max_pos]=max(ccc); % find the best match
    
    if (max_value==0)
        disp('stop');
    end;
    
    % check out shift according to the correction flag
    if (strcmp(correction_flag,'no Alignment')==1)
        tmp_sh=ccc_pos(max_pos,:)-middle_part;
    end;
    
    if (strcmp(correction_flag,'pre Alignment')==1)
        tmp_sh=tmp_sh_st(max_pos,:);
        rot_angle=rot_angle_st(max_pos);
    end;
    
    if (strcmp(correction_flag,'post Alignment'))
        best_ref=tom_emread(path_ref_stack,'subregion',[1 1 max_pos],[size_stack(1)-1 size_stack(2)-1 0]);
        best_ref.Value=tom_apply_filter(best_ref.Value,filter_param.filter.classify);
        best_ref.Value=double(best_ref.Value).*mask_classify1.*mask_classify2;
        [rot_angle tmp_sh ccc_o]=tom_av2_align(best_ref.Value,part,mask_align,mask_ccf_rot,mask_ccf_trans,filter_param.filter.align,alg_iter,disp_flag.show_alignment);
        max_value=ccc_o;
    end;
    
    
    
    class_st(1,zz)=max_pos;
    class_st(2,zz)=max_value;
    class_st(3,zz)=tmp_sh(1);
    class_st(4,zz)=tmp_sh(2);
    class_st(5,zz)=rot_angle(1);
    zz=zz+1;
    
    % display 
    if (disp_flag.developer==1)
        title{1}=['Particle Nr: ' num2str(i)];
        title{2}=['Matched Projection Nr.: ' num2str(max_pos) ' angle: ' num2str(disp_info(2,max_pos)) ' ' num2str(disp_info(1,max_pos)) ];
        title{3}=['Pos.: '  num2str(max_pos) ' CCC:' num2str(max_value)];
        p=tom_emread(path_ref_stack,'subregion',[0 0 max_pos],[size_stack(1)-1 size_stack(2)-1 0]);
        if(strcmp(correction_flag,'no Alignment'))
            tmp_sh=[0 0];
            rot_angle(1)=0;
        end;

        disp_classify(tom_shift(tom_rotate(part,rot_angle(1)),[tmp_sh(1) tmp_sh(2)]) ,p.Value,ccc,title,[max_pos max_value],'disp','down',disp_flag);
    end;
    
    %command window print
    store.i=i;
    [store]=tom_disp_estimated_time(store,'estimate_time');
    [store]=tom_disp_estimated_time(store,'progress');
    
    
end; % end of all particles loop





function [rot_angle shift]=calc_alignment(ref,part,mask_cc,mask_cc_rot,mask_align,filter,middle_part)

part_org=part;
part=part.*mask_align;
ref=ref.*mask_align;

part_polar=tom_cart2polar(part);
ref=tom_bandpass(ref,filter(1),filter(2));
ref_polar=tom_cart2polar(ref);

ccf=tom_corr(part_polar,ref_polar,'norm');
ccf=ccf.*mask_cc_rot; 
[ccc_pos ccc]=tom_peak(ccf,'spline');

rot_angle=(180-(360./size(ref_polar,2).*(ccc_pos(2)-1)));

part_rot=tom_rotate(part_org,-rot_angle);
part_rot=part_rot.*mask_align;

ccf=tom_corr(part,part_rot,'norm');
ccf=ccf.*mask_cc;
[ccc_pos ccc]=tom_peak(ccf,'spline');

shift=ccc_pos-middle_part;

%shift=shift-25; %danger !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





%disp_classify(avg_st_im,1,1,1,1,'end',1,disp_flag);
function disp_classify(im1,im2,im3,title_in,peak,flag,pos_flag,disp_flag);

if (disp_flag.developer==0)
    return;
end

if (strcmp(flag,'new'))
     if (isempty(findobj('tag','classify')))
        figure; set(gcf,'tag','classify');
        if (disp_flag.show_alignment==1)
            set(gcf,'Position',[897 105 968 392]);
        else
            set(gcf,'Position',[28 163 1202 416]);    
        end;
     end;
end;
    
if (strcmp(flag,'disp'))
   
    figure((findobj('tag','classify')));
    if (strcmp(pos_flag,'up')==1)
        pos=[1 2 3];
    else
        pos=[4 5 6];
    end;
    subplot(2,3,pos(1)); tom_imagesc(im1); %subimage(double(tom_norm(im1',1)));
    title(title_in{1});
    subplot(2,3,pos(2)); tom_imagesc(im2);  %subimage(double(tom_norm(im2',1)));
    title(title_in{2});
    if (strcmp(pos_flag,'up'))
        subplot(2,3,pos(3)); subimage(double(tom_norm(im3',1)));title(title_in{3});
    else
        subplot(2,3,pos(3)); plot(im3); title(title_in{3});
    end;
    hold on;plot(peak(1),peak(2),'ro');hold off;drawnow;

end;

if (strcmp(flag,'end'))
  if (isempty(findobj('tag','classes_over')))
 %       figure; set(gcf,'tag','classes_over');
        %set(gcf,'Position',[6    33   843   434]);    
   end;
  
%   figure(findobj('tag','classes_over'));
   tom_dspcub(tom_norm(im1,1)); drawnow;
   set(gcf,'Position',[634   661   593   434]);
   drawnow; drawnow;
   if (isempty(findobj('tag','classify'))==0)
        %close(findobj('tag','classify'));
   end;
end;





