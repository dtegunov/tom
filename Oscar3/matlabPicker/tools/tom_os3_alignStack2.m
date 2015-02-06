function [stack_out all_tmpl cc_out mean_ccc tmpl_nr al_param]= tom_os3_alignStack2(stack,reference,mask,binning,norm_flag,filter,iterations,mask_trans,mask_rot,demo_mode,parlallel_mode)
%TOM_OS3_alignStack2 performs iterative rotational, translational alignment of an imageStack (im) relative to reference (ref)
%
%   [stack_out all_tmpl cc_out mean_ccc]= tom_os3_alignStack2(stack,reference,mask,binning,norm_flag,iterations,mask_trans,mask_rot)
%
%PARAMETERS
%
%  INPUT
%   im                  stack to be aligned
%   reference           reference image or image stack
%   mask                mask for both, reference and stack
%   binning             vector for binned calculation handles.storage.selectedclassTemp
%                       if second element is 0 transformation is aplyied to or
%                       particles;
%   norm_flag           norming for stack and reference 
%   filter              ('') filter parameter for bandpass filter (low high smooth ...check tom_bandpass docu)
%                             use '' for no filter  
%   iterations          vector for iterations first element number of
%                                             refinements
%                                             second element number of trans rot iterations (in tom_av2_align)          
%   mask_trans          mask for cartesien correlation (to limit translation)
%   mask_rot            mask for polar correlation (to limit rotation)
%   demo_mode           dmoe mode switch 1=on 2=off
%   parallel_mode       (1) flag which disables parfor loop due to shitty matlab
%                       release use (0) in nested application 
%                       parfor parfor end; end; is not possible any more  FUCK!!!!!!!!!!!!!!!!!!!11
%                      
%
%        
%
%  OUTPUT
%   stack_out           aligned stack
%   all_tmpl            template for each iteration
%   cc_out              ccc for all particles (last iteration)   
%   mean_cc             mean_cc over all iterations
%   tmpl_nr             template with max_cc
%   al_param            alignmetn param
%
%EXAMPLE
%
%
%[stack_out all_tmpl cc_out cc_sum]=tom_os3_alignStack2(stack,template,mask,[1 0],'mean0+1std','',[5 3]);
%
%[stack_out all_tmpl cc_out cc_sum]=tom_os3_alignStack2(stack,template,'default',[1 0],'mean0+1std',[4 55 5],[5 3],'default');
%
% center only no rotation calculated
% m_polar=tom_rectanglemask(zeros(64,256),[65 1]);
% [stack_out all_tmpl cc_out cc_sum]=tom_os3_alignStack2(stack,template,'default',[1 0],'mean0+1std',[4 55 5],[5 3],'default',m_polar);
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by /FB late/08
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

  
if nargin < 3  || isempty(mask)
    mask=ones(size(reference));
end;

if (strcmp(mask,'default'))
    mask=tom_spheremask(ones(size(reference,1),size(reference,1)),round(size(reference,1)./2)-1);
end;

if (nargin < 4) || isempty(binning)
    binning=[0 0];
end;

if (length(binning)==1)
    binning(2)=binning(1);
end;

if (nargin < 5)
    norm_flag='mean0+1std';
end;

if (nargin < 6)
    filter='';
end;

if (nargin < 7)
    iterations=[1 1];
end;

if  nargin < 8 || (isempty(mask_trans)) 
    mask_trans=ones(size(reference,1),size(reference,2));
end;
  
if (strcmp(mask_trans,'default'))
    szt=round(size(reference,1).*0.18); 
    mask_trans=tom_spheremask(ones(size(reference,2),size(reference,2)),szt);
end;

if (nargin < 9 || isempty(mask_rot)) || strcmp(mask_rot,'default')
    mask_rot=tom_cart2polar(ones(size(reference,1),size(reference,2)));
end;
  
if  nargin < 10
    demo_mode=0;
end;

if  nargin < 11
    parlallel_mode=1;
end;


filter_param.Apply = 0;

if (binning(2)~=binning(1))
    mask_org=mask;
    stack_org=stack;
    sz_stack_org=size(stack_org);
end;

if (binning(1) > 0)
    mask=double(tom_bin(mask,binning))>0;
    mask_trans=double(tom_bin(mask_trans,binning))>0;
    if (exist('mask_rot','var'))
        mask_rot=tom_bin(mask_rot,binning(1))>0;
        if sum((size(mask_rot)== size(double(tom_cart2polar(ones(size(mask_trans,1),size(mask_trans,2))))) ));
            warning off;
            mask_rot=imresize(mask_rot,size(tom_cart2polar(ones(size(mask_trans)))) )>0;
            warning on;
        end;
            
    else
        mask_rot=double(tom_cart2polar(ones(size(mask_trans,1),size(mask_trans,2))));
    end;
    
    
    
    
    for i=1:size(reference,3)
        new_r(:,:,i)=tom_bin(reference(:,:,i),binning);
    end;
    reference=new_r; clear('new_r');
end;

reference=tom_norm(reference,norm_flag);
sz_reference=size(reference);
if (length(sz_reference)==2)
    sz_reference(3)=1;
end;
sz_stack=size(stack);

if (length(sz_stack)==2)
    sz_stack=[sz_stack(1) sz_stack(2) 1];
end;

%pre process Stack
new_stack=zeros(sz_reference(1),sz_reference(2),sz_stack(3));
disp('pre processing started !');
parfor i=1:sz_stack(3)
    new_stack(:,:,i)=tom_norm((tom_bin(stack(:,:,i),binning(1))+2).*2,'mean0+1std');
    if (isempty(filter)==0)
        new_stack(:,:,i)=tom_norm(tom_bandpass(new_stack(:,:,i),filter(1),filter(2),filter(3)),'mean0+1std');
    end;
end;
stack=new_stack;
disp('pre Processing done !');

%allocate memory
stack_out=zeros(sz_reference(1),sz_reference(2),sz_stack(3));
mean_ccc=zeros(iterations(1),1);
cc_out=zeros(sz_stack(3),iterations(1));
tmpl_nr=zeros(sz_stack(3),1);
angle=zeros(sz_stack(3),3);
shift=zeros(3,sz_stack(3));

warning off;
%do alignment
tmpl=double(reference);
for ii=1:iterations(1)
    
    parfor i=1:sz_stack(3) %loop over all particles
    %    for i=1:sz_stack(3) %loop over all particles
        part=stack(:,:,i);
        [angle(i,:) shift(:,i) cc_out(i,ii) aligned_part tmpl_nr(i)]=aling_parts(part,tmpl,mask,mask_rot,mask_trans,filter_param,iterations(2),demo_mode);
        stack_out(:,:,i)=tom_norm(aligned_part,norm_flag,mask).*mask;
        %stack_out(:,:,i)=tom_norm(aligned_part,norm_flag);
    end;
    
    parfor i=1:sz_reference(3)
        %for i=1:sz_reference(3)
        idx=find(tmpl_nr==i);
        num_per_class=length(idx);
        if (num_per_class==0)
            num_per_class=1;
        end;
        tmpl(:,:,i)=sum(stack_out(:,:,idx),3);
        tmpl(:,:,i)=tmpl(:,:,i)./num_per_class;
    end;
    
    all_tmpl{ii}=tmpl;
    mean_ccc(ii)=mean(cc_out(:,ii));
    disp(['Iteration nr:' num2str(ii) ' done']);
end;
warning on;

al_param.shift=shift;
al_param.angle=angle;

%apply to unbinned particles
if (binning(1)~=binning(2))
    shift=shift.*2^binning(1)+(1.*2^binning(1));
    stack_out=zeros(sz_stack_org);
    warning off;
    parfor i=1:sz_stack_org(3)
        aligned_part=tom_norm(stack_org(:,:,i),norm_flag);
        aligned_part=tom_shift(tom_rotate(aligned_part,angle(i,1)),[shift(1,i) shift(2,i)]);
        stack_out(:,:,i)=tom_norm(aligned_part,norm_flag,mask_org).*mask_org;
    end;
    
    
    
    ref_rs=zeros(sz_stack_org(1),sz_stack_org(2),sz_reference(3));
    for i=1:iterations(1)-1
        tmpp=all_tmpl{i};
        for ii=1:sz_reference(3)
           ref_rs(:,:,ii)=imresize(tmpp(:,:,ii),[sz_stack_org(1) sz_stack_org(2)]);
        end
        all_tmpl_big{i}=ref_rs;
    end;
    clear('tmpl');
    tmpl=zeros(sz_stack_org(1),sz_stack_org(2),sz_reference(3));
    parfor i=1:sz_reference(3)
        idx=find(tmpl_nr==i);
        tmpl(:,:,i)=sum(stack_out(:,:,idx),3);
        tmpl(:,:,i)=tmpl(:,:,i)./length(idx);
    end;
    warning on;
    all_tmpl_big{iterations(1)}=tmpl;
    all_tmpl=all_tmpl_big;
end;

    
    
    
function [angle shift cc aligned_part pos]=aling_parts(part,tmpl,mask,mask_rot,mask_trans,filter_param,iterations,demo)


sz_t=size(tmpl);
if (length(sz_t)==2)
    sz_t(3)=1;
end;
    

angle_tmp=zeros(sz_t(3),3);
shift_tmp=zeros(3,sz_t(3));
aligned_part_tmp=zeros(sz_t);
cc_tmp=zeros(sz_t(3),1);

for i=1:sz_t(3)
    ref=tmpl(:,:,i);
    [angle_tmp(i,:) shift_tmp(:,i) cc_tmp(i) aligned_part_tmp(:,:,i)]=tom_av2_align(ref,part,mask,mask_rot,mask_trans,filter_param,iterations,demo);
end;


%determine maximum of all references
[cc pos]=max(cc_tmp);
angle=angle_tmp(pos,:);
shift=shift_tmp(:,pos);
aligned_part=aligned_part_tmp(:,:,pos);

        
        
        
        
        