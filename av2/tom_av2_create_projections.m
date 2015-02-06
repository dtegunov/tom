function [align2d error_m]=tom_av2_create_projections(model,align2d,count_st,disp_st)
%TOM_AV2_CREATE_PROJECTIONS creates ...
%
%   [align2d error_m]=tom_av2_create_projections(model,align2d,count_st,disp_st)
%
%PARAMETERS
%
%  INPUT
%   model               ...
%   align2d             ...
%   count_st            ...
%   disp_st             ...
%  
%  OUTPUT
%   align2d             ...
%   error_m             ...
%
%EXAMPLE
%   ... = tom_av2_create_projections(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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




hist_count=count_st.hist;

error_m=0;


%build file structure
output_dir=align2d(hist_count,1).rec.file.Outputdir;

angles=tom_create_angles_st(align2d(hist_count,1).rec.project.angle(count_st.step,1:3),align2d(hist_count,1).rec.project.angle(count_st.step,4:6),align2d(hist_count,1).rec.project.scheme);

if (isempty(angles))
    error_m=1;
    return;
end;



%mask model
mask=tom_create_mask(align2d(hist_count,1).rec.project.filter.mask);
model=model.*mask;


%command window print
store.end=size(angles,2); store.num_of_mesure=1; store.mesured=0;
[store]=tom_disp_estimated_time(store,'start','projecting');


disp_create_projections(1,'new',disp_st.developer); 
im_store=zeros(size(model,1),size(model,2),1);

for angle_idx=1:size(angles,2) % loop over all angles
    
    %project
   % eu(angle_idx,:)=tom_sum_rotation([0 0 angles(1,angle_idx); 270 90 angles(2,angle_idx)],[0 0 0 ; 0 0 0]);
   [a,tmp_eu] = tom_eulerconvert_xmipp(angles(1,angle_idx),angles(2,angle_idx),0); 
   eu(angle_idx,:)=tmp_eu;
   %eu(angle_idx,:)=[angles(1,angle_idx) angles(2,angle_idx) 0];
   %eu(angle_idx,:)=[angles(2,angle_idx) angles(1,angle_idx)];
    rot_tmp=tom_rotate(single(model),eu(angle_idx,:),'linear');
    %rot_tmp=tom_xmipp_rotate(model,eu(angle_idx,:)');
    p(:,:,angle_idx)=double(squeeze(sum(rot_tmp,3))); 
    
    %command window print
    store.i=angle_idx;
    [store]=tom_disp_estimated_time(store,'estimate_time');
    [store]=tom_disp_estimated_time(store,'progress');
    
    %display 
    im_store=disp_create_projections(p(:,:,angle_idx),'disp',disp_st.developer,im_store);
end;

tom_emwrite([output_dir '/step' num2str(count_st.step) '/proj/proj_' num2str(count_st.iteration) '.em' ],p);
[a,b]=system(['chmod -R ugo+rwx ' output_dir '/step' num2str(count_st.step) '/proj/']);

%fill up alignment structure ...bookeeping
align2d(hist_count,1).rec.project.path=[output_dir '/step' num2str(count_st.step) '/proj/'];
align2d(hist_count,1).rec.project.ext=['.em'];
align2d(hist_count,1).rec.project.number=size(angles,2);
align2d(hist_count,1).rec.project.angular_scan=angles;
align2d(hist_count,1).rec.project.angular_scan_euler=eu;

%display
disp_create_projections(im_store,'end',disp_st.developer);





function im_stack=disp_create_projections(im,flag,disp_flag,im_stack);

if (disp_flag==0)
    return;
end;

if (strcmp(flag,'new'))
    if (isempty(findobj('tag','projections')))
        figure; set(gcf,'tag','projections');
        drawnow; set(findobj('tag','projections'),'pos',[1235 661 565 434]);
        drawnow; drawnow;
    end;
end;

if (strcmp(flag,'disp'))
    if (isempty(findobj('tag','projections')))
        figure; set(gcf,'tag','projections');
    end;
    figure((findobj('tag','projections'))); tom_imagesc(im,'noinfo'); drawnow;
    im_stack(:,:,size(im_stack,3))=im;
    im_stack(:,:,size(im_stack,3)+1)=im; % write it 2 times to save a counter var o' la la 
end;

if (strcmp(flag,'end'))
    if (isempty(findobj('tag','projections'))==0)
        %close(findobj('tag','projections'));
    end;
    if (isempty(findobj('tag','projections_overview')))
        figure; set(gcf,'tag','projections_overview');
    end;
    figure(findobj('tag','projections_overview'));
    tom_dspcub(im(:,:,(1:(size(im,3)-1))) ); drawnow;
    set(findobj('tag','projections_overview'),'pos',[28   659   597   436]);
    drawnow; drawnow;
end;


