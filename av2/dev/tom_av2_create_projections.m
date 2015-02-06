function [align2d error_m]=tom_av2_create_projections(model,angles,align2d,iter_num,disp_flag)

%if (exist('proj','dir')==0)
    mkdir proj;
%end;

error_m=0;
mean_mod=mean2(model);

file_path=align2d(iter_num,1).file_path;

%mask model
mask=tom_create_mask(align2d(iter_num,1).filter.mask.model_mask_sp);
model=model.*mask;


%command window print
store.end=size(angles,2); store.num_of_mesure=1; store.mesured=0;
[store]=tom_disp_estimated_time(store,'start','projecting');

%disp_create_projections(1,'new',disp_flag); im_store=zeros(size(model,1),size(model,2),1);

for angle_idx=1:size(angles,2) % loop over all angles
    
    %project
    eu(angle_idx,:)=tom_sum_rotation([0 0 angles(2,angle_idx); 270 90 angles(1,angle_idx)],[0 0 0 ; 0 0 0]);
    rot_tmp=tom_rotate(single(model),eu(angle_idx,:),'linear');
    p(:,:,angle_idx)=double(sum(rot_tmp,3)); 
    
    %command window print
    store.i=angle_idx;
    [store]=tom_disp_estimated_time(store,'estimate_time');
    [store]=tom_disp_estimated_time(store,'progress');
    
    %display 
    if disp_flag == 1
        disp_reconstruction2d('addprojection',p(:,:,angle_idx));
        disp_reconstruction2d('update_surface',angle_idx);
        disp_reconstruction2d('moviepicture');
    end
    %im_store=disp_create_projections(p(:,:,angle_idx),'disp',disp_flag,im_store);
end

tom_emwrite([file_path '/proj/proj_' num2str(iter_num) '.em' ],p);
[a,b]=system(['chmod -R ugo+rwx ' file_path '/proj/']);

%fill up alignment structure ...bookeeping
for i=1:size(align2d,2)
    align2d(iter_num,1).proj.path=['proj/proj_'];
    align2d(iter_num,1).proj.ext=['.em'];
    align2d(iter_num,1).proj.number=size(angles,2);
    align2d(iter_num,1).angular_scan=angles;
    align2d(iter_num,1).angular_scan_euler=eu;
end;

%display
%disp_create_projections(im_store,'end',disp_flag);

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


