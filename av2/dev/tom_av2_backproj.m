function [vol error_m]=tom_av2_backproj(align2d,iter_num,disp_flag)

%transfer often used data
error_m=0;
angular_scan_eu=align2d(iter_num,1).angular_scan_euler;
angular_scan=align2d(iter_num,1).angular_scan;
size_stack=align2d(iter_num,1).stack_size;
size_vol=align2d(iter_num,1).model.size;
size_part=align2d(iter_num,1).model.size_part;
weight_flag=align2d(iter_num,1).reconstruction.weight_flag;

vol=zeros(size_vol,'single'); %allocate some memory 

%command window print
store.end=size(angular_scan,2); store.num_of_mesure=1;  store.mesured=0;
[store]=tom_disp_estimated_time(store,'start','backproj');

%display
disp_back_proj(1,'new',disp_flag);

%find used angleclasses ...build theta vector for weighting 
zz=1;
for i=1:size(angular_scan,2);
     if (align2d(iter_num,1).avg.num_array(i) > 0 ) 
        eu_4_weight(zz,:)=angular_scan_eu(i,:);
        zz=zz+1;
     end;
end;


mask=tom_create_mask(align2d(iter_num,1).filter.mask.classify1);
part=zeros(size_stack(1),size_stack(2));

% perform weighting and backprojection
for i=1:size(angular_scan,2)
    
        eu=align2d(iter_num,1).angular_scan_euler(i,:);
        num_of_p=align2d(iter_num,1).avg.num_array(i);

        if (num_of_p~=0) 

            %read class avg
            part=tom_emreadc(align2d(iter_num,1).avg.path,'subregion',[1 1 i],[size_stack(1)-1 size_stack(2)-1 0]);
            part=double(part.Value);
            part=(part./num_of_p); % weight it according to number
            if (strcmp(weight_flag,'projection')==1)
                [w]=tom_calc_weight_functionc([size_stack(1) size_stack(2)],eu_4_weight,size_part(1),eu);
                part=tom_apply_weight_function(part,w);
            end;
            part=single(part.*mask);
            vol=single(vol);
            tom_backproj3d_euler(vol,part,eu(1),eu(2),eu(3),[0 0 0]);

        end;

        %display
        disp_back_proj(part,'progress',disp_flag);
        %command window print
        store.i=i;
        store=tom_disp_estimated_time(store,'estimate_time');
        store=tom_disp_estimated_time(store,'progress');
    
end;

if (strcmp(weight_flag,'volume')==1)
    [w]=tom_calc_weight_functionc(size_vol,eu_4_weight,size_part(1));
    vol=tom_apply_weight_function(vol,w);
end;





function disp_back_proj(im,flag,disp_flag)

if (disp_flag==0)
    return;
end

if strcmp(flag,'new')
    if (isempty(findobj('tag','backproj')))
        figure; set(gcf,'tag','backproj'); 
        drawnow; set(gcf,'pos',[1235 143 565 434]); drawnow; drawnow;
    end;
 end;

if strcmp(flag,'progress')
    if (isempty(findobj('tag','backproj')))
        figure; set(gcf,'tag','backproj'); 
        drawnow; set(gcf,'pos',[1235 143 565 434]); drawnow; drawnow;
    end;
    figure((findobj('tag','backproj')) ); tom_imagesc(im,'noinfo');
end;

if strcmp(flag,'end')
    if (isempty(findobj('tag','backproj'))==0)
        %close(findobj('tag','backproj'));
    end;
end;

