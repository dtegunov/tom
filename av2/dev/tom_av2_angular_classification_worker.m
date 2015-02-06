function class_st=tom_av2_angular_classification_worker(path_align_st,package,iter_num,disp_flag,task_number)

status=0;

if (nargin==4)
    task_number=0;
end;

%add path for paraellisation
addpath('/fs/bmsan/apps/tom/IOfun');
addpath('/fs/bmsan/apps/tom/Filtrans');
addpath('/fs/bmsan/apps/tom/Geom/');
addpath('/fs/bmsan/apps/tom/Analysis/');
addpath('/fs/bmsan/apps/tom_dev/Misc/');
addpath('/fs/bmsan/apps/tom_dev/Sptrans/');
rehash toolbox;

 
% transfer often used data
load(path_align_st);
angular_scan=align2d(iter_num,1).angular_scan;
proj=align2d(iter_num,1).proj;
size_stack=align2d(iter_num,1).stack_size;
middle_part(1)=(size_stack(1)./2)+1;
middle_part(2)=(size_stack(2)./2)+1;
shift_corr_flag=align2d(iter_num,1).shiftcorr;
file_path=align2d(iter_num,1).file_path;


mask=tom_spheremask(ones(size_stack(1)),size_stack(1)./2-19,2,[size_stack(1)./2+1 size_stack(1)./2+1 1]);

%command window print
store.end=size_stack(3); store.num_of_mesure=1;  store.mesured=0;
[store]=tom_disp_estimated_time(store,'start','classify');

%display   
disp_classify(1,1,1,1,1,'new',1,disp_flag); 

mask_cc=tom_spheremask(ones(size_stack(1)),round(size_stack(1)./3.5),5);

zz=1;
for i=package(1):package(2) % loop over all particles

    %read the particle of stack
    part=tom_emreadc([file_path '/' align2d(iter_num,1).filename],'subregion',[0 0 i],[size_stack(1)-1 size_stack(2)-1 0]);
    part=double(part.Value);
    part=tom_bandpass(part,1,32).*mask;
    
    for angle_idx=1:size(angular_scan,2) %loop over all projections
        
        ref=tom_emreadc([file_path '/' proj.path num2str(angle_idx) proj.ext]); 
        ref=double(ref.Value);
        ref=tom_bandpass(ref,1,25).*mask;
        
        %calculate correlation
        %ccf=tom_corr(tom_smooth(part,round(size_stack(1)./10)),tom_smooth(ref,round(size_stack(1)./10)),'norm');
        ccf=tom_corr(part,ref,'norm');
        ccf=ccf.*mask_cc;
        [ccc_pos(angle_idx,:) ccc(angle_idx)]=tom_peak(ccf);
        
        % display 
        title{1}=['Particle Nr.: ' num2str(i)]; 
        title{2}=['Projection Nr.: ' num2str(angle_idx) ' angle: ' num2str(angular_scan(2,angle_idx)) ' ' num2str(angular_scan(1,angle_idx)) ]; 
        title{3}=['Pos.: ' num2str(ccc_pos(angle_idx,1)) ' '  num2str(ccc_pos(angle_idx,2)) ' CCC:' num2str(ccc(angle_idx))];
        disp_classify(part,ref,ccf,title,ccc_pos(angle_idx,:),'disp','up',disp_flag);   
        
    end; % end of all projections loop
    
    [max_value max_pos]=max(ccc); % find the best match
    tmp_sh=ccc_pos(max_pos,:)-middle_part;
    
    class_st(1,zz)=max_pos;
    class_st(2,zz)=max_value;
    class_st(3,zz)=tmp_sh(1);
    class_st(4,zz)=tmp_sh(2);
    class_st(5,zz)=0;
    zz=zz+1;
    
    % display 
    title{1}=['Particle Nr: ' num2str(i)]; 
    title{2}=['Matched Projection Nr.: ' num2str(max_pos) ' angle: ' num2str(angular_scan(2,max_pos)) ' ' num2str(angular_scan(1,max_pos)) ]; 
    title{3}=['Pos.: '  num2str(max_pos) ' CCC:' num2str(max_value)];
    p=tom_emreadc([file_path '/' proj.path num2str(max_pos) proj.ext]);
    disp_classify(part,p.Value,ccc,title,[max_pos max_value],'disp','down',disp_flag);   
    
    %command window print
    store.i=i;
    [store]=tom_disp_estimated_time(store,'estimate_time');
    [store]=tom_disp_estimated_time(store,'progress');
    
    
end; % end of all particles loop



[a,b]=system(['chmod -R ugo+rwx ' file_path '/tmp/']);


%disp_classify(avg_st_im,1,1,1,1,'end',1,disp_flag);



function disp_classify(im1,im2,im3,title_in,peak,flag,pos_flag,disp_flag);

if (disp_flag==0)
    return;
end

if (strcmp(flag,'new'))
     if (isempty(findobj('tag','classify')))
        figure; set(gcf,'tag','classify');
        set(gcf,'Position',[28 163 1202 416]);    
     end;
end;
    
if (strcmp(flag,'disp'))
   
    figure((findobj('tag','classify')));
    if (strcmp(pos_flag,'up')==1)
        pos=[1 2 3];
    else
        pos=[4 5 6];
    end;
    subplot(2,3,pos(1)); subimage(double(tom_norm(im1',1)));title(title_in{1});
    subplot(2,3,pos(2)); subimage(double(tom_norm(im2',1)));title(title_in{2});
    if (strcmp(pos_flag,'up'))
        subplot(2,3,pos(3)); subimage(double(tom_norm(im3',1)));title(title_in{3});
    else
        subplot(2,3,pos(3)); plot(im3); title(title_in{3});
    end;
    hold on;plot(peak(1),peak(2),'ro');hold off;drawnow;

end;

if (strcmp(flag,'end'))
  if (isempty(findobj('tag','classes_over')))
        figure; set(gcf,'tag','classes_over');
        %set(gcf,'Position',[6    33   843   434]);    
   end;
  
   figure(findobj('tag','classes_over'));
   tom_dspcub(tom_norm(im1,1)); drawnow;
   set(gcf,'Position',[634   661   593   434]);
   drawnow; drawnow;
   if (isempty(findobj('tag','classify'))==0)
        %close(findobj('tag','classify'));
   end;
end;





