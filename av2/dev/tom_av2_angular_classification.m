function [align2d error_m]=tom_av2_angular_classification(align2d,iter_num,disp_flag)

if (exist('avg','dir')==0)
    mkdir avg;
end;

if (exist('align','dir')==0)
    mkdir align;
end;


error_m=0;
filtered_out=0;
correction_flag=align2d(iter_num,1).corr_flag;
angular_scan=align2d(iter_num,1).angular_scan;
num_of_tasks=align2d(iter_num,1).paraell.number_of_tasks;
max_lost_packages=align2d(iter_num,1).paraell.packageloss;
jobmanager=align2d(iter_num,1).paraell.jobmanager;
size_stack=align2d(iter_num,1).stack_size;
file_path=align2d(iter_num,1).file_path;
stack_path=[align2d(iter_num,1).filename];
ref_path=[file_path '/proj/proj_' num2str(iter_num) '.em'];
filter_param=align2d(iter_num,1).filter;                                                                                            
% initialize structure for class averages
avg_st_im=zeros(size_stack(1),size_stack(2),size(angular_scan,2));
avg_st_num=zeros(1,size(angular_scan,2));

save([file_path '/align/align2d'],'align2d');

if (num_of_tasks==1)
    %just rip the local horst hoe hoe
   [class_st]=tom_av2_multi_ref_alignment(stack_path,ref_path,filter_param,correction_flag,2,[1 size_stack(3)],0,disp_flag,angular_scan);
   lost=0;
else
    
    [a,b]=system(['chmod -R ugo+rwx ' file_path '/tmp/']);
    jm = findResource('jobmanager','name',jobmanager);
    set(jm.jobs,'MaximumNumberOfWorkers',40);
    result_p=tom_get_paraell_results(jm,'hosts');
    disp(['classifying on ' num2str(size(result_p,1)) ' hosts']); 
    for iii=1:size(result_p,1)
        disp(result_p{iii});
    end;    
    clear('result_p');
    j = createJob(jm,'Name','demojob');
    set(j,'FileDependencies',{'/fs/bmsan/apps/tom_dev/av2/tom_av2_multi_ref_alignment.m'})
    packages=tom_calc_packages(num_of_tasks,size_stack(3));
    for i=1:num_of_tasks
        createTask(j,@tom_av2_multi_ref_alignment,1,{stack_path,ref_path,filter_param,correction_flag,2,packages(i,1:2),i});
    end;
    submit(j);
    %tom_disp_paraell_progress(j,packages(:,3));
    waitForState(j);
    out = getAllOutputArguments(j);
    [result_p errorsum]=tom_get_paraell_results(j);
    if (errorsum > 0);
        for i=1:size(result_p,1)
            disp(['Hostname: ' result_p{i,1} '  Error: ' result_p{i,2}]);
        end;
    end;
    destroy(j);
    [class_st lost error_m]=built_class_st(out,packages,max_lost_packages);
    
    if (error_m==1)
        return;
    end;
    
    
end;

%command window print
store.end=size(class_st,2); store.num_of_mesure=1;  store.mesured=0;
[store]=tom_disp_estimated_time(store,'start','averaging');


for i=1:size(class_st,2)
    proj_nr=class_st(1,i);
    ccc=class_st(2,i);
    shift=[class_st(3,i) class_st(4,i)];
    rot=class_st(5,i);
    
    if  (proj_nr~=-1)

        %read particle
        part=tom_emreadc([align2d(1,1).filename],'subregion',[1 1 i],[size_stack(1)-1 size_stack(2)-1 0]);
        part=part.Value;
        
        %norm particle
        %part=(part-mean2(part))./mean2(part);
        
        
        
        align2d(iter_num,i).angleclass.proj_nr=proj_nr;
        align2d(iter_num,i).angleclass.angle=align2d(iter_num,1).angular_scan(:,proj_nr);
        align2d(iter_num,i).angleclass.angle_euler=align2d(iter_num,1).angular_scan_euler(proj_nr,:);
        align2d(iter_num,i).angleclass.ccc_max=ccc;
        align2d(iter_num,i).angleclass.shiftxy=shift;
        align2d(iter_num,i).angleclass.rot=rot;
        
        if (ccc > 0) %filter to be done !
        
            if (strcmp(correction_flag,'no_correction')==0)
                part=tom_shift(tom_rotate(part,rot),[shift]);
            end;
            avg_st_im(:,:,proj_nr)=avg_st_im(:,:,proj_nr)+part;
            avg_st_num(proj_nr)=avg_st_num(proj_nr)+1;
        else
            filtered_out=filtered_out+1;
        end;
    
    else
        %package got lost or classification failed
        align2d(iter_num,i).angleclass.proj_nr=-1;
        align2d(iter_num,i).angleclass.angle=[];
        align2d(iter_num,i).angleclass.angle_euler=[];
        align2d(iter_num,i).angleclass.ccc_max=0;
        align2d(iter_num,i).angleclass.shiftxy=[];
    end;

    %command window print
    store.i=i;
    [store]=tom_disp_estimated_time(store,'estimate_time');
    [store]=tom_disp_estimated_time(store,'progress');
end;
 
tom_emwrite([file_path '/avg/avg_' num2str(iter_num) '.em'],avg_st_im);
align2d(iter_num,1).avg.path=[file_path '/avg/avg_' num2str(iter_num) '.em'];
align2d(iter_num,1).avg.ext=['.em'];
align2d(iter_num,1).avg.number=size(angular_scan,2);
align2d(iter_num,1).avg.num_array=avg_st_num;

 if (sum(align2d(iter_num,1).avg.num_array)~=size_stack(3)-lost-filtered_out)
        error('classification corrupt');
 end;

[a,b]=system(['chmod -R ugo+rwx ' file_path '/avg/']);


function [class_st lost error_m]=built_class_st(out,packages,max_lost_packages)

error_m=0;
zz=1;
lost=0;

if (isempty(out)==1)
    error_m=1;
    disp(['ERROR: All Packages lost !']);
    disp(['...restarting current iteration !']);
    error_m=1;
    class_st='';
    return;
end;

for i=1:size(out)
    if (isempty(out{i})==1)
        for ii=1:packages(i,3);
            class_st(:,zz)=[-1 -1 -1 -1 -1];
            zz=zz+1;
            lost=lost+1;
            disp('warning: Package lost');
        end;
        continue;
    end;
    tmp_st=out{i};
    for ii=1:packages(i,3);
        class_st(:,zz)=tmp_st(:,ii);
        zz=zz+1;
    end;
end;

if (lost > max_lost_packages)
    error_m=1;
    disp(['ERROR: Maximum number of lost packages reached !']);
    disp(['...restarting current iteration !']);
   
end;




