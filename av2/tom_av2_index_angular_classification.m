function [align2d error_m]=tom_av2_index_angular_classification(align2d,count_st,disp_flag)
%TOM_AV2_INDEX_ANGULAR_CLASSIFICATION creates ...
%
%   [align2d error_m]=tom_av2_angular_classification(align2d,count_st,disp_flag)
%
%PARAMETERS
%
%  INPUT
%   align2d             ...
%   count_st            ...
%   disp_flag           ...
%  
%  OUTPUT
%   align2d             ...
%   error_m             ...
%
%EXAMPLE
%   ... = tom_av2_angular_classification(...);
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



%transfer often used data
error_m=0;
filtered_out=0;
hist_count=count_st.hist;
h=tom_reademheader(align2d(hist_count,1).rec.file.Stack_Path);

output_dir=align2d(count_st.hist,1).rec.file.Outputdir;
tmp=align2d(hist_count,1).rec.classify.alignment;
correction_flag=tmp{count_st.step};
angular_scan=align2d(hist_count,1).rec.project.angular_scan;
num_of_tasks=align2d(hist_count,1).rec.classify.parallel.number_of_tasks;
max_lost_packages=align2d(hist_count,1).rec.classify.parallel.packageloss;
jobmanager=align2d(hist_count,1).rec.classify.parallel.jobmanager;
size_stack=h.Header.Size;
file_path=[align2d(hist_count,1).rec.file.Outputdir '/' 'step' num2str(count_st.step)];
stack_path=align2d(hist_count,1).rec.file.Stack_Path;
ref_path=[output_dir '/step' num2str(count_st.step) '/proj/proj_' num2str(count_st.iteration) '.em'];
filter_param=align2d(hist_count,1).rec.classify.filter;                                                                                            
% initialize structure for class averages
avg_st_im=zeros(size_stack(1),size_stack(2),size(angular_scan,2));
avg_st_num=zeros(1,size(angular_scan,2));
stack_start=1;
stack_stop=align2d(hist_count,1).rec.classify.particles_vector(end);
stack_size=max(align2d(hist_count,1).rec.classify.particles_vector);

tree_split.Method= align2d(hist_count,1).rec.classify.indexs.treesplit;
tree_split.Values.Binning=align2d(hist_count,1).rec.classify.indexs.pca_binning; 
tree_split.Values.Eigenvalues=align2d(hist_count,1).rec.classify.indexs.pca_eigs; 


root_path_stack=fileparts(stack_path);
if (isempty(root_path_stack)==1)
    if (task_number > 1)
        error('stack path must be absolute!');
    end;
end;

root_path_ref=fileparts(ref_path);
if (isempty(root_path_ref)==1)
    if (task_number > 1)
        error('ref path must be absolute!');
    end;
end;


% %to be replaced by GUI-Values


%build up index structure
rest=strtok(ref_path,'.');
tree_split.num_of_classes=2;
tom_av2_index_calc2(ref_path,[rest '_index.em'],align2d(2,1).rec.classify.filter.mask.align,tree_split);



if (num_of_tasks==1)
    %just rip the local horst hoe hoe
   [class_st]=tom_av2_index_multi_ref_alignment(stack_path,ref_path,filter_param,correction_flag,2,[stack_start stack_stop],0,disp_flag,angular_scan);
   lost=0;
else
    
    [a,b]=system(['chmod -R ugo+rwx ' file_path '/tmp/']);
    %jm = findResource('jobmanager','name',jobmanager);
    jm = findResource('scheduler','type','jobmanager','lookupurl',jobmanager);
    % set(jm.jobs,'MaximumNumberOfWorkers',40);
    result_p=tom_get_paraell_results(jm,'hosts');
    disp(['classifying on ' num2str(size(result_p,1)) ' hosts']); 
    for iii=1:size(result_p,1)
        disp(result_p{iii});
    end;    
    clear('result_p');
    j = createJob(jm,'Name','demojob');

     pdd={'/fs/pool/pool-bmsan-apps/tom_dev/IOfun' '/fs/pool/pool-bmsan-apps/tom_dev/Filtrans/'  '/fs/pool/pool-bmsan-apps/tom_dev/Analysis/' '/fs/pool/pool-bmsan-apps/tom_dev/Misc/' '/fs/pool/pool-bmsan-apps/tom_dev/Sptrans/' ...
        '/fs/pool/pool-bmsan-apps/tom_dev/Geom/' '/fs/pool/pool-bmsan-apps/tom_dev/av2/'  };
    
    set(j,'FileDependencies',pdd)
    pdd2={[root_path_stack] [root_path_ref]};
    set(j,'PathDependencies',pdd2);
    rehash toolbox;

    %set(j,'FileDependencies',{'/fs/bmsan/apps/tom_dev/av2/tom_av2_multi_ref_alignment.m'})
    
    packages=tom_calc_packages(num_of_tasks,stack_stop);
    for i=1:num_of_tasks
        createTask(j,@tom_av2_index_multi_ref_alignment,1,{stack_path,ref_path,filter_param,correction_flag,2,packages(i,1:2),i,disp_flag});
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
        part=tom_emreadc([stack_path],'subregion',[1 1 i],[size_stack(1)-1 size_stack(2)-1 0]);
        part=part.Value;
        
        align2d(hist_count,i).rec.classify.angleclass.proj_nr=proj_nr;
        align2d(hist_count,i).rec.classify.angleclass.angle=align2d(hist_count,1).rec.project.angular_scan(:,proj_nr);
        align2d(hist_count,i).rec.classify.angleclass.angle_euler=align2d(hist_count,1).rec.project.angular_scan_euler(proj_nr,:);
        align2d(hist_count,i).rec.classify.angleclass.ccc_max=ccc;
        align2d(hist_count,i).rec.classify.angleclass.shiftxy=shift;
        align2d(hist_count,i).rec.classify.angleclass.rot=rot;
        
        if (ccc > 0) %filter to be done !
        
            if (strcmp(correction_flag,'no Alignment')==0)
                part=tom_shift(tom_rotate(part,rot),[shift]);
             end;
            avg_st_im(:,:,proj_nr)=avg_st_im(:,:,proj_nr)+part;
            avg_st_num(proj_nr)=avg_st_num(proj_nr)+1;
        else
            filtered_out=filtered_out+1;
        end;
    
    else
        %package got lost or classification failed
        align2d(hist_count,i).angleclass.proj_nr=-1;
        align2d(hist_count,i).angleclass.angle=[];
        align2d(hist_count,i).angleclass.angle_euler=[];
        align2d(hist_count,i).angleclass.ccc_max=0;
        align2d(hist_count,i).angleclass.shiftxy=[];
    end;

    %command window print
    store.i=i;
    [store]=tom_disp_estimated_time(store,'estimate_time');
    
end;
 
tom_emwrite([file_path '/avg/avg_' num2str(count_st.iteration) '.em'],avg_st_im);
align2d(hist_count,1).rec.classify.avg.path=[file_path '/avg/avg_' num2str(count_st.iteration) '.em'];
align2d(hist_count,1).rec.classify.avg.ext='.em';
align2d(hist_count,1).rec.classify.avg.number=size(angular_scan,2);
align2d(hist_count,1).rec.classify.avg.num_array=avg_st_num;

 if (sum(align2d(hist_count,1).rec.classify.avg.num_array)~=(stack_size)-lost-filtered_out)
        error('classification corrupt');
 end;

 save([file_path '/align/align2d'],'align2d');
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

