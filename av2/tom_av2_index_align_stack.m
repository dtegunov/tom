function [align2d new_ref]=tom_av2_index_align_stack(stack_path,ref_path,align2d,stack_alg_path,ref_alg_path,filter_param,paraell_param,iterations,demo,output_flag,ang_samp,inp_alg,tree_split)
%TOM_AV2_ALIGN_STACK performs iterative rotational, translational alignment of an image stack (im) relative to reference stack (ref)
%
%   [align2d new_ref]=tom_av2_align_stack(stack_path,ref_path,align2d,stack_alg_path,ref_alg_path,filter_param,paraell_param,iterations,demo,output_flag)
%
%PARAMETERS
%
%  INPUT
%   stack_path          path reference stack 
%   ref_path            path image stack to be aligned  
%   align2d             path alignment structure 
%   stack_alg_path      path for the aligned stack  ... use '' for
%   no output
%   ref_alg_path        path for the new reference stack ... use ''for no output
%   filter_param        contains filter and mask structures: 
%                                                filter_param.mask.classify1 
%                                                filter_param.mask.classify2
%                                                filter_param.mask.align 
%                                                filter_param.mask.ccf_rot
%                                                filter_param.mask.ccf_trans
%                                                filter_param.filter.classify
%                                                filter_param.filter.align
%                                                
%                                                
%                                                filter_st -structure
%                                                filter_st.Apply:  1 apply filter, 2 use default values, 0 not
%                                                filter_st.Value:  vector with parameters [low, high, smooth ...]
%                                                filter_st.Method: i.e. 'circ', 'quadr', 'bandpass'
%                                                filter_st.Space: 'real' or 'fourier'
%                                                filter_st.Times: 'apply filter n-times                                  
%
%
%                                                mask_st -structure
%                                                mask_st.Apply:  1 built mask according to the struct values, 2 use default values, 0 create ones
%                                                mask_st.Value:  vector with parameters [size,radius,smooth,center ...]
%                                                mask_st.Method: 'sphere' 'sphere3d' 'cylinder3d' 'rectangle'  
%                                                
%
%
%   paraell_param      structure for paraell processing ... use '' for one CU
%                                            
%                                                paraell_param.jobmanager:      Name of the jobmanager
%                                                paraell_param.packageloss:     maximum allowed package loss [0..1]
%                                                paraell_param.number_of_tasks: number of tasks in which the job will be split
%                                                paraell_param.workers.min:     minimum number of workers to use
%                                                paraell_param.workers.max:     maximum number of workers to use
%                                                paraell_param.timeout:         timeout value for each task in seconds
%                                                paraell_param.restart_workers: 1 = restart workers bef
%
%   iterations          [refinement alignment]   
%
%
%   demo                flag (1: show multiref, 2: show alignment  3: show multiref+alignment  0: off) for demo mode via graphical interface
%                        works only with on CU 
%                    
%   output_flag:(optional)  flag for substacks set 'sub_stacks' for
%   output by no no substacks are written !  
%  
%  OUTPUT
%   align2d             alignment structure with updated values for translation and rotation         
%   new_ref             aligned referece stack
%
%EXAMPLE
%   tom_amira_createisosurface(...);
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


%parse inputs!
h=tom_reademheader(stack_path);
max_num_iter=100000;

switch nargin
    case 2
        align2d='';
        stack_alg_path='';
        ref_alg_path='';
        filter_param=tom_av2_build_filter_param([h.Header.Size(1) h.Header.Size(2)],'default','tom_av2_multi_ref_alignment');
        paraell_param=tom_build_paraell_param('');
        iterations=[1 1];
        demo=0;
        output_flag='';
    case 3
        stack_alg_path='';
        ref_alg_path='';
        filter_param=tom_av2_build_filter_param([h.Header.Size(1) h.Header.Size(2)],'default','tom_av2_multi_ref_alignment');
        paraell_param=tom_build_paraell_param('');
        iterations=[1 1];
        demo=0;
        output_flag='';
    case 4
        ref_alg_path='';
        filter_param=tom_av2_build_filter_param([h.Header.Size(1) h.Header.Size(2)],'default','tom_av2_multi_ref_alignment');
        paraell_param=tom_build_paraell_param('');
        iterations=[1 1];
        demo=0;
        output_flag='';
    case 5
        filter_param=tom_av2_build_filter_param([h.Header.Size(1) h.Header.Size(2)],'default','tom_av2_multi_ref_alignment');
        paraell_param=tom_build_paraell_param('');
        iterations=[1 1];
        demo=0;
        output_flag='';
    case 6
        paraell_param=tom_build_paraell_param('');
        iterations=[1 1];
        demo=0;
        output_flag='';
    case 7
        iterations=[1 1];
        demo=0;
        output_flag='';
    case 8
        demo=0;
        output_flag'';
    case 9
         output_flag='';
    case 10
        ang_samp=2;
    case 11
        ang_samp=2;
    case 12
        tree_split.Metod='classify';
        tree_split.Values.Binning=1;
        tree_split.Values.Eigenvalues(1)=1;
        tree_split.Values.Eigenvalues(2)=10;
    case 13        
        
        
    otherwise
        error('wrong number of parameters!');
end;
  
if (isempty(filter_param)==1)
    filter_param=tom_av2_build_filter_param('',[h.Header.Size(1) h.Header.Size(2)],'default','tom_av2_multi_ref_alignment');
end;

filter_param=tom_av2_build_filter_param(filter_param,[h.Header.Size(1) h.Header.Size(2)],'default','tom_av2_multi_ref_alignment');

if (isempty(paraell_param)==1)
    paraell_param=tom_build_paraell_param('');
end;
paraell_param=tom_build_paraell_param(paraell_param);

if (isempty(iterations)==1)
    iterations=[1 1];
end;

if (isempty(align2d)==0)
    if (isstruct(align2d)==0)
        load(align2d);
    end;
end;




[path_tmp name_tmp ext_ref_alg_path]=fileparts(ref_alg_path);
ref_alg_path=[path_tmp '/' name_tmp];
if (isempty(ext_ref_alg_path)==1)
    ext_ref_alg_path='.em';
end;



error_m=0;


correction_flag=inp_alg;


shift_corr_flag=1;
h=tom_reademheader(stack_path);
size_stack=h.Header.Size;
h=tom_reademheader(ref_path);


size_ref_stack=h.Header.Size';
iter_num=1;
disp_flag=0;
angular_scan=0;
num_of_tasks=paraell_param.number_of_tasks;
max_lost_packages=round(paraell_param.packageloss.*num_of_tasks);
jobmanager=paraell_param.jobmanager;
number_of_alg_iter=iterations(1);
number_of_alg_steps=iterations(2);
mask=tom_create_mask(filter_param.mask.classify1);

%check for absolute path
tmp_path=fileparts(ref_alg_path);
if (isempty(tmp_path)==1)
    if (num_of_tasks > 1)
        ref_alg_path=[pwd '/' ref_alg_path];
    end;
end;

tmp_path=fileparts(stack_path);
if (isempty(tmp_path)==1)
    if (num_of_tasks > 1)
        stack_path=[pwd '/' stack_path];
    end;
end;


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


tmp_path=fileparts(ref_path);
if (isempty(tmp_path)==1)
    if (num_of_tasks > 1)
        ref_path=[pwd '/' ref_path];
    end;
end;



if (size_ref_stack(3)==1)
    output_flag'';
end;

step=1;



%pre align img Stack
disp('pre alg img-stack');
%  [ttt_align2d new_ref]=tom_av2_align_stack(stack_path,'/fs/sun15/lv01/pool/pool-nickell/tmp/index_test/mooped2/no_noise/model/tmp_ref.em',align2d,'/fs/sun15/lv01/pool/pool-nickell/tmp/index_test/mooped2/no_noise/model/stack_pre_alg.em','/fs/sun15/lv01/pool/pool-nickell/tmp/index_test/mooped2/no_noise/model/ref_outt.em',filter_param,paraell_param,[5 5],demo,output_flag);
%  stack_path='/fs/sun15/lv01/pool/pool-nickell/tmp/index_test/mooped2/no_noise/model/stack_pre_alg.em';

for iiii=1:max_num_iter

    
    % initialize structure for class averages
    avg_st_im=zeros(size_ref_stack);
    avg_st_num=zeros(1,size_ref_stack(3));

    
    if strcmp(correction_flag,'no Alignment')
        %sample it in-plane
        %mm2=tom_rectanglemask(zeros(40),[38 20],3);

        im=tom_emread(ref_path); 
        im=tom_apply_filter(im.Value,filter_param.filter.classify);
        im=tom_av2_index_inp_sample_refstack(im,ang_samp);
        rest=strtok(ref_path,'.'); tom_emwrite([rest '_sample.em'],im);
        t_ref_path=[rest '_sample.em'];
    else
        %pre align the stack 2 reduce inplane variance
         pre_align_flag=1;
        if (pre_align_flag==1)
            disp(['Pre align Stack ....']);
           
            ref_tmp=tom_emread(ref_path);
            tom_emwrite('tmp_ref.em',ref_tmp.Value(:,:,1));
            old_num=paraell_param.number_of_tasks;
            paraell_param.number_of_tasks=1;
         %   [tt_align2d new_ref]=tom_av2_align_stack(ref_path,'tmp_ref.em',align2d,'ref_pre_alg.em','/fs/sun15/lv01/pool/pool-nickell/tmp/index_test/mooped2/no_noise/model/ref_out.em',filter_param,paraell_param,[5 5],demo,output_flag);
            paraell_param.number_of_tasks=old_num;
            t_ref_path='ref_pre_alg.em';
        else
            t_ref_path=ref_path;
        end;

       
    end;
    
    %build up index structure
    rest=strtok(ref_path,'.');
    tree_split.num_of_classes=4;
    %tom_av2_index_calc2(t_ref_path,[rest '_index.em'],filter_param.mask.classify1,tree_split);
    disp('DANGER Index calculatio is commented!!');
    
    
    if (num_of_tasks==1)
        %just rip the local horst hoe hoe
        [class_st]=tom_av2_index_multi_ref_alignment(stack_path,ref_path,filter_param,correction_flag,number_of_alg_iter,[1 size_stack(3)],0,demo);
        lost=0;
    else

        [a,b]=system(['chmod ugo+rwx ' stack_path]);
        [a,b]=system(['chmod ugo+rwx ' ref_path]);
        
        jm = findResource('scheduler','type','jobmanager','lookupurl',jobmanager);
        
        
        result_p=tom_get_paraell_results(jm,'hosts');
        disp(['classifying on ' num2str(size(result_p,1)) ' hosts']);
        for iii=1:size(result_p,1)
            disp(result_p{iii});
        end;
        clear('result_p');
        j = createJob(jm,'Name','demojob');
        
        %add path of functions and data
        pdd={'/fs/pool/pool-bmsan-apps/tom_dev/IOfun' '/fs/pool/pool-bmsan-apps/tom_dev/Filtrans/'  '/fs/pool/pool-bmsan-apps/tom_dev/Analysis/' '/fs/pool/pool-bmsan-apps/tom_dev/Misc/' '/fs/pool/pool-bmsan-apps/tom_dev/Sptrans/' ...
            '/fs/pool/pool-bmsan-apps/tom_dev/Geom/' '/fs/pool/pool-bmsan-apps/tom_dev/av2/'  };
        set(j,'FileDependencies',pdd)
        pdd2={[root_path_stack] [root_path_ref]};
        set(j,'PathDependencies',pdd2);
        rehash toolbox;
        
        packages=tom_calc_packages(num_of_tasks,size_stack(3));
        for i=1:num_of_tasks
            createTask(j,@tom_av2_index_multi_ref_alignment,1,{stack_path,ref_path,filter_param,correction_flag,number_of_alg_iter,packages(i,1:2),i,demo});
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
            error_m=0;
            step=step-1;
            continue;
        end;

        %ref_path=ref_alg_path;
    end;

    
%     %start full search
%     thresh.std=0.9;
%     [merge_v stack_path_fine stack_alg_path_fine ref_alg_path_fine]=find_particles_below_threshold(class_st,thresh,stack_path,stack_alg_path,ref_alg_path);
%     disp(['Starting full Search  with ' num2str(round((length(merge_v)./size_stack(3) ).*10000)./100  ) '%  ('  num2str(length(merge_v)) ' of ' num2str(size_stack(3))  ') of all particles ']);
%     [align2d new_ref]=tom_av2_align_stack(stack_path_fine,ref_path,align2d,stack_alg_path_fine,ref_alg_path_fine,filter_param,paraell_param,iterations,demo,output_flag);
%     class_st=merge_class_st(class_st,align2d,merge_v);
% %     

    
    
     %command window print
    store.end=size(class_st,2); store.num_of_mesure=1;  store.mesured=0;
    [store]=tom_disp_estimated_time(store,'start','averaging');

    [h_out_stack_path,h_out_stack_name,h_out_stack_ext]=fileparts(stack_alg_path);
    h_out_stack_path=[h_out_stack_path '/'];
    
    % allocate space on Harddisk
    if (isempty(stack_alg_path)==0)
          %num=get_class_numbers(class_st,size_ref_stack(3));
          num=1; %to do !!!!!
          if (strcmp(output_flag,'sub_stacks'))
            for ii=1:size_ref_stack(3)
                tom_emwritec([h_out_stack_path h_out_stack_name '_' num2str(ii) h_out_stack_ext],[size_stack(1) size_stack(1) num(ii)],'new');
            end;
          end;
          tom_emwritec([stack_alg_path],[size_stack(1) size_stack(2) size_stack(3)],'new');
    end;

    out_stack_counter=1;
    for i=1:size(class_st,2)
        ref_nr=class_st(1,i);
        ccc=class_st(2,i);
        shift=[class_st(3,i) class_st(4,i)];
        rot=round(class_st(5,i));
        ratio=class_st(6,i);
        tmp_ratio_vect=class_st(7:18,i);
        ratio_vect=tmp_ratio_vect(find(tmp_ratio_vect));

        
        if  (ref_nr~=-1)

            %read particle
            part=tom_emreadc([stack_path],'subregion',[1 1 i],[size_stack(1)-1 size_stack(2)-1 0]);
            part=part.Value;

            if strcmp(correction_flag,'no Alignment')
                ref_nr=floor((class_st(1,i)-1)./(360./ang_samp))+1;
                rot=-((class_st(1,i)-1).*ang_samp);

                if (isempty(stack_alg_path)==0 || nargout < 2 || isempty(ref_alg_path) )
                    part=tom_rotate(tom_shift(part,shift),rot);
                    avg_st_im(:,:,ref_nr)=avg_st_im(:,:,ref_nr)+part;
                    avg_st_num(ref_nr)=avg_st_num(ref_nr)+1;
                    if (isempty(stack_alg_path)==0)
                        tom_emwritec([stack_alg_path],part,'subregion',[1 1 out_stack_counter],[size_stack(1) size_stack(1) 1]);
                        out_stack_counter=out_stack_counter+1;
                        if (strcmp(output_flag,'sub_stacks'))
                            tom_emwritec([h_out_stack_path h_out_stack_name '_' num2str(ref_nr) h_out_stack_ext],part,'subregion',[1 1 avg_st_num(ref_nr)],[size_stack(1) size_stack(1) 1]);
                        end;
                    end;
                end;
            end;

            if strcmp(correction_flag,'pre Alignment')
                if (isempty(stack_alg_path)==0 || nargout < 2 || isempty(ref_alg_path) )
                    part=tom_shift(tom_rotate(part,rot),[shift]);
                    avg_st_im(:,:,ref_nr)=avg_st_im(:,:,ref_nr)+part;
                    avg_st_num(ref_nr)=avg_st_num(ref_nr)+1;
                    if (isempty(stack_alg_path)==0)
                        tom_emwritec([stack_alg_path],part,'subregion',[1 1 out_stack_counter],[size_stack(1) size_stack(1) 1]);
                        out_stack_counter=out_stack_counter+1;
                        if (strcmp(output_flag,'sub_stacks'))
                            tom_emwritec([h_out_stack_path h_out_stack_name '_' num2str(ref_nr) h_out_stack_ext],part,'subregion',[1 1 avg_st_num(ref_nr)],[size_stack(1) size_stack(1) 1]);
                        end;
                    end;
                end;
            end;

            align2d(1,i).ref_class=ref_nr;
            warning off;
            align2d(1,i).shift.x = shift(1);
            align2d(1,i).shift.y = shift(2);
            warning on;
            align2d(1,i).angle=rot;
            align2d(1,i).ccc=ccc;
            align2d(1,i).ratio_min=ratio;
            align2d(1,i).ratio=ratio_vect;
        else
            %package got lost or classification failed
            align2d(1,i).ref_class=-1;
            align2d(1,i).shift=[];
            align2d(1,i).angle=0;
            align2d(1,i).ccc=0;
            align2d(1,i).ratio_min=ratio;
            align2d(1,i).ratio=ratio_vect;
        end;

        %command window print
        store.i=i;
        [store]=tom_disp_estimated_time(store,'estimate_time');
        %       [store]=tom_disp_estimated_time(store,'progress');
    end;


    %norm new ref stack according to the number particles
    for n=1:size_ref_stack(3)
        if (avg_st_num(n)~=0)
            avg_st_im(:,:,n)=avg_st_im(:,:,n)./avg_st_num(n);
        end;
    end;



    ref_path=ref_alg_path;
    if (isempty(ref_alg_path)==0 && (strcmp(ref_alg_path,'tmpXX')==0))
        ref_path=[ref_alg_path '_' num2str(step) ext_ref_alg_path];
        tom_emwrite([ref_alg_path '_' num2str(step) ext_ref_alg_path],avg_st_im);
    else
        tom_emwrite([pwd '/tmpXX'],avg_st_im);
        ref_path=[pwd '/tmpXX'];
    end;

    new_ref=avg_st_im;
    step=step+1;
    if (step > number_of_alg_steps)
        break;
    end;

    save('align_tmp','align2d');

    
    
    


end;


%final bookkeeping for writing sub alignment files
if (strcmp(output_flag,'sub_stacks'))
    align_old=align2d;
    for t=1:size(num,1)
        zz=1;
        clear('align2d');
        for tt=1:size(align_old,2)
            if (align_old(size(align_old,1),tt).ref_class ==t )
                align2d(1,zz)=align_old(size(align_old,1),tt);
                zz=zz+1;
            end;

        end;
        save([h_out_stack_path h_out_stack_name '_' num2str(t) '.mat'],'align2d');
    end;
    align2d=align_old;
end;



disp('end');

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
end;



function num=get_class_numbers(class_st,num_of_classes)

num=zeros(num_of_classes,1);

for i=1:size(class_st,2)
    ref_nr=class_st(1,i);
    num(ref_nr)=num(ref_nr)+1;
end;


function [merge_vect,stack_path_fine,stack_alg_path_fine,ref_alg_path_fine]=find_particles_below_threshold(class_st,thresh_st,stack_path,stack_alg_path,ref_alg_path)

all_min_ratios=class_st(6,:);
mean_min_ratio=mean(all_min_ratios);
std_min_ratio=std(all_min_ratios).*thresh_st.std;

merge_vect=find(class_st(6,:) < (mean_min_ratio-std_min_ratio) );

% %hack
%  disp('DANGER HACK ........random particles');
% merge_vect=round(rand(1,206).*959);
% % 
%  disp('DANGER HACK ........random particles');
% %end Hack


h=tom_reademheader(stack_path);
size_stack=h.Header.Size;

%create new filenames
[rootf namef extf]=fileparts(stack_path);
stack_path_fine=[rootf '/' namef '_fine'  extf];

[rootf namef extf]=fileparts(stack_alg_path);
stack_alg_path_fine=[rootf '/' namef '_fine'  extf];

[rootf namef extf]=fileparts(ref_alg_path); extf='.em';
ref_alg_path_fine=[rootf '/' namef '_fine'  extf];


sub_stack=zeros(size_stack(1),size_stack(2),length(merge_vect));
%build new stack
for i=1:length(merge_vect)
     part=tom_emreadc([stack_path],'subregion',[1 1 merge_vect(i)],[size_stack(1)-1 size_stack(2)-1 0]);
     part=part.Value;
     sub_stack(:,:,i)=part;
end;

tom_emwritec(stack_path_fine,sub_stack);



%db plot
tmp1(1:size(class_st,2))=mean_min_ratio;
tmp2(1:size(class_st,2))=mean_min_ratio-std(all_min_ratios).*thresh_st.std;
figure; plot(class_st(6,:)); hold on; plot(merge_vect,all_min_ratios(merge_vect),'ro'); 
plot(tmp1,'g-'); plot(tmp2,'r-');  
hold off;
%end db plot


% for i=1: size(class_st,2)
%     
%     
%     
% end;

disp('end');


function class_st=merge_class_st(class_st,align2d,merge_v)



for i=1:length(merge_v)
    class_st(1,merge_v(i))=align2d(1,i).ref_class;
    class_st(2,merge_v(i))=align2d(1,i).ccc;
    class_st(3,merge_v(i))=align2d(1,i).shift.x;
    class_st(4,merge_v(i))=align2d(1,i).shift.y;
    class_st(5,merge_v(i))=round(align2d(1,i).angle);
    class_st(6,merge_v(i))=0;
    class_st(7:18,merge_v(i))=0;
end;


disp('end');




