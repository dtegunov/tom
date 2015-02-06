function [vol_var,vol_avg,vol_var_even,vol_var_odd]=tom_xmippvariancemap3d_penczek_use_xmipp(path_doc,path_angles,num_of_models,proj_per_model,limit_cc,sym,max_shift,filter_kernel,mask,tmp_filepath,outputlevel)
%tom_xmippvariancemap3d_penczek_use_xmipp calculates a 3d variance Map 
%
%   [vol_var,var_stack,vol,stack]=tom_xmippvariancemap3d_penczek_use_xmipp(doc_filename)
%
%   calculates 3d variance volume by the Penzcek bootstrapping method
%   ...doc file is randomized  xmipp_angular_class_average
%   xmipp_reconstruct_fourier is called
%
%PARAMETERS
%
%  INPUT
%   path_doc          name of the xmipp *.doc file from ml3d or xmipp projMatch
%   path_angles       ref_angles file from xmipp use '' to generate (check examples or NOTE)
%   num_of_models     number of calculated 3d model
%   proj_per_model    number of particles per model
%   limit_cc          (10) in percent cc particles filter the lowes 10% are thrown away by default   
%   sym               (C1) symetry of model xmipp covention  (C1,Cn,dn ...)
%   max_shift         (% of length) max allowed shift in pix  
%   filter_kernel     (2) filter 
%   mask              (no mask) mask   
%   tmp_filepath      (rec_files4var) filepath for output of var files    
%   outputlevel       (1) debut outputlevel 
%
%
%  OUTPUT
%   vol_var           variance volume
%   vol_avg           average over all volumes
%   vol_var_even      even variance            
%   vol_var_odd       odd variance
%
%EXAMPLE
%  matlabpool open local 8; 
%  
% [vol_var,vol_avg,vol_var_even,vol_var_odd]=tom_xmippvariancemap3d_penczek _use_xmipp('ProjMatch/run1/Iter_10/Iter_10_current_angles.doc','ProjMatch/run1/Iter_10/ReferenceLibrary/ref_angles.doc',1000,4000)
%
% [vol_var,vol_avg,vol_var_even,vol_var_odd]=tom_xmippvariancemap3d_penczek_use_xmipp('Iter_10_current_angles.doc','ref_angles.doc',3000,4000);
%
% [vol_var,vol_avg,vol_var_even,vol_var_odd]=tom_xmippvariancemap3d_penczek_use_xmipp('Iter_10_current_angles.doc','ref_angles.doc',500,3000,10,'C1',10,2);
%
% %example without refangles
% 
% [vol_var,vol_avg,vol_var_even,vol_var_odd]=tom_xmippvariancemap3d_penczek_use_xmipp('Iter_10_current_angles.doc','',400,10000,10,'C1',10,2);
% 
%
%NOTE
%
% proj_per_model should be that the angular space is sampled once (26S ~4 000)
% num_of_models (26S ~10 000) ...check even odd correlation should be > 0.99
% 
% use for 64x64x64 a filter_kernel of 2
%
% ref_angles file (contains class nr and projangle) ...switch off clean up in xmipp protocols to get it
% path in proj-match is: ProjMatch/run1/Iter_10/ReferenceLibrary/ref_angles.doc
% or use tom_av2_xmipp_extract_ref_ang to obtain ref_angles.doc   
%
%!tail -5 ref_angles.doc
%   192 3    42.00000   90.00000    0.00000
%   193 3    30.00000   90.00000    0.00000
%   194 3   -36.00000   80.81005    0.00000
%   195 3   -30.00000   90.00000    0.00000
%   196 3   -42.00000   90.00000    0.00000
%
% use tail -f rec_files4var/logs/log.txt on unix shell 2 view xmipp output
%
%
%
%REFERENCES
%
%   Penczek Pawel A; Yang Chao; Frank Joachim; Spahn Christian M T
%   Estimation of variance in single-particle reconstruction using the bootstrap technique.
%   Journal of structural biology 2006;154(2):168-83.
%
%SEE ALSO
%   tom_av2_xmipp_extract_ref_ang,tom_xmippdocread,tom_av3_calc_variance,
%
%   created by FB (feat) Heinz Schenk 27/07/10
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


if (isempty(path_angles))
    disp('gen. ref_angles')
    tom_av2_xmipp_extract_ref_ang(path_doc,'tmp_angXXX.doc');
    path_angles='tmp_angXXX.doc';
    disp('saving 2 tmp_angXXX.doc');
end;


%read first entry
[a b]=unix(['head -2 ' path_doc ' | tail -1']);
name_first=strtrim(strrep(b(1:end-1),' ; ',''));
im_tmp=tom_spiderread(name_first);

if nargin < 3
    num_of_models=3000;
end;

if nargin < 4
    proj_per_model=5000;
end;

if nargin < 5
   limit_cc=10;
end;

if nargin < 6
    sym='C1';  
end;

if nargin < 7
    max_shift=10;
end;

if nargin < 8
   filter_kernel=2;
end;

if nargin < 9 || strcmp(mask,'no_mask')
    mask=ones(size(im_tmp.Value,1),size(im_tmp.Value,1),size(im_tmp.Value,1));
end;

if nargin < 10
    tmp_filepath=[pwd '/rec_files4var'];
end;

if nargin < 11
    outputlevel=1;
end;


if exist(tmp_filepath,'dir')==0
    mkdir(tmp_filepath);
    unix(['chmod ugo+rwx ' tmp_filepath]);
end;

warning off;
mkdir([tmp_filepath '/logs/']);
unix(['chmod ugo+rwx ' tmp_filepath '/logs/']);
log_p=[tmp_filepath '/logs/log.txt'];

mkdir([tmp_filepath '/classes/']);
unix(['chmod ugo+rwx ' tmp_filepath '/classes/']);

mkdir([tmp_filepath '/vols/']);
unix(['chmod ugo+rwx ' tmp_filepath '/vols/']);
warning on;


start_model=1;
dd=dir([tmp_filepath '/vols/*.spi']);

if (isempty(dd)==0)
    answer=questdlg('Volumes found. Append ?');
    if (strcmp(answer,'Cancel'))
        vol_var='';
        vol_avg='';
        vol_var_odd='';
        vol_var_even='';
        return;
    end;
    if (strcmp(answer,'Yes'))
        for i=1:length(dd)
            [a rest]=strtok(dd(i).name,'_');
            rest=strrep(rest,'_','');
            num(i)=str2double(strtok(rest,'.'));
        end;
        start_model=max(num);
    end;
    if (strcmp(answer,'No'))
        disp(['Cleaning up: rm -Rf ' tmp_filepath ]);
        unix(['rm -Rf ' tmp_filepath]);
    end;
end;



%unix(['echo ' datestr(now) ' > ' log_p]);





[a program_name_cl_avg]=unix('which xmipp_angular_class_average');
program_name_cl_avg=program_name_cl_avg(1:end-1);
program_name_cl_avg=deblank(strtrim(program_name_cl_avg));

max_sh_tmp=round((max_shift./100)*size(im_tmp.Value,1));
paramstr_cl_avg=['-limit0 -1 -limitR ' num2str(limit_cc) ' -iter 1 -Ri 0 -Ro ' num2str(round(size(im_tmp.Value,1)./2)) ' -max_shift ' num2str(max_sh_tmp) ' -max_shift_change ' num2str(max_sh_tmp-1) ' -max_psi_change 1'];


[a program_name_rec_four]=unix('which xmipp_reconstruct_fourier');
program_name_rec_four=program_name_rec_four(1:end-1);
program_name_rec_four=deblank(strtrim(program_name_rec_four));

clean_up_flag=0;


disp('Volume Reconstruction started..');
dsp_param(start_model,num_of_models,proj_per_model,limit_cc,sym,max_sh_tmp,filter_kernel,mask,tmp_filepath,clean_up_flag,log_p);

parfor i=start_model:start_model+(num_of_models-1)
%for i=start_model:start_model+(num_of_models-1)

    try
        
        fpp=fopen(log_p,'a+');
        tic;
        disp(['Processing Vol Nr: ' num2str(i)]);
        pr_cl=['/ProjMatchClasses' num2str(i)];
        
        warning off;
        mkdir([tmp_filepath '/classes' pr_cl]);
        warning on;
        
        [a b]=unix(['chmod ugo+rwx ' tmp_filepath '/classes' pr_cl]);
        
        create_random_doc(path_doc,pr_cl,tmp_filepath,250,proj_per_model);
        
        lib_str=path_angles;
        doc_str=[tmp_filepath '/classes' pr_cl '/input_doc_cut.doc'];
        output_str=[tmp_filepath '/classes' pr_cl '/proj_match'];
        
        %call=[program_name_cl_avg ' -i ' doc_str ' -lib ' lib_str ' -dont_write_selfiles '  paramstr_cl_avg ' -o ' output_str ' >> ' log_p ' 2>&1'];
        
       % call=[program_name_cl_avg ' -i ' doc_str ' -lib ' lib_str ' -dont_write_selfiles '  paramstr_cl_avg ' -o ' output_str ' > /dev/null 2>&1'];
        
        %fix for speed 2>&1 takes ages since update !!!!
        call=[program_name_cl_avg ' -i ' doc_str ' -lib ' lib_str ' -dont_write_selfiles '  paramstr_cl_avg ' -o ' output_str ];
        disp(call);
        
        %fix 2 pipe 2 /dev/null
        %call=[program_name_cl_avg ' -i ' doc_str ' -lib ' lib_str ' -dont_write_selfiles '  paramstr_cl_avg ' -o ' output_str ' > test.txt 2>&1'];
        
        error_out=unix(call);
        if (error_out > 0)
            error('error building classes');
        end;
        
        %build up reconstruction sel
        xmp_str=[tmp_filepath '/classes' pr_cl '/proj_match'];
        output_str=[tmp_filepath  '/classes' pr_cl '/reconstruction.sel'];
        program_name='ls';
        call=[program_name ' ' xmp_str '_class*.xmp  | awk ''' '{print $1 " 1"}'' > '  output_str];
        error_out=unix(call);
        if (error_out > 0)
            error('error building .sel');
        end;
        
        rec_sel_str=[tmp_filepath '/classes' pr_cl '/reconstruction.sel'];
        output_str=[tmp_filepath '/vols' '/rec_' num2str(i) '.spi'];
        %call=[program_name_rec_four ' -i ' rec_sel_str ' -o ' output_str ' -sym ' sym ' -thr 1 -weight  -max_resolution 0.36562528 >  /dev/null 2>&1'];
        call=[program_name_rec_four ' -i ' rec_sel_str ' -o ' output_str ' -sym ' sym ' -thr 1 -weight  -max_resolution 0.49'];
        disp(call);
        
        error_out=unix(call);
        if (error_out > 0)
            error(['error reconstruction Vol Nr: ' num2str(i)']);
            fprintf(fp,'error reconstruction Vol Nr:  %s \n',num2str(i));
        end;
        disp(['Vol Nr: ' num2str(i) ' done!']);
        fprintf(fpp,'Vol Nr: %s  done!\n',num2str(i));
        toc;
        fclose(fpp);
    catch ME
        disp(['Error reconstruction Vol Nr '  num2str(i)]);
        disp(ME.message);
    end;
    
end;

[vol_var vol_avg vol_var_even vol_var_odd]=tom_av3_calc_variance({[tmp_filepath '/vols/rec_']},{'.spi'},'all_you_can',mask,filter_kernel,outputlevel);

if (clean_up_flag==1)
    unix(['rm -R ' tmp_filepath '/classes']);
    disp(['clean up: rm -R ' tmp_filepath '/classes']);
end;

warning off;
mkdir([tmp_filepath '/var_vols/']);
warning on;


%build model containing all unfiltered


tom_spiderwrite([tmp_filepath '/var_vols/vol_var.spi'],vol_var);
tom_spiderwrite([tmp_filepath '/var_vols/vol_even.spi'],vol_var_even);
tom_spiderwrite([tmp_filepath '/var_vols/vol_odd.spi'],vol_var_odd);
tom_spiderwrite([tmp_filepath '/var_vols/vol_avg.spi'],vol_avg);
unix(['chmod -R ugo+rwx ' tmp_filepath '/var_vols &']);


function dsp_param(start_model,num_of_models,proj_per_model,limit_cc,sym,max_shift,filter_kernel,mask,tmp_filepath,clean_up_flag,log_p)

disp('---------->');
disp(' ');
disp(['Start Nr: ' num2str(start_model)]);
disp(['Num of models: ' num2str(num_of_models)]);
disp(['Proj per model: ' num2str(proj_per_model)]);
disp(['Limit cc: ' num2str(limit_cc)]);
disp(['Sym: ' sym]);
disp(['Max shift: ' num2str(max_shift)]);
disp(['Filter kernel: ' num2str(filter_kernel)]);
disp(['Mask std2: ' num2str(std2(mask)) ]);
disp(['Tmp filepath: ' tmp_filepath]);
disp(['Clean up: (db flag)' num2str(clean_up_flag)]);
disp(' ');
disp('<----------');


fp=fopen(log_p,'wt');
fprintf(fp,['%s \n'],datestr(now));
fprintf(fp,'----------> \n');
fprintf(fp,'Start Nr:  %s \n',num2str(start_model));
fprintf(fp,'Num of models:  %s \n',num2str(num_of_models));
fprintf(fp,'Proj per model:  %s \n',num2str(proj_per_model));
fprintf(fp,'Limit cc:  %s \n',num2str(limit_cc));
fprintf(fp,'Sym:   %s \n',sym);
fprintf(fp,'Max shift:  %s \n',num2str(max_shift));
fprintf(fp,'Filter kernel:  %s \n',num2str(filter_kernel));
fprintf(fp,'Mask std2:  %s \n',num2str(std2(mask)) );
fprintf(fp,'Tmp filepath:  %s \n',tmp_filepath);
fprintf(fp,'Clean up: (db flag)  %s \n',num2str(clean_up_flag));
fprintf(fp,'<---------- \n');
fclose(fp);

disp(' ');

function create_random_doc(path_doc,pl,tmp_filepath,package_size,proj_per_model)

[a num_of_particles]=unix(['wc -l ' path_doc]);
num_of_particles=str2num(strtok(num_of_particles,' '));
num_of_particles=(num_of_particles-1)./2;

if (package_size > proj_per_model)
    package_size=proj_per_model;
end;

used_part=randperm(num_of_particles);

% if (rand(1)>0.5)
%     used_part=1:615;
% else
%     used_part=615:1200;
% end;

packages=tom_calc_packages(round(proj_per_model./package_size),proj_per_model);

doc_str=path_doc;
output_str=[tmp_filepath '/classes' pl '/input_doc_cut.doc'];
unix(['head -1 ' doc_str ' > ' output_str]);
tic;

for ii=1:size(packages,1)
    
    pattern=['"' num2str(used_part(packages(ii,1))) ' 8|'];
    
    for i=packages(ii,1)+1:packages(ii,2)
        pattern=[pattern num2str(used_part(i)) ' 8|'];
    end;
    pattern=[pattern(1:end-1) '"'];
    
    
    %cut input doc file
    program_name='egrep';
    
    call=[program_name ' -w -B1 ' pattern ' ' doc_str ' | sed ''/--/d'''  ' >> ' output_str];
    unix(call);
end;



