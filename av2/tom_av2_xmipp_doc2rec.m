function vol=tom_av2_xmipp_doc2rec(path_doc,sym,discard_p,max_res,np,outputfilename,verbose,outputfolder,clean_up,mpi_prefix,ref_angles,add_param_cl_av,add_param_rec,fsc_param)
%TOM_AV2_XMIPP_DOC2REC calculates 3d rec out of doc file
%   
%
%  vol=tom_av2_xmipp_doc2rec(path_doc,sym,discard_p,max_res,np,outputfilename,verbose,outputfolder,clean_up,mpi_prefix,ref_angles,add_param_cl_av,add_param_rec,fsc_param)
%
%  TOM_AV2_XMIPP_DOC2REC creates 3d volume from a give .doc file using
%  xmipp
%  
%
%PARAMETERS
%
%  INPUT
%   path_doc            xmipp doc filename
%   sym                 (C1) symmetry C1..Cn 
%   discard_p           (10) dispcard particles according 2 cc in percent  
%                            10 means the best 90% are used
%   max_res             (0.48) filter cut off for reconstrution in digital freq
%                             0.5   means no filter
%                             0.25  means half nyquist   
%   np                  (7) number of processors   
%   outputfilename      ('') filename 2 write the reconstructed volume 
%                            use '' for no output (format acc 2 extension use .spi 4 spider and .em for em) 
%   verbose             (1)  verbose level use 0 to switch off xmipp output
%   outputfolder        ('./rec_xmipp') folder 4 xmipp output   
%   clean_up            (1) flag for deleting the outputfolder                         
%   mpi_prefix          ('mpirun -np')
%   ref_angles          path to ref_angles file (speeds up a little) 
%                       (use '' to generat on the fly) 
%   add_param_cl_av     ('-dont_write_selfiles  -limit0 -1 ') 
%   add_param_rec       ('-thr 1 -weight  -pad_proj 2.0 -pad_vol 2.0 -mpi_job_size 7')  
%   fsc_param           ([1 50]) vector with pixels size in ang and max res 
%                       
%
%  OUTPUT
%
%  vol                 reconstructed volume
%
%EXAMPLE
%     
%  %get volume,discard 10% and filter at 0.3
%  vol=tom_av2_xmipp_doc2rec('Iter_8_current_angles.doc','C2',10,0.3,8);
%
%  %write volume discard 40% and filter at 0.3 and machinefile
%  tom_av2_xmipp_doc2rec('Iter_8_current_angles.doc','C2',40,0.3,8,'my_rec.em',1,'./rec_tmp',1,'mpirun -machinefile  ~/machinefile.dat -np');
%
%
%  %reduce pad to 1.2 ...speed up
%  add_param_cl_av='-dont_write_selfiles  -limit0 -1 ';
%  add_param_rec='-thr 1 -weight  -pad_proj 1.2 -pad_vol 1.2 -mpi_job_size 7';  
%  tom_av2_xmipp_doc2rec('Iter_8_current_angles.doc','C2',40,0.3,8,'my_rec.em',0,'./rec_tmp',0,'mpirun -np','',add_param_cl_av,add_param_rec);
% 
%  %reduce pad to 1.2 and calc fsc with pixelsize of 7.2 Angstroem and maxres of 25 (25 means everything smaller than 16 is set to 0)
%  add_param_cl_av='-dont_write_selfiles  -limit0 -1 -split';
%  add_param_rec='-thr 1 -weight  -pad_proj 1.2 -pad_vol 1.2 -mpi_job_size 7';  
%  tom_av2_xmipp_doc2rec('Iter_8_current_angles.doc','C2',40,0.3,8,'my_rec.em',0,'./rec_tmp',0,'mpirun -np','',add_param_cl_av,add_param_rec,[7.2 16]);
%
% 
%
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_xmipp_extract_ref_ang,tom_av2_em_classify3d
%
%   created by FB 08/09/09
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


if nargin < 2
    sym='C1';  
end;

if nargin < 3
   discard_p=10;  
end;

if nargin < 4
    max_res=0.48;  
end;

if nargin < 5
    np=7;  
end;

if nargin < 6
    outputfilename='';
end;

if nargin < 7
    verbose=1;
end;    

if nargin < 8
   outputfolder='./rec_xmipp';
end;    

if nargin < 9
   clean_up=1;
end;

if nargin < 10
   mpi_prefix='mpirun -np';
end;

if nargin < 11
   ref_angles='';
end;

if nargin < 12
   add_param_cl_av=' -dont_write_selfiles';
end;

if nargin < 13
   add_param_rec='';
end;

if nargin < 14
   fsc_param=[1 2];
end;

%number of reconstructed volumes
if (isempty(findstr(add_param_cl_av,'-split')) )
    sp_nr=1;
else
    sp_nr=3;
end;

if (np==1)
    mpi_prefix_str='';
else
    mpi_prefix_str=[mpi_prefix ' ' num2str(np)];
end;

%generate output folder struct
[log_p]=gen_fold_struct(outputfolder);


%generate ref-angles from input doc file
tic;
if (isempty(ref_angles))
     disp('generating refangles ...');
     path_angles=[outputfolder '/refangles.doc'];
     tom_av2_xmipp_extract_ref_ang(path_doc,path_angles);
     disp('done!');
end;


if (verbose==0)
    console_output_str=['>> ' log_p ' 2>&1'];
end;

if (verbose==1)
    console_output_str='';
end;

    
%gen log file
unix(['echo ' datestr(now) ' > ' log_p]);

disp('Volume Reconstruction started..');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate class averages
%%%%%%%%%%%%%%%%%%%%%%%%%%%

program_name=get_abs_prog_name('class_avg',np);
output_str=[outputfolder '/classes/ProjMatchClasses/proj_match'];
parameters=[' -i ' path_doc ' -lib ' path_angles ' -o ' output_str ' -limitR ' num2str(discard_p) ' ' add_param_cl_av];


call=[mpi_prefix_str ' ' program_name ' ' parameters ' ' console_output_str];
write_log(log_p,call);
disp(call);
error_out=unix(call);

if (error_out > 0)
    error('error building classes');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate reconstruction sel
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:sp_nr
    
    if (i==1)
        xmp_str=[outputfolder '/classes/ProjMatchClasses/proj_match'];
        output_str=[outputfolder '/classes/ProjMatchClasses/reconstruction.sel'];
    end;
    if (i==2)
        xmp_str=[outputfolder '/classes/ProjMatchClasses/proj_match_split_1'];
        output_str=[outputfolder '/classes/ProjMatchClasses/reconstruction_sp1.sel'];
    end; 
    if (i==3)
        xmp_str=[outputfolder '/classes/ProjMatchClasses/proj_match_split_2'];
        output_str=[outputfolder '/classes/ProjMatchClasses/reconstruction_sp2.sel'];
    end;
    program_name='ls';
    call=[program_name ' ' xmp_str '_class*.xmp  | awk ''' '{print $1 " 1"}'' > '  output_str];
    write_log(log_p,call);
    disp(call);
    error_out=unix(call);
    if (error_out > 0)
        error('error building .sel');
    end;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reconstruct volume(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:sp_nr
    if (i==1)
        output_str=[outputfolder '/vols/rec.spi'];
        reconstruction_sel_str=[outputfolder '/classes/ProjMatchClasses/reconstruction.sel'];
    end;
    if (i==2)
        output_str=[outputfolder '/vols/rec_sp1.spi'];
        reconstruction_sel_str=[outputfolder '/classes/ProjMatchClasses/reconstruction_sp1.sel'];
    end;
    if (i==3)
        output_str=[outputfolder '/vols/rec_sp2.spi'];
        reconstruction_sel_str=[outputfolder '/classes/ProjMatchClasses/reconstruction_sp2.sel'];
    end;
    program_name=get_abs_prog_name('rec',np);
    parameters=[' -i ' reconstruction_sel_str ' -o ' output_str ' -sym ' sym ' -max_resolution ' num2str(max_res) ' ' add_param_rec];
    
    call=[mpi_prefix_str ' ' program_name ' ' parameters ' ' console_output_str];
    write_log(log_p,call);
    disp(call);
    error_out=unix(call);
    if (error_out > 0)
        error('error reconstruction volume');
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calc fsc
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (sp_nr>1)
    output_str=[outputfolder '/fsc/resolution.fsc'];
    program_name=get_abs_prog_name('fsc',np);
    parameters=['-ref ' outputfolder '/vols/rec_sp1.spi' ' -i ' outputfolder '/vols/rec_sp2.spi' ' -sam ' num2str(fsc_param(1))  ' -max_sam ' num2str(fsc_param(2)) ];
    
    call=[program_name ' ' parameters ' ' console_output_str];
    write_log(log_p,call);
    disp(call);
    error_out=unix(call);
    if (error_out > 0)
        error('error calculating fsc');
    end;
    unix(['mv ' outputfolder '/vols/rec_sp2.spi.frc ' output_str ]);
end;

disp(['Reconstruction done!']);

if (nargout > 0 || isempty(outputfilename)==0)
    vol=tom_spiderread([outputfolder '/vols/rec.spi']);
    vol=vol.Value;
end;

if (isempty(outputfilename)==0)
    [a b c]=fileparts(outputfilename);
    if (strcmp(c,'.em'))
        tom_emwrite(outputfilename,vol);
    end;
    if (strcmp(c,'.spi'))
        tom_spiderwrite(outputfilename,vol);
    end;
    if (strcmp(c,'.hdf'))
        tom_eman2_write(outputfilename,vol,'volume');
    end;
end;

if (clean_up>0)
    disp('cleaning up ...');
    [a b]=unix(['rm -R ' outputfolder]);
    if (clean_up>1)
        [a b]=unix(['rm -R ' outputfolder]);
    end;
    disp('done!');
end;
toc;

function [log_p]=gen_fold_struct(tmp_filepath)


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

mkdir([tmp_filepath '/fsc/']);
unix(['chmod ugo+rwx ' tmp_filepath '/fsc/']);

pr_cl='/ProjMatchClasses';

mkdir([tmp_filepath '/classes' pr_cl]);
warning on;

[a b]=unix(['chmod ugo+rwx ' tmp_filepath '/classes' pr_cl]);


function name_abs=get_abs_prog_name(prog,np)

if (strcmp(prog,'class_avg'))
    if (np == 1)
        [a name_abs]=unix('which xmipp_angular_class_average');
        mpi_prefix='';
    else
        [a name_abs]=unix('which xmipp_mpi_angular_class_average');
    end;
end;

if (strcmp(prog,'rec'))
    if (np == 1)
        [a name_abs]=unix('which xmipp_reconstruct_fourier');
    else
        [a name_abs]=unix('which xmipp_mpi_reconstruct_fourier');
    end;
end;

if (strcmp(prog,'fsc'))
    [a name_abs]=unix('which xmipp_resolution_fsc');
end;

name_abs=deblank(strtrim(name_abs(1:end-1)));

function write_log(log,str)

fid=fopen(log,'a');
if (fid==-1)
    error('cannot write log!!');
end;
fprintf(fid,'%s\n',str);
fclose(fid);
