function avg=tom_k2_average_stack(wkFold,wkImg,target_name,cut_size,sub_idx,head_st,cc_thr,bin,param4alg,verbose)
%tom_av2_k2_align_stack aligns image stack
%
%    avg=tom_av2_k2_align_stack(wk,max_shift)
%
%PARAMETERS
%
%  INPUT
%   wkFold            wildcard for folders
%   wkImg             wildcard for images 
%   target_name       basename 2 put the averages of the stacks
%   cut_size          ('') size 2 make images qudratic  
%   sub_idx           ('all') index for using only a subset of all frames
%   head_st           ('') em-header struct ...check example 
%   cc_thr            (-1) threshold for removing frames by cc
%   bin               (0) binning 
%   param4alg         ('') structure containing the parameters 4
%                         tom_k2_align_stack (use '' to switch off alignment)
%   verbose           (1) verbose flag    
%
%EXAMPLE
%   
%   %generate example data
%   mkdir low_stacks
%   mkdir low_stacks/stack1
%   mkdir low_stacks/stack2 
%
%   all_sh{1}=[-50 80];
%   all_sh{2}=[32 -103];
%   all_sh{3}=[102 98]; 
%
%   img{1}=tom_spheremask(ones(2048,2060),400,5)+tom_spheremask(ones(2048,2060),200,2,[1025 1025 1]-[400 400 0]);   
%   img{2}=tom_spheremask(ones(2048,2060),400,5)+tom_spheremask(ones(2048,2060),200,2,[1025 1025 1]+[400 400 0])+tom_spheremask(ones(2048,2060),150,2,[1025 1025 1]-[305 305 0]);
%   
%   for i=1:2
%      for ii=1:3   
%           tom_mrcwrite(tom_move(img{i},all_sh{ii})+(0.4.*randn(2048,2060)),'name',['low_stacks/stack' num2str(i) '/K' sprintf('%03d',ii) '.mrc']);      
%      end; 
%   end;
%   
%   tmp=tom_emheader(ones(2048,2048)); %use the final size after cutting!!
%   tmp.Header.Defocus=-3.*10000;
%   tmp.Header.Objectpixelsize=1.3;
%   tmp.Header.Microscope='Titan2';
%   tmp.Header.Voltage=300000;
%   tmp.Header.Cs=2;   
%   
%   matlabpool open local 4; 
%  
%   %example 4 initial avg  (no alignment used for sorting and viewing)
%   tom_k2_average_stack('low_stacks/stack*','K*.mrc','low_ini/low_',[2048 2048],'all',tmp.Header);
% 
%   %example with alignment 
%   param4alg.cut_size='square';
%   param4alg.max_shift=200;
%   param4alg.sl_window=1;
%   param4alg.sl_weights=1;
%   param4alg.filt=[4 4 3 3 2 2 1 1];
%   param4alg.num_iter=15;
%   param4alg.ref='sum';
%   param4alg.memory='read_all_frames';
%   param4alg.remove_cc_mid=1;
%   param4alg.remove_auto_corr=1;
%
%   tom_k2_average_stack('low_stacks/stack*','K*.mrc','low_alg/low_',[2048 2048],'all',tmp.Header,-1,0,param4alg); 
%
%  %example 4 using only the first 2 images with global alignment
%  tom_k2_average_stack('low_stacks/stack*','align.log','low_algsub/low_',[2048 2048],[1:2],tmp.Header); 
%  
%
%SEE ALSO
%   tom_k2_average_stack
%
%   created by FB 03/18/13
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


if (nargin < 4)
    cut_size='';
end;

if (nargin < 5)
    sub_idx='all';
end;

if (nargin < 6)
    head_st='';
end;

if (nargin < 7)
    cc_thr=-1;
end;

if (nargin < 8)
    bin=0;
end;

if (nargin < 9)
    param4alg='';
end;

if (nargin < 10)
    verbose=1;
end;


warning off;
mkdir(fileparts(target_name));
warning onf;

%initialize log
log_name=[fileparts(target_name) '_log.txt'];
[a b]=unix(['rm ' log_name]);

%scan stack folders
all_folds=get_folders(wkFold);
[sz_img,all_stack_len]=get_stack_information(all_folds,wkImg);
if (std(all_stack_len)==0)
    disp_it(['Found ' num2str(length(all_folds)) ' folders  with ' num2str(all_stack_len(1)) ' files '],verbose,'');
else
    disp_it(['Found ' num2str(length(all_folds)) ' folders with ' num2str(min(all_stack_len))  ' to  '  num2str(max(all_stack_len)) ' files '],verbose,''); 
end;

%allocating memory for output avg
n_outputargs=nargout;
if (n_outputargs > 0)
    try
        avg=zeros(sz_img(1),sz_img(2),length(all_folds));
    catch Me
        disp('Memory allocation problem');
    end;
        
end;

%initial screen output
disp_it(' ',verbose,'');
disp_it('Averaging starts ...',verbose,'');
disp_it(' ',verbose,'');

%working (par)loop over all folders
parfor i=1:length(all_folds)
    disp_it(['Processing ' all_folds{i}],verbose,'');
    [framesInFold,shifts,cc]=get_frames(all_folds{i},wkImg);
    if (strcmp(sub_idx,'all'))
        tmp_idx=1:all_stack_len(i);
    else
        tmp_idx=sub_idx;
    end;
    avg_tmp=proj_st(framesInFold(tmp_idx),sz_img,shifts,cc,cc_thr,verbose);
    if (n_outputargs > 0)
        avg(:,:,i)=avg_tmp;
    end;
    name_out=[target_name num2str(i) '.em'];
    if (isempty(head_st)==0)
        avg_tmp=tom_emheader(avg_tmp);
        avg_tmp.Header=head_st;
    end;
    if (isempty(cut_size)==0)
        avg_tmp.Value=tom_cut_out(avg_tmp.Value,'center',cut_size);
    end;
    
    if (bin>0)
        avg_tmp.Value=tom_bin(avg_tmp.Value,bin);
        avg_tmp.Header.Objectpixelsize=avg_tmp.Header.Objectpixelsize./(2.^bin);
    end;
    if (isempty(param4alg)==0)
        [avg_tmp.Value,alg_param]=tom_k2_align_stack(framesInFold(tmp_idx),param4alg.cut_size,param4alg.max_shift,...
                                                     param4alg.sl_window,param4alg.sl_weights,param4alg.filt,...
                                                     param4alg.num_iter,avg_tmp.Value,param4alg.memory,...
                                                     param4alg.remove_cc_mid,param4alg.remove_auto_corr,verbose); 
        gen_alg_log(alg_param,framesInFold(tmp_idx),all_folds{i}); 
    end;
    if (isstruct(avg_tmp))
        avg_tmp.Value=tom_xraycorrect2(avg_tmp.Value);
    else
        avg_tmp=tom_xraycorrect2(avg_tmp);
    end;
    tom_emwrite(name_out,avg_tmp);
    disp_it([all_folds{i} ' ' name_out],0,log_name);
    disp_it(['to ' target_name num2str(i) '.em'],verbose,'');
    disp_it([' '],verbose,'');
end;

function gen_alg_log(alg_param,framesInFold,fold)
fid=fopen([fold '/align.log' ],'wt');
for i=1:length(framesInFold)
    fprintf(fid,'%s %f %f %f\n',framesInFold{i},alg_param(i,1),alg_param(i,2),alg_param(i,3));
end;
fclose(fid);





function avg=proj_st(names_stack,sz_img,shifts,cc,cc_thr,verbose)

avg=zeros(sz_img);
for i=1:length(names_stack);
    img_tmp=tom_mrcread(names_stack{i});
    if (cc > cc_thr)
        avg=avg+tom_shift(single(img_tmp.Value),shifts(i,:));
        disp_it(['   adding: ' names_stack{i} ' shift: ' num2str(shifts(i,:)) ],verbose,'');
    else
        disp_it(['   skipped cc < thr: ' names_stack{i}],verbose,'');
    end;
end;

function [names,shifts,cc]=get_frames(fold,wkImg)

dd=dir([fold '/' wkImg]);
if (tom_ismrcfile([fold '/' dd(1).name])==1)
    shifts=zeros(length(dd),2);
    cc=ones(length(dd),1).*1;
    for i=1:length(dd);
        names{i}=[fold '/' dd(i).name];
    end;
else
    logf=importdata([fold '/' dd(1).name]);
    names=logf.textdata;
    shifts=logf.data(:,1:2);
    cc=logf.data(:,3);
end;


function all_folds=get_folders(wkFold)
%resolve wk
d=dir(wkFold);
for i=1:length(d)
    if (d(i).isdir)
        tmp_part=fileparts(wkFold);
        if (isempty(tmp_part))
            all_folds{i}=d(i).name;
        else
            all_folds{i}=[tmp_part '/' d(i).name];
        end;
    end;
end;


function [sz_vol,all_len]=get_stack_information(folds,wkImg)


dd=dir([folds{1} '/' wkImg]);
name=[folds{1} '/' dd(1).name];
if (tom_ismrcfile(name))
    t_tmp=tom_mrcread(name);
else
    logf=importdata(name);
    t_tmp=tom_mrcread(logf.textdata{1});
end;
sz_vol=size(t_tmp.Value);


for i=1:length(folds)
    
    dd=dir([folds{i} '/' wkImg]);
    if (tom_ismrcfile([folds{i} '/' dd(1).name]) )
        all_len(i)=length(dd); 
    else
        logf=importdata([folds{i} '/' dd(1).name]);
        all_len(i)=length(logf.textdata);
    end;
end;

function disp_it(message,verbose,log)

if (verbose==1)
    disp(message);
end;

if (isempty(log)==0)
    fid_log=fopen(log,'a');
    fprintf(fid_log,'%s\n',message);
    fclose(fid_log);
end;

