function [avg,alg_param,all_avg]=tom_k2_align_stack(frames,cut_size,max_shift,sl_window,sl_weights,filt,num_iter,ref,memory,remove_cc_mid,remove_auto_corr,verbose)
%tom_av2_k2_align_stack aligns image stack
%
%   avg=tom_k2_align_stack(frames,max_shift,filt,num_iter)
%
%PARAMETERS
%
%  INPUT
%   frames                  matlab cell of frames
%   cut_size                ('square') size 2 cut the images                             
%   max_shift               (300) max shift 4 masking cc in pixel (use Inf 2 switch off)
%   sl_window               (5) sliding window 5 avg over time 
%   sl_weights              (0.6 0.8 1 0.8 0.6) weights 4 frames of sl.window 
%                                                (..middle should be biggest!)
%   filt                    (4) filter 4 crosscorrelation  
%   num_iter                (4) number of alg iter
%   ref                     ('sum') initial reference  
%   memory                  ('read_all_frames') memory buffer behaviour 
%                            or 'single frame' to have only 1 frame in memory
%   remove_cc_mid           (1) set mid of cc func to -1
%   remove_auto_corr        (1) subtract the act particel from mean 
%                               to avoid autocc
%   verbose                 (1) verbose flag 
%
%  OUTPUT
%    avg                  final average 
%    alg_param            shifts and cc-values for ervery frame
%    all_avg              stack with average 4 each iteration   
% 
%EXAMPLE
%
%  mkdir low_stacks
%  mkdir low_stacks/stack1
% 
%  all_sh{1}=[60 133];
%  all_sh{2}=[-88 -78];
%  all_sh{3}=[77 123]; 
%  all_sh{4}=[-40 61];   
%  all_sh{5}=[22 11]; 
%  all_sh{6}=[0 0]; 
%
%  img{1}=tom_spheremask(ones(2048,2060),400,5)+tom_spheremask(ones(2048,2060),200,2,[1025 1025 1]-[400 400 0]);   
%
%  for ii=1:6   
%     tom_mrcwrite(tom_move(img{1},all_sh{ii})+(0.8.*randn(2048,2060)),'name',['low_stacks/stack' num2str(1) '/K' sprintf('%03d',ii) '.mrc']);      
%     frames{ii}=['low_stacks/stack' num2str(1) '/K' sprintf('%03d',ii) '.mrc'];
%  end; 
%
%  
% disp('test vs reference');
% [avg2gt,alg_param2gt,allavg2gt]=tom_k2_align_stack(frames,'square',200,1,1,1,1,img{1});
%
% disp('check weighted sliding window ');  
% [avg2gt,alg_param2gt,allavg2gt]=tom_k2_align_stack(frames,'square',200,3,[0.3 1 0.3],1,1,img{1});
%
% disp('checking optimisation');
% [avg,alg_param,allavg]=tom_k2_align_stack(frames,'square',200,1,1,[4 4 3 3 3 2 2 2 1 1 1],20,'sum'); 
%
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

if (nargin < 2)
    cut_size='square';
end;

if (nargin < 3)
    max_shift=300;
end;

if (nargin < 4)
    sl_window=5;
end;

if (nargin < 5)
    sl_weights=[0.6 0.8 1 0.8 0.6];
end;

if (nargin < 6)
    filt=[3 3 2 2 1 1];
end;

if (isempty(find(filt==0))==0)
    idx_tmp=find(filt==0);
    filt(idx_tmp)=1;
end;

if (nargin < 7)
    num_iter=6;
end;

if (nargin < 8)
    ref='sum';
end;

if (nargin < 9)
    memory='read_all_frames';
end;

if (nargin < 10)
    remove_cc_mid=1;
end;

if (nargin < 11)
    remove_auto_corr=1;
end;

if (nargin < 11)
    verbose=1;
end;


num_out_arg=nargout;

%extract size from dataset
tmp=tom_mrcread(frames{1});
sz_img=size(tmp.Value);
clear('tmp');
if (strcmp(cut_size,'square'))
    cut_size=[min(sz_img) min(sz_img)];
end;
mid=floor(cut_size./2)+1;
num_of_frames=length(frames);

%extend filt struct
if (num_iter > length(filt) )
    filt=[filt ones(1,num_iter-length(filt)).*filt(end)];
end;

%buld cc-mask index
if (isinf(max_shift))
    mask_cc_trans=ones(cut_size);
else
    mask_cc_trans=tom_spheremask(ones(cut_size),max_shift);
end;
if (remove_cc_mid==1)
    mask_cc_trans(mid(1),mid(2))=0;
end;
idx_ccf=find(mask_cc_trans==0);

%generate reference
disp_it('generating inital reference',verbose);
if (strcmp(ref,'sum'))
    if (strcmp(memory,'read_all_frames'))
        [ref,frames]=proj_st(frames,cut_size);
    else
        ref=proj_st(frames,cut_size);
    end;
    act_remove_auto_corr=remove_auto_corr;
else
    if (strcmp(memory,'read_all_frames'))
        [~,frames]=proj_st(frames,cut_size);
    end;
    ref=tom_cut_out(ref,'center',cut_size);
    act_remove_auto_corr=0;
end;
disp_it('done!',verbose);

%allocate memory 4 all averages and shifts
if (num_out_arg>2)
    all_avg=zeros(cut_size(1),cut_size(2),num_iter);
    all_avg(:,:,1)=ref;
end;
all_shifts=zeros(num_of_frames,2);
all_shifts_old=zeros(num_of_frames,2);
alg_param='';
%loop over all iterations
for itr_count=1:num_iter
    %working loop over all frames
    for i=1:num_of_frames
        [slWindowIdx,slweightIdx]=get_slWindowIdx(i,num_of_frames,sl_window);
        if (iscell(frames))
            tmp_frames=frames(slWindowIdx);
        else
            tmp_frames=frames(:,:,slWindowIdx);
        end;
        [all_shifts(i,:),all_cc(i)]=align2ref(ref,tmp_frames,cut_size,filt(itr_count),mid,idx_ccf,act_remove_auto_corr,all_shifts(i,:),sl_weights(slweightIdx),i,verbose);
    end;
    disp_it(['iter: ' num2str(itr_count) ' mean cc: ' num2str(mean(all_cc)) ' (filt=' num2str(filt(itr_count)) ') mean shift change: ' num2str(mean(abs(all_shifts-all_shifts_old))) ],verbose);
    ref=proj_st(frames,cut_size,all_shifts);
    %tom_emwrite(['db_ref_itr_' num2str(itr_count) '.em'],ref);
    alg_param=cat(2,all_shifts,all_cc');
    if (num_out_arg > 2)
        all_avg(:,:,itr_count+1)=ref;
    end;
    if ((mean(mean(abs((all_shifts-all_shifts_old))))<0.02) && (filt(itr_count) == min(filt))) 
        disp_it('Alignment converged!',verbose);
        break;
    end;
    all_shifts_old=all_shifts; 
    clear('all_cc');
end;

%transfer output variables
if (num_out_arg > 2) 
    all_avg=all_avg(:,:,1:itr_count);
end;
avg=tom_xraycorrect2(single(ref)./num_of_frames);

function [shift,ccc]=align2ref(ref,frames,cut_size,filt,mid,idx_ccf,remove_auto_corr,shifts_old,sl_weights,count,verbose)

img_tmp=zeros(size(ref),'single');
if (iscell(frames))
    for i=1:length(frames)
        img_tmpf=tom_mrcread(frames{i});
        img_tmp=img_tmp+tom_cut_out(img_tmpf.Value,'center',cut_size);
    end;
else
    for i=1:size(frames,3)
        img_tmp=img_tmp+single((frames(:,:,i).*sl_weights(i)));
    end;
end;

if (remove_auto_corr==1)
    ref_tmp=tom_filter(ref-tom_shift(single(img_tmp),shifts_old),filt);
else
    ref_tmp=tom_filter(ref,filt);
end;
img_tmp=tom_filter(img_tmp,filt);
ccf=tom_corr(img_tmp,ref_tmp,'norm');
ccf(idx_ccf)=-1;
[ccc_pos_trans,ccc]=tom_peak(ccf,'spline');
shift=ccc_pos_trans-mid;

%[a b c]=fileparts(frame);
disp_it(['frame nr: ' num2str(count) ' cc: ' num2str(ccc) ' shift: ' num2str(shift(1)) ' ' num2str(shift(2)) ],verbose);


function [sl_Idx,sl_weiIdx]=get_slWindowIdx(pos,num_frames,window_sz)

half_len=((window_sz-1)/2);
if (half_len<1)
    half_len=0;
end;

start=pos-half_len;
if (start<1)
    start=1;
end;
stop=pos+half_len;
if (stop>num_frames)
    stop=num_frames;
end;

sl_Idx=start:stop;

sl_weiIdx=1:window_sz;
if (pos <= half_len)
   sl_weiIdx=sl_weiIdx(median(sl_weiIdx)-(pos-1):end);
end;

function [avg,st]=proj_st(frames,sz_img,shifts)

if (iscell(frames))
    num_of_frames=length(frames);
else
    num_of_frames=size(frames,3);
end;

sh_flag=1;
if (nargin < 3)
    shifts=zeros(length(num_of_frames),2);
    sh_flag=0;
end;

avg=zeros(sz_img,'single');
if (nargout>1)
    st=zeros(sz_img(1),sz_img(2),num_of_frames,'int8');
end;
for i=1:num_of_frames
    if (iscell(frames))
        img_tmp=tom_mrcread(frames{i});
        img_tmp=tom_cut_out(img_tmp.Value,'center',sz_img);
    else
        img_tmp=frames(:,:,i);
    end;
    
    if (sh_flag==1)
        avg=avg+tom_shift(single(img_tmp),shifts(i,:));
    else
        avg=avg+single(img_tmp);
    end;
    if (nargout>1)
        st(:,:,i)=img_tmp;
    end;
end;

function disp_it(message,verbose,log)

if (nargin < 3)
    log='';
end;

if (verbose==1)
    disp(message);
end;

if (isempty(log)==0)
    fid_log=fopen(log,'a');
    fprintf(fid_log,'%s\n',message);
    fclose(fid_log);
end;




