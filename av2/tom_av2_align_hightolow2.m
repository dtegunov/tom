function tom_av2_align_hightolow2(inalignfile,post_fix,outalignfile,corr_flag,max_translation,Filter_kernel,binning,xray_corr,verboseflag)
%TOM_AV2_ALIGN_HIGHTOLOW2 aligns focal pairs
%
%   tom_av2_align_hightolow(inalignfile,post_fix,outalignfile,corr_flag,max_translation,Filter_kernel,binning,xray_corr,verboseflag)
%
%  TOM_AV2_ALIGN_HIGHTOLOW2 transforms the picklist (align2d mat-file) from
%  the high defocus images into a picklist for the low defocus images
%  
%
%PARAMETERS
%
%  INPUT
%   inalignfile         filename of the input alignment file
%   post_fix            extension for high and low folder for e.g. 'corr'
%                       use '' for high and low    
%   outalignfile        filename of the output alignment file
%   corr_flag           ('subregions') or 'micrograph'
%   max_translation     max translation between the focal pairs 
%   Filter_kernel       size of the circular filter cernel (0==off)
%   binning             binning of the images
%   xray_corr           flag for x-ray correction (0==off/1==on)    
%   verboseflag         1 for output messages, 0 for quiet operation 
%  OUTPUT
%
%EXAMPLE
%    
%   % use feature tracking for registration
%   tom_av2_align_hightolow2('log/pick/high/pick.mat','','log/pick/low/pick.mat','subregions',300,0,2,1,1);
%   creates the picklist log/pick/low/pick.mat 
%
%   % use correlation of complete images ...classic function
%   tom_av2_align_hightolow2('log/pick/high/pick.mat','corr','log/pick/low/pick.mat','micrograph',300,1,2,1,1);
%   creates the picklist log/pick/low/pick.mat 
%   
%   NOTE: 1: All Parameters refer to the original images Size !!!!!
%         2: The folder for the low images must be at the same level as the
%            high folder and named low !!
%         3: The path in the align2d-file (align2d.filename) must be correct !
%
%
%REFERENCES
%
%SEE ALSO
%   tom_aling2d
%
%   created by AK 1/19/06 pimped by fb
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




try
    s = load(inalignfile);
catch
    error('Input alignment file could not be loaded.');
end


[aa bb]=fileparts(outalignfile);
log_name=[aa '/' bb '.txt'];
fid=fopen(log_name,'w');
inalign = s.align2d;


outalign = inalign;
h = tom_reademheader(inalign(1).filename);
sz = h.Header.Size(1);


if (nargin < 3)
    corr_flag='subregions';
end;

if (nargin < 4)
    max_translation=round(sz./2);
end;

if (nargin < 5)
    Filter_kernel=3;
end;

if (nargin < 6)
    binning=1;
end;

if (nargin < 7 )
    xray_corr=1;
end;

if (nargin < 8)
    verboseflag=1;
end;


if (strcmp(corr_flag,'subregions')==0 && strcmp(corr_flag,'micrograph')==0)
    error(['corr_flag: ' corr_flag ' not implemented ! use subregions or micrograph']);
    return;
end;


if (strcmp(corr_flag,'subregions'))
    fprintf(fid,'%s %s %s %s %s %s\n','Filename','shift_x','shift_y','std of shift','num of feature','mean cc');
else
    fprintf(fid,'%s %s %s %s\n','Filename','shift_x','shift_y','cc');
end;
base_num_of_feat=40;
std_thresh=0.2; % in percent of mean
start_method='feature';
subregion_size=256./(2^binning);
sub_m=start_method;
num_of_feat=base_num_of_feat;
min_num_of_backm_feat=5;
min_mean_cc=0.1;
auto_corr_rad=4;
std_clean=2;
demo=0;
num_of_try=2;


mask_trans=tom_bin(tom_spheremask(ones(sz,sz),max_translation,round(sz./20)),binning);


for i=1:20
    try
        h_tmp=tom_reademheader(inalign(i).filename);
        h_tmp=h_tmp.Header;
        break;
    catch ME
    end;
end;
h_tmp.Size=[h_tmp.Size(1).*2^binning h_tmp.Size(2).*2^binning 0];

for i=1:size(inalign,2)
    all_filenames_h{i}=inalign(i).filename;
end;
all_f_high=unique(all_filenames_h);


for i=1:size(inalign,2)
    all_filenames_l{i}=regexprep(inalign(i).filename,['/high' post_fix '/'], ['/low' post_fix '/']);
    [a b c]=fileparts(all_filenames_l{i});
    if (isempty(findstr(b,'high_'))==0)
        all_filenames_l{i}=[a '/' strrep(b,'high','low') c];
    end;
end;
all_f_low=unique(all_filenames_l);


if (strcmp(all_f_low{1},all_f_high{1})==1)
    disp('file organisation not tom style or wrong post_fix !!');
    disp('use high and low or high_corr and low_corr foldernames only');
    disp('use post_fix for high_corr and low_corr');
    error(' ');
    return;
end;

auto_cc_mask=tom_bin(tom_spheremask(ones(sz,sz),8,1),binning)==0;

parfor i=1:length(all_f_high)
    corr_images(log_name,all_f_high{i},all_f_low{i},binning,corr_flag,xray_corr,sub_m,base_num_of_feat,std_thresh,num_of_try,subregion_size,num_of_feat,start_method,min_num_of_backm_feat,min_mean_cc,std_clean,auto_corr_rad,auto_cc_mask,max_translation,mask_trans,Filter_kernel,verboseflag,demo);
end;

log=importdata(log_name);


for i=1:size(inalign,2)
    outalign(i).filename = all_filenames_l{i};
    idx=find(ismember(log.textdata,{inalign(i).filename}))-1;
    shift_out=[log.data(idx,1) log.data(idx,2)];
    pos_tmp=[inalign(i).position.x inalign(i).position.y];
    pos_tmp=pos_tmp-shift_out(1:2);
    outalign(i).position.x = round(pos_tmp(1));
    outalign(i).position.y = round(pos_tmp(2));
end;

if (strcmp(corr_flag,'subregions'))
    for i=1:size(log.data,1)
        [a b c]=fileparts(log.textdata{i+1});
        particlepicker.filelist{i}=[b c];
        shift_out=[log.data(idx,1) log.data(idx,2)];
        shifts_std=log.data(i,3);
        hit=log.data(i,4);
        cc_mea=log.data(i,5);
        if ( max(shifts_std) < min(abs(std_thresh.*shift_out)) && hit > min_num_of_backm_feat && cc_mea > min_mean_cc)
            particlepicker.filefilter{i}=0;
        else
            particlepicker.filefilter{i}=1;
        end;
    end;
    particlepicker.maskcell={};
    [a idx_sort]=sort(particlepicker.filelist);
    particlepicker.filefilter=particlepicker.filefilter(idx_sort);
    particlepicker.filelist=particlepicker.filelist(idx_sort);
end;



align2d = outalign;
save(outalignfile,'align2d');
if (strcmp(corr_flag,'subregions'))
    save([aa bb '_filefilter.mat'],'particlepicker');
end;
fclose(fid);


function corr_images(log_name,f_name_h,f_name_l,binning,corr_flag,xray_corr,sub_m,base_num_of_feat,std_thresh,num_of_try,subregion_size,num_of_feat,start_method,min_num_of_backm_feat,min_mean_cc,std_clean,auto_corr_rad,auto_cc_mask,max_translation,mask_trans,Filter_kernel,verbose,demo)

try
    im_h = tom_emreadc(f_name_h,'binning',[binning binning 0]);
    im_h=im_h.Value;
catch
    warning(['File ' f_name_h ' not found']);
    return;
end

if  (Filter_kernel > 0)
    im_h=tom_filter(im_h,Filter_kernel);
end;

try
    im_l = tom_emreadc(f_name_l,'binning',[binning binning 0]);
    h_tmp=im_l.Header;
    im_l=im_l.Value;
catch
    warning(['File ' f_name_l ' not found']);
    return;
end;
if  (Filter_kernel > 0)
     im_l=tom_filter(im_l,Filter_kernel);
end;

h_tmp.Size=[h_tmp.Size(1).*2^binning h_tmp.Size(2).*2^binning 1];

if (xray_corr==1)
    im_h = tom_xraycorrect2(im_h);
    im_h=tom_smooth(im_h,30);
end;
if (xray_corr==1)
    im_l = tom_xraycorrect2(im_l);
    im_l=tom_smooth(im_l,30);
end;
im_sz=size(im_h);
middle_im=floor(im_sz./2)+1;

if (verbose==1)
tic;
end;

try
    if (strcmp(corr_flag,'subregions'))
        for iii=1:num_of_try
            [shifts shifts_std statistic]=tom_corr_subregions(im_l,im_h,sub_m,subregion_size,num_of_feat,0,0,'xcf',1,max_translation,auto_corr_rad,std_clean,demo,0);
            if (isempty(shifts)==0 && (max(shifts_std) < min(abs(std_thresh.*shifts)) || mean(abs(shifts)) < 5 )  && statistic.hit > min_num_of_backm_feat && mean(statistic.ccc_hit) > min_mean_cc)
                num_of_feat=base_num_of_feat;
                sub_m=start_method;
                break;
            else
                
                if (isempty(shifts)==0)
                    if (max(shifts_std) >= min(abs(std_thresh.*shifts)) || mean(abs(shifts)) < 5)
                        disp([f_name_h ' shift std: ' num2str(max(shifts_std)) ' > ' num2str(min(abs(std_thresh.*shifts)))  ]);
                    end;
                    if (statistic.hit <= min_num_of_backm_feat)
                        disp([f_name_h ' number of remaining features too small (' num2str(statistic.hit) ') ' ]);
                    end;
                    if (isempty(shifts))
                        disp([f_name_h ' no features found']);
                    end;
                    if (mean(statistic.ccc_hit) <= min_mean_cc)
                        disp([f_name_h ' mean ccc too small ' num2str(mean( statistic.ccc_hit)) ]);
                    end;
                else
                     disp([f_name_h ' no shifts after cleaning!']);
                end;
                
                num_of_feat=base_num_of_feat.*(iii.*2);
                if (mod(iii,2)==1)
                    sub_m='grid';
                    %sub_m='feature';
                else
                    sub_m='feature';
                end;
                
                if (iii < num_of_try)
                    if (strcmp(sub_m,'grid'))
                        disp(['Retrying with: ' sub_m]);
                    else
                        disp(['Retrying with: ' sub_m ' and ' num2str(num_of_feat) ' features' ]);
                    end;
                        
                end;
                if (isempty(shifts))
                    shifts=[0 0];
                    shifts_std=[10000 10000];
                end;
            end;
        end;
        shift_out=round(shifts');
        ccc=shifts_std;
        outtext='std ';
    else
        ccf_trans=tom_corr(im_l,im_h,'norm');
        ccf_trans=ccf_trans.*mask_trans.*auto_cc_mask;
        [ccc_pos_trans ccc]=tom_peak(ccf_trans,'spline');
        shift_out=ccc_pos_trans-middle_im;
        shift_out=round(shift_out)';
        outtext='ccc';
    end;
    
catch ME
    shifts=[0 0];
    shifts_std=[10000 10000];
    shift_out=round(shifts');
    ccc(1)=0;
    statistic.ccc_hit=[0 0];
    outtext='err_ccc';
    disp('Error calculating shifts!');
    disp(ME.message);
end;


%mange shifts
fid=fopen(log_name,'a');

if (strcmp(corr_flag,'subregions'))
   fprintf(fid,'%s  %f %f %f %f %f\n',f_name_h,shift_out(1).*2^binning,shift_out(2).*2^binning,ccc(1),statistic.hit,mean(statistic.ccc_hit));
else
   fprintf(fid,'%s  %f %f %f\n',f_name_h,shift_out(1).*2^binning,shift_out(2).*2^binning,ccc(1));
end;

fclose(fid);
h_tmp.Comment=[num2str(shift_out(1).*2^binning) '  ' num2str(shift_out(2).*2^binning)];
tom_writeemheader(f_name_l,h_tmp);

if (verbose==1)
    if (strcmp(corr_flag,'subregions'))
        disp([f_name_h ' shift: ' num2str(shift_out(1).*2^binning) ' ' num2str(shift_out(2).*2^binning) '  ' outtext ' ' num2str(ccc(1)) ' num of feat: ' num2str(statistic.hit) '  mean ccc: ' num2str(mean(statistic.ccc_hit)) ]);
    else
        disp([f_name_h ' shift: ' num2str(shift_out(1).*2^binning) ' ' num2str(shift_out(2).*2^binning) '  ' outtext ' ' num2str(ccc(1)) ]);
    end;
    toc;
end;












%old code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% for i=1:size(inalign,2)
%     outalign(i).filename = regexprep(inalign(i).filename,['/high' post_fix '/'], ['/low' post_fix '/']);
% 
%     %only align files if particle is on a new file
%     if ~isequal(filename_cache,inalign(i).filename)
%         try
%             im_h = tom_emreadc(inalign(i).filename,'binning',[binning binning 0]);
%         catch
%             warning(['File ' inalign(i).filename ' not found']);
%             continue;
%         end
%         if (xray_corr==1)
%             im_h.Value = tom_xraycorrect2(im_h.Value);
%             im_h.Value=tom_smooth(im_h.Value,30);
%         end;
%         try
%             im_l = tom_emreadc(outalign(i).filename,'binning',[binning binning 0]);
%         catch
%             warning(['File ' outalign(i).filename ' not found']);
%             continue;
%         end;
%         
%         if (xray_corr==1)
%             im_l.Value = tom_xraycorrect2(im_l.Value);
%             im_l.Value=tom_smooth(im_l.Value,30);
%         end;
%         im_sz=size(im_h.Value);
%         middle_im=floor(im_sz./2)+1;
%         
%         im_h =im_h.Value;
%         im_l =im_l.Value;
%         angle_out=0;
%         
%         if (strcmp(corr_flag,'subregions'))
%             for iii=1:num_of_try
%                 [shifts shifts_std statistic]=tom_corr_subregions(im_l,im_h,sub_m,subregion_size,num_of_feat,0,0,'xcf',1,max_translation,auto_corr_rad,std_clean,demo,0);
%                 if (isempty(shifts)==0 && (max(shifts_std) < min(abs(std_thresh.*shifts)) || mean(abs(shifts)) < 5 )  && statistic.hit > min_num_of_backm_feat && mean(statistic.ccc_hit) > min_mean_cc)
%                    num_of_feat=base_num_of_feat;
%                    sub_m=start_method;
%                    break;
%                 else
%                     if (max(shifts_std) >= min(abs(std_thresh.*shifts)) || mean(abs(shifts)) < 5)
%                         disp([inalign(i).filename ' shift std: ' num2str(max(shifts_std)) ' > ' num2str(min(abs(std_thresh.*shifts)))  ]);
%                     end;
%                     if (statistic.hit <= min_num_of_backm_feat)
%                         disp([inalign(i).filename ' number of remaining features too small (' num2str(statistic.hit) ') ' ]);
%                     end;
%                     if (isempty(shifts))
%                         disp([inalign(i).filename ' no features found']);
%                     end;
%                     if (mean(statistic.ccc_hit) <= min_mean_cc)
%                         disp([inalign(i).filename ' mean ccc too small ' num2str(mean( statistic.ccc_hit)) ]);
%                     end;
%                     
%                     num_of_feat=base_num_of_feat.*(iii.*2);
%                     if (mod(iii,2)==1)
%                         sub_m='random';
%                     else
%                         sub_m='feature';
%                     end;
%                     
%                     if (iii < num_of_try)
%                         disp(['Retrying with: ' sub_m ' and ' num2str(num_of_feat) ' features' ]);
%                     end;
%                     if (isempty(shifts))
%                         shifts=[0 0];
%                         shifts_std=[10000 10000];
%                     end;
%                 end;
%             end;
%             shift_out=round(shifts');
%             ccc=shifts_std;
%             outtext='std ';    
%             [a b c]=fileparts(outalign(i).filename);
%             img_count=img_count+1;
%             particlepicker.filelist{img_count}=[b c];
%             if (isempty(shifts)==0 && max(shifts_std) < min(abs(std_thresh.*shifts)) && statistic.hit > min_num_of_backm_feat)
%                 particlepicker.filefilter{img_count}=1;
%             else
%                 particlepicker.filefilter{img_count}=0;
%             end;
%         
%         else
%             ccf_trans=tom_corr(im_l,im_h,'norm');
%             ccf_trans=ccf_trans.*mask_trans.*auto_cc_mask;
%             [ccc_pos_trans ccc]=tom_peak(ccf_trans,'spline');
%             shift_out=ccc_pos_trans-middle_im;
%             shift_out=round(shift_out)';
%             outtext='ccc';
%         end;
%         
%         %mange shifts
%         fprintf(fid,'%s  %f %f %f %f\n',inalign(i).filename,shift_out(1).*2^binning,shift_out(2).*2^binning,ccc(1),statistic.hit);
%         disp([inalign(i).filename ' shift: ' num2str(shift_out(1).*2^binning) ' ' num2str(shift_out(2).*2^binning) '  ' outtext ' ' num2str(ccc(1)) ' num of feat: ' num2str(statistic.hit) '  mean ccc: ' num2str(mean(statistic.ccc_hit)) ]);
%         h_tmp.Comment=[num2str(shift_out(1).*2^binning) '  ' num2str(shift_out(2).*2^binning)];
%         tom_writeemheader(outalign(i).filename,h_tmp);
%         
%         filename_cache = inalign(i).filename;
%         
%     end
% 
%     pos_tmp=[inalign(i).position.x inalign(i).position.y];
%     pos_tmp=pos_tmp-shift_out(1:2)'.*2^binning;
%     pos_tmp = tom_pointrotate2d(pos_tmp,-angle_out(1),[center center]);
% 
%     outalign(i).position.x = round(pos_tmp(1));
%     outalign(i).position.y = round(pos_tmp(2));
% 
%     if verboseflag == 1 && mod(i,10) == 0
%         disp([num2str(i) ' of ' num2str(size(inalign,2)) ' particles done']);
%     end
% 
% end
% 
% 
% align2d = outalign;
% save(outalignfile,'align2d');
% if (strcmp(corr_flag,'subregions'))
%     save([aa bb '_filefilter.mat'],'particlepicker');
% end;
% 
% fclose(fid);
% 
% if verboseflag == 1
%     disp('Finished');
% end