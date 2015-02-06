function tom_av2_align_hightolow_beate(assing_high_low,outpath,out_ext,max_translation,max_rotation,Filter_kernel,binning,xray_corr,verboseflag)
%TOM_AV2_ALIGN_HIGHTOLOW2 aligns focal pairs
%
%   tom_av2_align_hightolow_beate(inalignfile,outalignfile,max_translation,max_ro tation,Filter_kernel,binning,xray_corr,verboseflag)
%
%  TOM_AV2_ALIGN_HIGHTOLOW transforms the picklist (align2d mat-file) from
%  the high defocus images into a picklist for the low defocus images
%  
%
%PARAMETERS
%
%  INPUT
%   assing_high_low     text-file with high 2 low assignment (use tom_av2_align_hightolow_list_beate)
%   outpath             base path for aligned high images
%   outext              extension for high images
%   max_translation     max translation between the focal pairs
%   max_rotation        max rotation between the focal pairs    
%   Filter_kernel       size of the circular filter cernel (0==off)
%   binning             binning of the images
%   xray_corr           flag for x-ray correction (0==off/1==on)    
%   verboseflag         1 for output messages, 0 for quiet operation 
%  OUTPUT
%
%EXAMPLE
%     
% 
% for paralell work open  matlab-pool
%
% matlabpool open local 6  
%
% tom_av2_align_hightolow_beate('list.txt','new_high/testal_','.em',550,0,0,2,1,1);
%
% 
%  
%
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_align_hightolow_list_beate
%
%   created  by fb
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



disp('Starting Alignment ...');
disp(' ');
disp(['max_translation: ' num2str(max_translation)]);
disp(['Filter_kernel: ' num2str(Filter_kernel)]);
disp(['binning: ' num2str(binning)]);
disp(['xray_corr: ' num2str(xray_corr)]);
disp(['xray_corr: ' num2str(xray_corr)]);
disp(['output: ' outpath '*' out_ext]);
disp(' ');


disp(['reading ' assing_high_low]);

try
    s = importdata(assing_high_low,' ');
catch
    error('High 2 low assignment file could not be read!!');
end


[aa bb]=fileparts(outpath);

if (exist(aa,'dir')==0)
    mkdir(aa);
end;


fid=fopen([aa bb '.txt'],'w');
fprintf(fid,'%s %s %s %s\n','Filename','shift_x','shift_y','cc');
fclose(fid);

h = tom_reademheader(strtrim(deblank(strtok(s{1},' '))));
sz = h.Header.Size(1);

if (nargin < 4)
    max_translation=round(sz./2);
end;

if (nargin < 5)
    max_rotation=360;
end;

if (nargin < 6)
    Filter_kernel=2;
end;

if (nargin < 7 )
     binning=1;
end;

if (nargin < 8 )
     xray_corr=1;
end;

if (nargin < 9 )
    verboseflag=1;
end;


ang_pix_polar=(max_rotation.*sz)./360;


disp('building masks ... ' );

mask_im =tom_bin(tom_spheremask(ones(sz,sz),(sz./2-round(sz./10)),round(sz./20)),binning);
mask_rot=tom_bin(tom_rectanglemask(ones(round(sz./2),2.*sz),[round(sz./2) ang_pix_polar]),binning);
mask_trans=tom_bin(tom_spheremask(ones(sz,sz),max_translation,round(sz./20)),binning);
auto_cc_mask=tom_bin(tom_spheremask(ones(sz,sz),8,1),binning)==0;

disp(' ');

tic;

parfor i=1:length(s)
    
    [filename_high filename_low]=strtok(s{i},' ');
    filename_high=strtrim(deblank(filename_high));
    filename_low=strtrim(deblank(filename_low));
    
    [a tmp]=fileparts(filename_high);
    [a tmp]=strtok(tmp,'_');
    num=str2double(strrep(tmp,'_',''));
    
    try
        im_h_org = tom_emreadc(filename_high);
        im_h_org.Value = double(im_h_org.Value);
        im_h.Value = tom_bin(im_h_org.Value,binning);
    catch
        disp(['File ' filename_high ' not found']);
        continue;
    end
    
    if (xray_corr==1)
        im_h.Value = tom_xraycorrect2(im_h.Value);
        im_h.Value=tom_smooth(im_h.Value,30);
    end;
    
    if (Filter_kernel>0)
        im_h.Value=tom_filter(im_h.Value,Filter_kernel);
    end;
    
    
    try
        im_l = tom_emreadc(filename_low,'binning',[binning binning 0]);
    catch
        disp(['File ' filename_low ' not found']);
        continue;
    end;
    
    if (xray_corr==1)
        im_l.Value = tom_xraycorrect2(im_l.Value);
        im_l.Value=tom_smooth(im_l.Value,30);
    end;
    
    if (Filter_kernel>0)
        im_l.Value=tom_filter(im_l.Value,Filter_kernel);
    end;
    
    
    im_sz=size(im_h.Value);
    middle_im=floor(im_sz./2)+1;
    
    im_h = tom_norm(double(im_h.Value),'mean0+1std');
    im_l = tom_norm(double(im_l.Value),'mean0+1std');

    ccf_trans=tom_corr(im_l,im_h,'norm');
    ccf_trans=ccf_trans.*mask_trans.*auto_cc_mask;
    [ccc_pos_trans ccc]=tom_peak(ccf_trans,'spline');
    shift_out=ccc_pos_trans-middle_im;
    shift_out=shift_out';
    if verboseflag == 1
        fprintf('%s <==> %s shift: %f %f cc: %f ==> %s \n',filename_high,filename_low,shift_out(1).*2^binning,shift_out(2).*2^binning,ccc(1),[outpath num2str(num) out_ext]);
    end;
    
    
    fid=fopen([aa bb '.txt'],'w');
    fprintf(fid,'%s %f %f %f\n',filename_high,shift_out(1).*2^binning,shift_out(2).*2^binning,ccc(1));
    fclose(fid);
    
    im_h_org.Value=tom_shift(im_h_org.Value,-shift_out);
    
    tom_emwritec([outpath num2str(num) out_ext],im_h_org);
    
    if verboseflag == 1 && mod(i,50) == 0
        disp([num2str(i) ' of ' num2str(length(s)) ' micrographs done']);
        toc;
    end; 
    
end;


if verboseflag == 1
    disp('Finished');
end