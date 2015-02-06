 function tom_seg_26S(model_name,pixs,sym_use,outputfolder,use_seg,ini_rot,thresh2,cut_size)
% TOM_SEG_26S segments atpase and 20s of 26S
%  
%     tom_seg_26S(model_name,pixs,outputfolder,use_seg,thresh2,cut_size)
%  
%  PARAMETERS
%  
%    INPUT
%     model_name          name of 26S input struct
%     pixs                pixelsize
%     sym_use             symmetry used 
%     outputfolder        (opt.) outputfoldername (default=out_seg)
%     use_seg             (opt.) vector describing which structures should be
%                         fitted first=20s second=upper atpase third=loweratpase
%                         (default [1 1 1])     
%     ini_rot             (opt.) initial rotation of model (default=[0 0 0]) to have sym-axis along 7-fold axis  
%     thresh2             (opt.) thesh for 26S (default=calib20s) 
%     cut_size            (opt.) size for cutting out sym. struct (default=size(model)./2)
%   
%    
%    OUTPUT
%
%  
%  EXAMPLE
%   tom_seg_26S('11__25_ctf_r1i73_reconstruction.vol',4.42,'C1','output_seg');
%  
%   %segment only 20S and use initial rotation
%   tom_seg_26S('11__25_ctf_r1i73_reconstruction.vol',4.42,'C2','output_seg',[1 0 0],[270 90 90]);
%
%  REFERENCES
%  
%  SEE ALSO
%     ...
%  
%     created by SN/FB 01/24/06
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom

 
try
    im=tom_spiderread(model_name);
catch ME
    im=tom_emread(model_name);
end;



im=im.Value;
sz_im=size(im); 
 
 
 if (nargin<2)
    pixs=4.42;
    disp('Warning no pixelsize given assuming 4.42 Ang');
 end;

if (nargin<4 || isempty(outputfolder)) 
    outputfolder='out_seg'; 
end;

if (nargin<5)
    use_seg=[1 1 1];
end;
  
if (nargin<6 || isempty(ini_rot) ) 
    ini_rot=[0 0 0];
end;

if (nargin<7 || isempty(thresh2) ) 
    thresh2='calib20s'; 
end;


if (nargin<8)
    cut_size=round(sz_im./2);
end;


im=tom_rotate(im,ini_rot);

drosophila=1;

if (drosophila)
    mass_20s=720;
    mass_Atpase=350;
    mass_Atpase_withCC=300;
else
    mass_20s=755;
    mass_Atpase=225;
    mass_Atpase_withCC=280;
end;



use_err_dil=1;

do_20s=use_seg(1);
do_atpase1=use_seg(2);
do_atpase2=use_seg(3);

warning off;
mkdir(outputfolder);
warning on;

warning off;
mkdir([outputfolder '/log' ]);
warning on;


if (sum(use_seg)~=0)
     fid=fopen([outputfolder '/log/log.txt' ],'wt');
end;


if (do_20s==1)
    
     %get 20s axis
    disp('Warning Filter hack before get sym profile !');
    [profile minima]=tom_get_sym_profile(tom_filter(im,1),7,'cent_proj',0,1);
    [sub coord]=get_subvolume(im,minima,[-2 2],cut_size,sym_use);
    
    load_fl=0;
    
    if (load_fl==0)
        [ang_back shift_back sym_struct aligned_sym_struct]=tom_get_sym_axis(sub,7,3,[-8 4 8],[-8 4 8],[-8 4 8],5,[0 0 0],[0 0 0],'default',[outputfolder '/20s/debug_sym_alg']);
        fprintf(fid,'20s: ang: %f %f %f  shift: %f %f %f \n',ang_back(1),ang_back(2),ang_back(3),shift_back(1),shift_back(2),shift_back(3));
        im_al=tom_shift(tom_rotate(im,[-ang_back(2) -ang_back(1) -ang_back(3)]),-shift_back); 
    else
        im_al=tom_emread([outputfolder '/26S_aling.em']);
        im_al=im_al.Value;
    end;
    
    
    tom_emwrite([outputfolder '/26S_aling.em'],im_al);
    
    
    %do it again to get a better z segmentation
    [ang_back,shift_back,thresh20s]=do_segment(im_al,7,cut_size,[-2 2],[outputfolder '/20s'],pixs,mass_20s,mass_20s,thresh2,sym_use,0);
    
else
   im_al=tom_emread([outputfolder '/26S_aling.em']);
   im_al=im_al.Value;
end;

%segment ATPase
if (do_atpase1)
    
    [ang_back,shift_back]=do_segment(im_al,6,cut_size,[-4 -2],[outputfolder '/atpase1'],pixs,mass_Atpase_withCC,mass_Atpase,thresh2,sym_use,1);
   % [euler_7_al shift_out rott]=tom_sum_rotation([euler_y_ax; -ang_back(1) -ang_back(2) -ang_back(3); -180 90 -90 ],[0 0 0; -shift_back' ;0 0 0]);
    
    fprintf(fid,'Atpase1: ang: %f %f %f  shift: %f %f %f \n',ang_back(1),ang_back(2),ang_back(3),shift_back(1),shift_back(2),shift_back(3));
    
    %fprintf(fid,'Atpase1: ang: %f %f %f  shift: %s %s %s \n',euler_7_al(1),euler_7_al(2),euler_7_al(3),shift_out(1),shift_out(2),shift_out(3));
end;

if (do_atpase2)
    
    [ang_back,shift_back]=do_segment(im_al,6,cut_size,[2 4],[outputfolder '/atpase2'],pixs,mass_Atpase_withCC,mass_Atpase,thresh2,sym_use,1);
    %[euler_7_al shift_out rott]=tom_sum_rotation([euler_y_ax; -ang_back(1) -ang_back(2) -ang_back(3); -180 90 -90 ],[0 0 0; -shift_back' ;0 0 0]);
    
    fprintf(fid,'Atpase2: ang: %f %f %f  shift: %f %f %f \n',ang_back(1),ang_back(2),ang_back(3),shift_back(1),shift_back(2),shift_back(3));
    
    %fprintf(fid,'Atpase2: ang: %f %f %f  shift: %s %s %s \n',euler_7_al(1),euler_7_al(2),euler_7_al(3),shift_out(1),shift_out(2),shift_out(3));
   
end;

%calculate rest
warning off;
    mkdir([outputfolder '/all']);
warning on;

a1=tom_emread([outputfolder '/atpase1/mask_half.em']);
a2=tom_emread([outputfolder '/atpase2/mask_half.em']);
p=tom_emread([outputfolder '/20s/mask_full.em']);
%p=tom_emread([outputfolder '/20s/mask.em']);

if (strcmp(thresh2,'calib20s'))
  cyl=tom_cylindermask(ones(128,128,128),12);
  sig=find(squeeze(sum(sum(p.Value,2),1))>0);
  mi=min(sig); ma=max(sig);  
  cyl(:,:,1:mi)=0; cyl(:,:,ma:end)=0;
  [thresh2 mass vol vol_bin]=tom_calc_isosurface(cyl.*im_al,mass_20s,pixs,.01);
  thresh2=1.9;
end;

mask_all=(a1.Value+a2.Value+p.Value);
rest_all=(im_al>thresh2)-mask_all;

r_proc=rest_all;

if (use_err_dil)
    %m_a1=tom_circumscr_mask(a1.Value,'cyl',[1 1]);
    %m_a2=tom_circumscr_mask(a2.Value,'cyl',[1 1]);
    bsz=4;
    se=strel('ball',bsz,bsz);
    m_a1=imdilate(a1.Value,se)>bsz;
    m_a2=imdilate(a2.Value,se)>bsz;
    m_p=tom_circumscr_mask(p.Value,'cyl',[2 2]);
    m3=(m_a1+m_a2+m_p)>0;
    %se=strel('disk',2);
    se=strel('ball',4,4);
    r_erod=imerode(rest_all,se);
    r_dill=imdilate(r_erod,se);
    r_dill=r_dill==1;
    r_proc=(r_dill.*m3+rest_all.*(m3==0))>0;
end;


tom_emwrite([outputfolder '/all/rest.em'],r_proc);
tom_emwrite([outputfolder '/all/rest_org.em'],rest_all);

function [ang_back,shift_back,thresh]=do_segment(vol,sym,cut_size,mins,outpath,pixs,mass_half,mass_full,th,sym_use,app_half_sym)


warning off;
mkdir(outpath);
warning on;

sz_im=size(vol);

disp('Warning Filter hack before get sym profile !');
[profile minima]=tom_get_sym_profile(tom_filter(vol,1),sym,'cent_proj',0,0);
%[profile minima]=tom_get_sym_profile(vol,sym,'cent_proj',20,0); old stuff

[sub coord]=get_subvolume(vol,minima,[mins(1) mins(2)],cut_size,sym_use);


[ang_back shift_back sym_struct aligned_sym_struct]=tom_get_sym_axis(sub,sym,3,[-16 8 16],[-16 8 16],[-16 8 16],5,[0 0 0],[0 0 0],'default',[outpath '/debug_sym']);

%n fold mask
tmp=zeros(sz_im); 
mask_org=tom_paste(tmp,aligned_sym_struct,coord);


%n/2 fold mask
sub_al=tom_shift(tom_rotate(sub,[-ang_back(2) -ang_back(1) -ang_back(3)]),-shift_back); 
sym_struct_3=tom_symref(sub_al,3);
aligned_sym_struct_3=tom_shift(tom_rotate(sym_struct_3,[ang_back(1) ang_back(2) ang_back(3)]),shift_back);

tmp=zeros(sz_im); 
mask3=tom_paste(tmp,aligned_sym_struct_3,coord);
mask3_org=mask3;

[mask6,mask3]=tom_seg_26S_clean_up(mask_org,mask3_org,pixs,2000,80,18,2.3,mass_half,mass_full);

tom_emwrite([outpath '/mask_full.em'],mask6);
tom_emwrite([outpath '/full.em'],mask_org);
tom_emwrite([outpath '/mask_half.em'],mask3);
tom_emwrite([outpath '/half.em'],mask3_org);


thresh=1;


function [sub coord]=get_subvolume(vol,minima,offset,sz,sym_use)

sz_im=size(vol);
mid_im=floor(sz_im./2)+1;

cut_start=round((sz_im-sz)./2);
cut_stop=sz_im-round((sz_im-sz)./2)-1;

if (strcmp(sym_use,'C2'));
    [ind_start_t1 start_val_t1]=tom_nearestpoint(mid_im(1)-1,minima);
    [ind_start_t2 start_val_t2]=tom_nearestpoint(mid_im(1)+1,minima);
    tmp=round((start_val_t1+start_val_t2)./2);
    minima(ind_start_t1)=tmp;
    minima(ind_start_t2)=tmp;
    minima=unique(minima);
end;

[ind_start start_val]=tom_nearestpoint(mid_im(1),minima);  


stop=minima(ind_start-offset(1));
start=minima(ind_start-offset(2));

org=vol(cut_start:cut_stop,cut_start:cut_stop,start:stop);
pos=[1 1 round((sz(1)-size(org,3) )./2) ];
sub=tom_paste(zeros(sz),org,pos);

coord=[cut_start(1) cut_start(2)  (start-round((sz(1)-size(org,3) )./2)+1)];
