 function tom_seg_26S(model_name,pixs,outputfolder,use_seg,th2,cut_size)


 
im=tom_spiderread(model_name);
im=im.Value;
sz_im=size(im); 
 
 
 if (nargin<2)
    pixs=4.42;
    disp('Warning no pixelsize given assuming 4.42 Ang');
 end;

if (nargin<3) 
    outputfolder='out_seg'; 
end;

if (nargin<4 || isempty(th2) ) 
    th2='equal'; 
end;

if (nargin<5)
    cut_size=round(sz_im./2);
end;



mass_1ryp=755;
mass_Atpase=300;




warning off;
mkdir(outputfolder);
warning on;



mid_im=floor(sz_im./2)+1;

do_20s=0;

if (do_20s==1)
    
    [profile minima]=tom_get_sym_profile(im,7,'2proj',0,1);
    
    [sub coord]=get_subvolume(im,minima,[-2 2],cut_size);
    
    [euler_y_ax shift_out rott]=tom_sum_rotation([270 90 90; 90 0 0],[0 0 0; 0 0 0]);
    
    sub_rot=tom_rotate(sub,euler_y_ax);
    
    [ang_back shift_back sym_struct aligned_sym_struct]=tom_get_sym_axis(sub_rot,7,3,[-8 4 8],[-8 4 8],[-8 4 8],5);
    
    %align 26s along 7-fold axis
    [euler_7_al shift_out rott]=tom_sum_rotation([euler_y_ax; -ang_back(1) -ang_back(2) -ang_back(3); -180 90 -90 ],[0 0 0; -shift_back' ;0 0 0]);
    im_al=tom_shift(tom_rotate(im,euler_7_al),shift_out);
    
    tom_emwrite([outputfolder '/26S_aling.em'],im_al);
    
    %do it again to get a better z segmentation
    [profile minima]=tom_get_sym_profile(im_al,7,'2proj',0,1);
    
    [sub coord]=get_subvolume(im_al,minima,[-2 2],cut_size);
    
    [euler_y_ax shift_out rott]=tom_sum_rotation([270 90 90; 90 0 0],[0 0 0; 0 0 0]);
    
    sub_rot=tom_rotate(sub,euler_y_ax);
    
    [ang_back shift_back sym_struct aligned_sym_struct]=tom_get_sym_axis(sub_rot,7,3,[-1 0.5 1],[-1 0.5 1],[-1 0.5 1],2);
    
    
    tmp=zeros(sz_im);
    
    mask20s=tom_paste(tmp,tom_rotate(aligned_sym_struct,[-180 90 -90]),coord);
    
    [thresh mass vol mask20s]=tom_calc_isosurface(mask20s,mass_1ryp,pixs,0.001);
    
    %write file output
    
    warning off;
    mkdir([outputfolder '/20s']);
    warning on;
    
    tom_emwrite([outputfolder '/20s/mask20s.em'],mask20s);
    struct=tom_paste(tmp,tom_rotate(sym_struct,[-180 90 -90]),coord);
    tom_emwrite([outputfolder '/20s/sym20s.em'],struct);
    [bfactor decay_restore struct_corr]=tom_fit_bfactor(sym_struct,4.42,[5 20],[10 Inf],2500000,1,0);
    tom_emwrite([outputfolder '/20s/sym20s_corr.em'],tom_paste(tmp,tom_rotate(struct_corr,[-180 90 -90]),coord));
    [bfactor decay_restore struct_corr]=tom_fit_bfactor(aligned_sym_struct,4.42,[5 20],[10 Inf],2500000,1,0);
    tom_emwrite([outputfolder '/20s/sym20s_alg_corr.em'],tom_paste(tmp,tom_rotate(struct_corr,[-180 90 -90]),coord));
    
    rest=(im_al>th2)-mask20s;
    tom_emwrite([outputfolder '/20s/rest.em'],rest);
    
else
    im_al=tom_emread('mod1/26S_aling.em');
    im_al=im_al.Value;
end;

%segment ATPase

[profile minima]=tom_get_sym_profile(im_al,6,'2proj',0,1);

[sub coord]=get_subvolume(im_al,minima,[2 4],cut_size);

[euler_y_ax shift_out rott]=tom_sum_rotation([270 90 90; 90 0 0],[0 0 0; 0 0 0]);

sub_rot=tom_rotate(sub,euler_y_ax);

[ang_back shift_back sym_struct aligned_sym_struct]=tom_get_sym_axis(sub_rot,6,3,[-8 4 8],[-8 4 8],[-8 4 8],5);


tmp=zeros(sz_im); 

maskAtase1=tom_paste(tmp,tom_rotate(aligned_sym_struct,[-180 90 -90]),coord);

[thresh mass vol maskAtase1]=tom_calc_isosurface(maskAtase1,mass_Atpase,pixs,0.001);

warning off;
mkdir([outputfolder '/atpa1']);
warning on;

tom_emwrite([outputfolder '/atpa1/maskAtpa1.em'],maskAtase1);
struct=tom_paste(tmp,tom_rotate(sym_struct,[-180 90 -90]),coord);
tom_emwrite([outputfolder '/atpa1/atpa1.em'],struct);
[bfactor decay_restore struct_corr]=tom_fit_bfactor(sym_struct,4.42,[5 20],[10 Inf],2500000,1,0);
tom_emwrite([outputfolder '/atpa1/atpa1_corr.em'],tom_paste(tmp,tom_rotate(struct_corr,[-180 90 -90]),coord));
[bfactor decay_restore struct_corr]=tom_fit_bfactor(aligned_sym_struct,4.42,[5 20],[10 Inf],2500000,1,0);
tom_emwrite([outputfolder '/atpa1/atpa1_alg_corr.em'],tom_paste(tmp,tom_rotate(struct_corr,[-180 90 -90]),coord));



rest=(im_al>th2)-maskAtase1;
tom_emwrite([outputfolder '/atpa1/rest.em'],rest);

%segment ATPase2

[profile minima]=tom_get_sym_profile(im_al,6,'2proj',0,1);

[sub coord]=get_subvolume(im_al,minima,[-4 -2],cut_size);

[euler_y_ax shift_out rott]=tom_sum_rotation([270 90 90; 90 0 0],[0 0 0; 0 0 0]);

sub_rot=tom_rotate(sub,euler_y_ax);

[ang_back shift_back sym_struct aligned_sym_struct]=tom_get_sym_axis(sub_rot,6,3,[-8 4 8],[-8 4 8],[-8 4 8],5);


tmp=zeros(sz_im); 

maskAtase1=tom_paste(tmp,tom_rotate(aligned_sym_struct,[-180 90 -90]),coord);

[thresh mass vol maskAtase1]=tom_calc_isosurface(maskAtase1,mass_Atpase,pixs,0.001);

warning off;
mkdir([outputfolder '/atpa2']);
warning on;

tom_emwrite([outputfolder '/atpa2/maskAtpa2.em'],maskAtase1);
struct=tom_paste(tmp,tom_rotate(sym_struct,[-180 90 -90]),coord);
tom_emwrite([outputfolder '/atpa2/atpa2.em'],struct);
[bfactor decay_restore struct_corr]=tom_fit_bfactor(sym_struct,4.42,[5 20],[10 Inf],2500000,1,0);
tom_emwrite([outputfolder '/atpa2/atpa2_corr.em'],tom_paste(tmp,tom_rotate(struct_corr,[-180 90 -90]),coord));
[bfactor decay_restore struct_corr]=tom_fit_bfactor(aligned_sym_struct,4.42,[5 20],[10 Inf],2500000,1,0);
tom_emwrite([outputfolder '/atpa2/atpa2_alg_corr.em'],tom_paste(tmp,tom_rotate(struct_corr,[-180 90 -90]),coord));

th2=0.0186;

rest=(im_al>th2)-maskAtase1;
tom_emwrite([outputfolder '/atpa2/rest.em'],rest);








function [sub coord]=get_subvolume(vol,minima,offset,sz)

sz_im=size(vol);
mid_im=floor(sz_im./2)+1;

cut_start=round((sz_im-sz)./2);
cut_stop=sz_im-round((sz_im-sz)./2)-1;

[ind_start start_val]=tom_nearestpoint(mid_im(1),minima);  


stop=minima(ind_start-offset(1));
start=minima(ind_start-offset(2));

org=vol(cut_start:cut_stop,cut_start:cut_stop,start:stop);
pos=[1 1 round((sz(1)-size(org,3) )./2) ];
sub=tom_paste(zeros(64,64,64),org,pos);

coord=[cut_start(1) cut_start(2)  (start-round((sz(1)-size(org,3) )./2)+1)];
