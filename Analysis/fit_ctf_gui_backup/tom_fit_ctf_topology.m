function [topo statistic]=tom_fit_ctf_topology(img,Fit,EM,Search,grid_size,region_size,ctf_size,step)


img_size= size(img)./grid_size;


Search_tmp=Search;

il=0;
iil=0;

for i=1:step:size(img,1)-region_size
        il=il+1
    for ii=1:step:size(img,2)-region_size

        iil=iil+1;
        
region=img(i:i+region_size-1,ii:ii+region_size-1);

ps=tom_calc_periodogram(double(region),ctf_size);
ps=(log(fftshift(ps)));

inner_rad_bg=8;
outer_rad_bg=(ctf_size)/2-1;



[decay decay_image]=calc_decay(ps,inner_rad_bg,outer_rad_bg,16);
background_corrected_ps=double(ps-decay_image);


mask_in_radius=9;
mask_out_radius=64;
ps_size=size(ps);
mask_in = tom_spheremask(ones(ps_size),mask_in_radius,0,[ps_size(1)./2+1 ps_size(2)./2+1 1]);
mask_out = tom_spheremask(ones(ps_size),mask_out_radius,0,[ps_size(1)./2+1 ps_size(2)./2+1 1]);
mask=mask_out-mask_in;



ps=background_corrected_ps.*mask;
tom_imagesc(ps); drawnow;
warning off;
[Fit]=tom_fit_ctf(ps,EM,Search_tmp);
warning on;
Fit.Dz_det
topo(il,iil)=Fit.Dz_det;
[a1 a2 a3 a4 ]=tom_dev(Fit.corr_all);
statistic.cc_std(il,iil)=a4;
statistic.cc_mean(il,iil)=a1;


end;
iil
iil=0;
end;