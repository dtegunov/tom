function [result quality]=tom_check_image_quality(in)

load 'EM_for_26S_at_Titan.mat'
load 'Search_for_26S_at_Titan.mat'
% power spectrum, log!
ps=tom_calc_periodogram_parallel(single(in.Value),256,0,floor(256./16));
ps=(fftshift(ps));
fit_range=[3.5 7];
bfactor=fit_bfactor(ps,in.Header.Objectpixelsize,fit_range);

quality.bfactor=bfactor;
ps=log(ps);
% correct PS for background
correctbackground_inner_radius=4;
correctbackground_outer_radius=64;
[decay decay_image]=calc_decay(ps,correctbackground_inner_radius,correctbackground_outer_radius,32);
background_corrected_ps=double(ps-decay_image);

% mask PS
img_size=size(background_corrected_ps);
mask_in_radius=5;
mask_out_radius=32;
mask_in = tom_spheremask(ones(img_size),mask_in_radius,0,[img_size(1)./2+1 img_size(2)./2+1 1]);
mask_out = tom_spheremask(ones(img_size),mask_out_radius,0,[img_size(1)./2+1 img_size(2)./2+1 1]);
mask=mask_out-mask_in;
background_corrected_masked_ps=background_corrected_ps.*mask;

% adapt defocus search range, 1 mue up and down
dz=in.Header.Defocus.*(1e-10);
Search.Dz_search=[dz-1e-6:0.1e-6:dz+1e-6];
[Fit]=tom_fit_ctf(background_corrected_masked_ps,EM,Search);
quality.ps=background_corrected_masked_ps;
quality.Fit=Fit;
sz=size(quality.Fit.corr_all);
r=reshape(quality.Fit.corr_all,[sz(1).*sz(2).*sz(3).*sz(4) 1]);
quality.test.ccc_min=min(r);
quality.test.ccc_max=max(r);
quality.test.ccc_mean=mean(r);
quality.test.ccc_std=std(r);
quality.test.confidence=(quality.test.ccc_max-quality.test.ccc_mean)./quality.test.ccc_std;
if quality.test.confidence>1.5 && quality.test.ccc_max>0.1
    quality.result=1;
    result=1;
else
    quality.result=0;
    result=0;
end;


 



function bfactor=fit_bfactor(ps,objectpixelsize,fit_range)

lln=calc_fourier_shell(sqrt(ps));
lln=lln./max(lln);
decay=log(lln)';
x=(2*objectpixelsize.*length(decay))./([1:1:length(decay)]);
x_plot=[1:1:length(decay)];
idx=1;
for i=1:length(decay)
    if x(i)>=fit_range(1) && x(i)<=fit_range(2)
        decay_area(idx)=decay(i);
        x_area(idx)=x(i);
        x_area_plot(idx)=x_plot(i);
        idx=idx+1;
    end;
end;
q_square_area=1./((x_area.^2));
var = polyfit(q_square_area,decay_area,1);
a = var(2); % y-intercept of the fitted line
b = var(1); % slope of fitted line
bfactor=4.*var(1);


function [l nr]=calc_fourier_shell(in)

l=0;
sz=size(in,1);
sz_2=sz./2+1;
idx=1;
[x,y]=ndgrid(0:size(in,1)-1,0:size(in,2)-1);
v = sqrt((x+1-sz_2).^2+(y+1-sz_2).^2);
for r=0:sz_2-1
    ind = find(round(v)==r);
    l(idx)=mean(in(ind));
    nr(idx)=length(ind);
    idx=idx+1;
end
