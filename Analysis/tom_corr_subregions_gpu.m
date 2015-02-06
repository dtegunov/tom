function [shifts shifts_std statistic]=tom_corr_subregions(img_1,img_2,area_type,area_size,area_nr,area_seed,area_seed_region,method,search_radius,nr_std,demo);

%TOM_CORR_SUBREGIONS calculates a shift between two images based on
%subregion-correlations.
%
%   [shifts shifts_std statistic]=tom_corr_subregions(img_1,img_2,area_type,area_size,area_nr,area_seed,area_seed_region,method,search_radius,nr_std,demo);
%
%
%PARAMETERS
%
%  INPUT
%   img_1, img_2        shifted input images.
%   area_type           defines different type of region definitions:
%                       'random', 'grid', 'feature', 'tiltaxis'
%                       'random': create random postions, spread over
%                       the image.
%                       'grid': create regular grid postions evenly spread over
%                       the image.
%                       'feature': apply an edge-filter then a bandpass filter
%                       and extract subregions.
%                       'tiltaxis': extract subregions from an area along
%                       the tiltaxis. Use the tiltaxis definition by
%                       tom_setmark. Positive y-direction is 0 degree with positive
%                       angles in right-turn direction. Use degrees.
%   area_size           x-y-dimension of subregions.
%   area_nr             number of subregions.
%   area_seed           x-y coordinates of seed positions used for
%                       subregion-correlation. If area_type is set to
%                       'tiltaxis' use this parameter as an tiltangle in
%                       degrees.
%   area_seed_region    area [x1 y1; x2 y2] defining the range of randomly
%                       generated subregions. If area_type is set to
%                       'tiltaxis' use this parameter as an area definition
%                       for the width of the subregion. 0<a<=.5. (0.3 is 30%
%                       of the image size).
%   method              only 'xcf' supported.
%   search_radius       not supported, yet.
%   nr_std              times standard deviation, used for cleaning the
%                       statistics (typically 1, 2, 3).
%   demo                demo mode on=1/off=0.
%                     
%  
%  OUTPUT
%   shifts              shift vector of images.
%   shifts_std          standard deviation of shift vector of images.
%   statistic           complete statistic of subregion-correlation.
%
%EXAMPLE
%   Calculate the shift of image 'in1b' and 'in2b' with randomly spread
%   subregions of size 64x64 pixels. Use 200 subregions, 1std for cleaning
%   and turn the demo mode on;
%   [shifts shifts_std statistic]=tom_corr_subregions(in1b,in2b,'random',64,200,0,0,'xcf',512,1,1);
%
%   Calculate the shift of image 'in1b' and 'in2b' with randomly spread
%   subregions along the tiltaxis with 76 degrees of size 64x64 pixels. The
%   width of the area is 20% of the image.
%   Use 50 subregions, 2std for cleaning and turn the demo mode on;
%   [shifts shifts_std statistic]=tom_corr_subregions(in1b,in2b,'tiltaxis',64,50,76,.2,'xcf',512,2,1);
%
%REFERENCES
%
%SEE ALSO
%   TOM_CORR, TOM_CCC
%
%   created by SN 30/09/08
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


if ~exist('nr_std', 'var') 
    nr_std=3;
end;
if ~exist('area_type', 'var') 
    area_type='random';
end;
if ~exist('demo', 'var') 
    demo=0;
end;
if ~exist('search_radius', 'var') 
    search_radius=size(img_1,1)/2;
end;
if ~exist('method', 'var') 
    method='xcf';
end;
if ~exist('area_size', 'var') 
    area_size = round(size(img_1,1)./8);    
end;
if ~exist('area_nr', 'var') 
    area_nr = 20;
end;

area_size_2 = area_size/2;

idx=1;
img_1=gforce(gsingle(img_1));
img_2=gforce(gsingle(img_2));

switch lower(area_type)
    case 'feature'
        disp('area_seed is feature defined.');        
        area_seed=extract_area_seed_by_features(img_2,area_nr,area_size_2);
    case 'tiltaxis'
        tiltaxis=area_seed-90;
        disp(['area_seed is randomly defined along tiltaxis: ' num2str(area_seed) '.']);        
        area_seed=generate_area_seed_along_tiltaxis(img_2,area_size_2,area_nr,tiltaxis,area_seed_region);
        area_seed_region=0;
    case 'user'
        disp('area_seed is user defined.');
    case 'grid'
        disp('area_seed is defined as a regular grid.');
        clear area_seed;
        d=[area_size./2:area_size:size(img_1,1)-area_size./2];
        for i=1:size(d,2)
            for ii=1:size(d,2)
                area_seed(idx,1)=d(i);
                area_seed(idx,2)=d(ii);
                idx=idx+1;
            end;
        end;
        area_nr=idx-1;
        disp([num2str(area_nr) ' subregions are correlated.']);
    otherwise
        disp('area_seed is randomly defined.');
        if ~exist('area_seed', 'var')
            rand('twister',sum(100*clock));
            area_seed = round(rand(area_nr,2).*(size(img_1,1)-area_size-1));
            area_seed = area_seed + area_size./2;
        elseif area_seed==0
            rand('twister',sum(100*clock));
            area_seed = round(rand(area_nr,2).*(size(img_1,1)-area_size-1));
            area_seed = area_seed + area_size./2;
        end;
end;

if exist('area_seed_region', 'var')
if area_seed_region~=0
    rand('twister',sum(100*clock));
    area_seed_region(1,1)=area_seed_region(1,1)+area_size;
    area_seed_region(1,2)=area_seed_region(1,2)-area_size;
    area_seed_region(2,1)=area_seed_region(2,1)+area_size;
    area_seed_region(2,2)=area_seed_region(2,2)-area_size;
    area_seed_1 =  area_seed_region(1,1) + (area_seed_region(1,2)-area_seed_region(1,1)).*rand(area_nr,1);
    area_seed_2 =  area_seed_region(2,1) + (area_seed_region(2,2)-area_seed_region(2,1)).*rand(area_nr,1);
    area_seed=[];
    area_seed(:,1) = round(area_seed_1);
    area_seed(:,2) = round(area_seed_2);
end;
end;

if demo==1;
    tom_corr_subregions_handle=figure;set(gcf,'Position',[5         609        1392         510]);
    set(gcf,'Doublebuffer','on');
    m1=mean2(img_1);
    m2=mean2(img_2);
    subplot(1,3,1), imagesc(img_1',[m1-m1/4 m1+m1/4]); axis image; axis off; colormap gray;title('img1')
    subplot(1,3,2), imagesc(img_2',[m2-m2/4 m2+m2/4]); axis image; axis off; colormap gray;title('img2')
    subplot(1,3,3), imagesc(img_2',[m1-m1/4 m1+m1/4]); axis image; axis off; colormap gray;title('ccf')
end;
    
    

img_2_area_paste=zeros(size(img_1));
img_1_area_paste=zeros(size(img_2));
s1=size(img_1,1)./2;

img_1=tom_norm(img_1,'mean0+1std');
img_1_fft=fft2(single(img_1));
img_2=tom_norm(img_2,'mean0+1std');
img_2_fft=fft2(single(img_2));

% run for shifts
idx=1;

% statistic structure
statistic.area_seed=area_seed;

hit=1;
nohit=1;
ccfs=0;
for i=1:area_nr
    img_2_area=img_2(area_seed(i,1)-area_size_2+1:area_seed(i,1)+area_size_2,area_seed(i,2)-area_size_2+1:area_seed(i,2)+area_size_2);
    s2=size(img_2_area,1)./2;

    if demo==1;
        subplot(1,3,2);
        hold on;
        if exist('rect_ccf','var')
            set(rect_ccf,'EdgeColor',[0 0 0]);
        end;
        rect_ccf=rectangle('Position', [area_seed(i,1)-s2 area_seed(i,2)-s2 size(img_2_area,1) size(img_2_area,2)]);
        set(rect_ccf,'EdgeColor',[1 0 0]);
        hold off;
    end;
    img_2_area_paste(s1-s2+1:s1+s2,s1-s2+1:s1+s2)=img_2_area-mean2(img_2_area);
    img_2_area_paste_fft=single(fft2(img_2_area_paste));
    ccf=corr(single(img_2_area_paste_fft), single(img_1_fft), method,s2)./(s1.*2).^2;
    [me,ma,mi,st,va]=tom_dev(ccf,'noinfo');
    [c val]=tom_peak(ccf,'spline');
%    [c val]=tom_peak(ccf);
    c=c-1;
    shift_1_x=area_seed(i,1)-c(1);
    shift_1_y=area_seed(i,2)-c(2);
    statistic.shift_forward(i,1)=shift_1_x;
    statistic.shift_forward(i,2)=shift_1_y;
%    disp(['X_d: ' num2str(shift_1_x) ' Y_d: ' num2str(shift_1_y) ' Length: ' num2str(sqrt(shift_1_x^2+shift_1_y^2))]);            
    if demo==1
        subplot(1,3,3), imagesc(ccf',[me-st me+st]);axis image; axis off; colormap gray;
        title(['Max:' num2str(ma) ' Mean:' num2str(me) ' Std:' num2str(st)]); drawnow;
        hold on;
        p=plot(c(1),c(2),'ro');set(p,'LineWidth',5);
        hold off;
        subplot(1,3,1);
        hold on;
        if exist('rect_img1','var')
            if ~get(rect_img1,'EdgeColor')==[0 1 0];
                set(rect_img1,'EdgeColor',[0 0 0]);
            end;
        end;
        rect_img1=rectangle('Position', [c(1)-s2 c(2)-s2 size(img_2_area,1) size(img_2_area,2)]);
        set(rect_img1,'EdgeColor',[1 1 0]);
        hold off;
    end;
    cr=round(c);
    img_1_area=img_1(cr(1)-area_size_2+1:cr(1)+area_size_2,cr(2)-area_size_2+1:cr(2)+area_size_2);
    img_1_area_paste(s1-s2+1:s1+s2,s1-s2+1:s1+s2)=img_1_area-mean2(img_1_area);
    img_1_area_paste_fft=single(fft2(img_1_area_paste));
    ccf_b=corr_back(single(img_1_area_paste_fft), single(img_2_fft), method,s2)./(s1.*2).^2;
    [me_b,ma_b,mi_b,st_b,va_b]=tom_dev(ccf_b,'noinfo');
    [c_b val_b]=tom_peak(ccf_b,'spline');
    %        [c_b val_b]=tom_peak(ccf_b);
    c_b=c_b-1;
    shift_2_x=c_b(1)-cr(1);
    shift_2_y=c_b(2)-cr(2);
    statistic.shift_back(i,1)=shift_2_x;
    statistic.shift_back(i,2)=shift_2_y;
    %        disp(['length_1: ' num2str(sqrt(shift_1_x^2+shift_1_y^2)) 'length_2: ' num2str(sqrt(shift_2_x^2+shift_2_y^2)) ]);

    %        if abs(sqrt(shift_1_x^2+shift_1_y^2)-sqrt(shift_2_x^2+shift_2_y^2))<2
    %            disp(['length_1: ' num2str(sqrt(shift_1_x^2+shift_1_y^2)) 'length_2: ' num2str(sqrt(shift_2_x^2+shift_2_y^2)) ]);
    %        end;

    if sqrt((c_b(1)-area_seed(i,1))^2 + (c_b(2)-area_seed(i,2))^2 )<2
        shift_12_x=mean([shift_1_x shift_2_x]);
        shift_12_y=mean([shift_1_y shift_2_y]);
        statistic.shift_hit(hit,1)=shift_12_x;
        statistic.shift_hit(hit,2)=shift_12_y;
        statistic.pos_hit(hit,1)=area_seed(i,1);
        statistic.pos_hit(hit,2)=area_seed(i,2);
        hit=hit+1;
    else
        shift_12_x=mean([shift_1_x shift_2_x]);
        shift_12_y=mean([shift_1_y shift_2_y]);
        statistic.shift_nohit(nohit,1)=shift_12_x;
        statistic.shift_nohit(nohit,2)=shift_12_y;
        statistic.pos_nohit(nohit,1)=area_seed(i,1);
        statistic.pos_nohit(nohit,2)=area_seed(i,2);
        nohit=nohit+1;
    end;

    if demo==1
        subplot(1,3,3), imagesc(ccf_b',[me_b-st_b me_b+st_b]);axis image; axis off; colormap gray;
        title(['Backmatch: Max:' num2str(ma_b) ' Mean:' num2str(me_b) ' Std:' num2str(st_b)]); drawnow;
        hold on;
        plot(c_b(1),c_b(2),'ro');
        hold off;
        if sqrt((c_b(1)-area_seed(i,1))^2 + (c_b(2)-area_seed(i,2))^2 )<1
            if exist('rect_img1','var')
                set(rect_img1,'EdgeColor',[0 1 0]);
            end;
        end;
    end;
    %ccfs=ccfs+ccf;
end;
statistic.hit=hit-1;
statistic.nohit=nohit-1;

% clean hits

if statistic.hit~=0
    statistic.shift_hit_mean(1)=mean(statistic.shift_hit(:,1));
    statistic.shift_hit_mean(2)=mean(statistic.shift_hit(:,2));
    statistic.shift_hit_std(1)=std(statistic.shift_hit(:,1));
    statistic.shift_hit_std(2)=std(statistic.shift_hit(:,2));

    old_hit_nr=statistic.hit;

%    if statistic.shift_hit_std(1)~=0 && statistic.shift_hit_std(2)~=0 && statistic.nohit~=0
    statistic=tom_corr_subregions_clean_statistic(statistic,nr_std);
%    end;

    while old_hit_nr~=statistic.hit
%        if statistic.shift_hit_std(1)~=0 && statistic.shift_hit_std(2)~=0 
            old_hit_nr=statistic.hit;
            statistic=tom_corr_subregions_clean_statistic(statistic,nr_std);
            disp(['clean statistic with ' num2str(nr_std) 'std.'])
%        end;
    end;


    shifts=statistic.shift_hit_mean;
    shifts_std=statistic.shift_hit_std;
else
    % no hit!!!
    shifts=[];
    shifts_std=[];
end;


function ccf=corr(a_fft,b_fft,method, template_size_2)    

ccf=real(fftshift(ifft2( (conj(a_fft).*(b_fft)) )));
ccf(1:template_size_2,1:end)=0;
ccf(1:end,1:template_size_2)=0;
ccf(size(ccf,2)-template_size_2+2:end,1:end)=0;
ccf(1:end,size(ccf,2)-template_size_2+2:end)=0;

function ccf=corr_back(a_fft,b_fft,method, template_size_2)    

ccf=real(fftshift(ifft2( (conj(a_fft).*(b_fft)) )));
ccf(1:template_size_2,1:end)=0;
ccf(1:end,1:template_size_2)=0;
ccf(size(ccf,2)-template_size_2+1:end,1:end)=0;
ccf(1:end,size(ccf,2)-template_size_2+1:end)=0;


% function statistic=tom_corr_subregions_clean_statistic(statistic,nr_std)
%     
% hit=1;
% nohit=1;
% hits=statistic.shift_hit;
% nohits=statistic.shift_nohit;
% allhits=[hits; nohits];
% statistic.shift_nohit=[];
% statistic.shift_hit=[];
% pos_hit=statistic.pos_hit;
% statistic.pos_hit=[];
% pos_nohit=statistic.pos_nohit;
% statistic.pos_nohit=[];
% allpos=[pos_hit; pos_nohit];
% 
% for i=1:size(allhits,1)
%     devx_min=statistic.shift_hit_mean(1)-statistic.shift_hit_std(1).*nr_std;
%     devx_max=statistic.shift_hit_mean(1)+statistic.shift_hit_std(1).*nr_std;
% 
%     devy_min=statistic.shift_hit_mean(2)-statistic.shift_hit_std(2).*nr_std;
%     devy_max=statistic.shift_hit_mean(2)+statistic.shift_hit_std(2).*nr_std;
% 
%     
%     if allhits(i,1)<devx_min ||  allhits(i,1)>devx_max ...
%             || allhits(i,2)<devy_min || allhits(i,2)>devy_max 
%                 statistic.shift_nohit(nohit,:)=allhits(i,:);
%                 statistic.pos_nohit(nohit,:)=allpos(i,:);
%                 nohit=nohit+1;
%     else
%                 statistic.shift_hit(hit,:)=allhits(i,:);
%                 statistic.pos_hit(hit,:)=allpos(i,:);
%                 hit=hit+1;
%     end;
% end;
% statistic.hit=hit-1;
% statistic.nohit=nohit-1;
% statistic.shift_hit_mean(1)=mean(statistic.shift_hit(:,1));
% statistic.shift_hit_mean(2)=mean(statistic.shift_hit(:,2));
% statistic.shift_hit_std(1)=std(statistic.shift_hit(:,1));
% statistic.shift_hit_std(2)=std(statistic.shift_hit(:,2));

function area_seed=extract_area_seed_by_features(img,area_nr,template_size_2)

img_fil=tom_bandpass(edge(img),2,size(img,1)./16);
img_fil(1:template_size_2,1:end)=0;
img_fil(1:end,1:template_size_2)=0;
img_fil(size(img,2)-template_size_2+1:end,1:end)=0;
img_fil(1:end,size(img,2)-template_size_2+1:end)=0;

for i=1:area_nr;
    [c val img_fil] = tom_peak(img_fil,template_size_2./2);
    area_seed(i,1)=c(1);
    area_seed(i,2)=c(2);
end;
disp('area_seed created')

function area_seed=generate_area_seed_along_tiltaxis(img_2,area_size_2,area_nr,tiltaxis,area_seed_region);

l=zeros(size(img_2),'single');
s1=size(img_2,1);

area=round((s1-2).*area_seed_region);

l(s1./2-area:s1./2+area,:)=1;
lr=tom_rotate(single(l),[tiltaxis])>0.5;
lr(1:area_size_2,1:end)=0;
lr(1:end,1:area_size_2)=0;
lr(size(lr,2)-area_size_2+1:end,1:end)=0;
lr(1:end,size(lr,2)-area_size_2+1:end)=0;

i=1;
rand('twister',sum(100*clock));
 while i<area_nr+1
    t=round(rand(1,2).*(size(img_2)-1))+1;
    if lr(t(1),t(2))>0 
        area_seed(i,:)=t;
        i=i+1;
    end;
 end;
 disp('area_seed created');



