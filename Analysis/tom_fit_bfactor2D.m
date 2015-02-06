function [bfactor bfit]=tom_fit_bfactor2D(in, objectpixelsize, fit_range,ps_size,info_flag)

% [bfactor ]=tom_fit_bfactor2D(in, objectpixelsize, fit_range, info_flag)
%
% Input:
%           in:                         input 3D volume.
%           objectpixelsize:            objectpixelsize in Angstrom
%           fit_range:                  fit range for bfactor determination in Angstrom.
%           ps_size:                    size of calculated power spectrum.
%           info_flag:                  verbose flag, show Guinier graph and fit.
% Output:
%           bfactor:                    determined Bfactor in 1/(A^2).
%
% Example:
%   [bfactor]=tom_fit_bfactor2D(img_2.Value, 1.018, [3 8],256,1);
%
%   % with bfit, fit-structure
%   [bfactor bfit]=tom_fit_bfactor2D(img_2.Value, 1.018, [3 8],256,1);
%
%
%SEE ALSO
%   tom_fit_bfactor
%
%   created by SN 02/16/11
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2010
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

if nargin<5
    info_flag=0;
end;
if nargin<4
    ps_size=256;
end;
if info_flag==1
    disp(['------------------------------------------------------------------------------------']);
    disp(['Electron micrograph with:']);
    disp(['Objectpixelsize: ' num2str(objectpixelsize) ' Angstrom.']);
    disp(['Bfactor Fit Range: ' num2str(fit_range(1)) ' Angstrom to ' num2str(fit_range(2)) ' Angstrom.']);
    disp('calculate Guinier graph and fit line ...');
end;

bfit.objectpixelsize=objectpixelsize;
bfit.fit_range=fit_range;
bfit.ps_size=ps_size;


bfactor=0;
fit_range_pixel=(2.*objectpixelsize./fit_range).*ps_size./2;
sz=size(in,1);
sz_2=sz./2;
ps=fftshift(tom_calc_periodogram_parallel(single(in),ps_size,0,floor(ps_size./16)));
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

q_square=1./((x.^2));
q_square_area=1./((x_area.^2));
var = polyfit(q_square_area,decay_area,1);
a = var(2); % y-intercept of the fitted line
b = var(1); % slope of fitted line
fitline = a+b*q_square_area;
bfactor=4.*var(1);
bfit.bfactor=bfactor;

if info_flag==1
    plot_decay_curve(decay,objectpixelsize,decay_area,x_area_plot,fitline,q_square_area,q_square,ps_size); drawnow;
    disp(['Line Fit: y=a+b*x, a: ' num2str(var(2)) ', b: ' num2str(var(1))]);
    disp(['B-factor: ' num2str(bfactor)]);
end;

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

function plot_decay_curve(decay,objectpixelsize,decay_area,x_area_plot,fitline,q_square_area,q_square,ps_size)
x=([0:q_square./64:q_square]);
fig=figure;set(fig,'Position',[680   487   802   613]);
drawnow;
p1=plot(q_square,decay,'Linewidth',2); hold on;
p2=plot(q_square_area,decay_area,'r-','Linewidth',2);
p3=plot(q_square_area,fitline,'g-','Linewidth',2);
p=get(p1,'Parent');
xlabel('Resolution, 1/q^2 [1/(A^2)]');
set(p,'Ytick',[round((min(decay)-1)):1:round((max(decay)+1))]);
ylabel('log(spherically averaged structure factor curve)');
grid on;
p=get(p1,'Parent');
legend('Guinier graph','Fit area','Line fit');
title(['Guinier Graph, Spherically Averaged Structure Factor Curve, Nyquist @ ' num2str(objectpixelsize.*2) ' A. Power Spectrum size: ' num2str(ps_size) ' pixel.']);
hold off;
