function [fsc f v_05 v_03 v_01]=tom_fsc(vol1,vol2,nr_shells,objectpixelsize,display,wedge)

%TOM_FSC calculates a Fourier Shell Correlation Function of two 3D volumes.
%
%   [fsc f v_05 v_03 v_01]=tom_fsc(vol1,vol2,nr_shells,objectpixelsize,display,wedge)
%
%PARAMETERS
%
%  INPUT
%   vol1, vol2          input volumes
%   nr_shells           number of shells in Fourier space
%   objectpixelsize     object pixelsize of volumes
%   display             0/1 flag for plot display
%   wedge               missing wedge of volumes, optional
%
%  OUTPUT
%   fsc                 Fourier Shell Correlation function
%   f                   frequencies
%   v_05, v_03, v_01    0.5, 0.3 and 0.1 threshould criteria
%
%
%EXAMPLE
%   Example:
%   [fsc f]=tom_fsc(vol1,vol2,32,3.6,1);
%   % with missing wedge
%   wedge=f3d_create_wedge([64 64 64],[-60 60]);
%   [fsc_wedge f v05 v03 v01]=tom_fsc(vol1,vol2,32,3.6,1,wedge);
%
%REFERENCES
%
%SEE ALSO
%
%   created by SN, 12/27/10
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

if nargin < 2
    error('tom_fsc:argChk','Wrong number of input arguments\nSyntax: [fsc]=tom_fsc(vol1,vol2,nr_shells,objectpixelsize,display).')
end
if size(vol1)~=size(vol2)
    error('tom_fsc:argChk','Input volumes have different dimensions.\nSyntax: [fsc]=tom_fsc(vol1,vol2,nr_shells,objectpixelsize,display).');
end
if nargin<6
    wedge=1;
end
if nargin<5
    display=0;
end
if nargin<4
    display=0;objectpixelsize=1;
end
if nargin<3
    nr_shells=size(vol1,1);
end
if nargin>2
    nr_shells=nr_shells.*2;
end
vol1_fft = fftshift(fftn(vol1,[nr_shells nr_shells nr_shells]));
vol2_fft_conj = fftshift(conj(fftn(vol2,[nr_shells nr_shells nr_shells])));
sz=size(vol1_fft,1);
sz_2=floor(sz./2);
vol1_fft =vol1_fft./(sz*sz*sz);
vol2_fft_conj =vol2_fft_conj./(sz*sz*sz);
[x,y,z]=ndgrid(1:sz,1:sz,1:sz);
v = sqrt((x-(sz_2+1)).^2+(y-(sz_2+1)).^2+(z-(sz_2+1)).^2);
v = v.*wedge;
fsc=zeros(1,sz_2);
idx=1;
for r=0:sz_2-1
    ind = find(floor(v)==r);
    nr_idx=length(ind);
    if nr_idx>1
        fsc_n=sum((vol1_fft(ind)-mean(vol1_fft(ind))).*(vol2_fft_conj(ind)-mean(vol2_fft_conj(ind))))./(nr_idx-1);
        fsc(idx)=real(fsc_n./(std(vol1_fft(ind)).*std(vol2_fft_conj(ind))));
    else
        fsc(idx)=1;
    end;
    idx=idx+1;
end
f=1./(length(fsc)./[1:length(fsc)].*(objectpixelsize.*2));
v_03=0;
v_05=0;
v_01=0;
for i=1:length(fsc)
    if fsc(i)<0.3
        v1=f(i-1);
        v2=f(i);
        v=(v1+v2)./2;
        v_03=1/v;
        break;
    end;
end;
for i=1:length(fsc)
    if fsc(i)<0.5
        v1=f(i-1);
        v2=f(i);
        v=(v1+v2)./2;
        v_05=1/v;
        break;
    end;
end;
for i=1:length(fsc)
    if fsc(i)<0.1
        v1=f(i-1);
        v2=f(i);
        v=(v1+v2)./2;
        v_01=1/v;
        break;
    end;
end;

if display==1
    figure;
    ny = objectpixelsize.*2;
    p1=plot(f,fsc,'LineWidth',2,'Color','black');
    grid on;
    p=get(p1,'Parent');
    set(p,'Ylim',[(floor(min(fsc).*10))./10 1.1]);
    set(p,'Xlim',[min(f) max(f)]);
    xlabel('Resolution [1/Ang]');
    ylabel('FSC');    
    set(p,'Ytick',(floor(min(fsc).*10))/10:0.1:1.1);
    set(p,'Xtick',(0:max(f)/4:max(f)));
    xtick=get(p,'Xtick');
    xtick_A=zeros(size(xtick));
    xtick_A(2:end)= ((xtick(end)./xtick(2:end)*ny));
    for i=2:size(xtick,2)
        text(xtick(i),(floor(min(fsc).*10))./10+0.02,sprintf('%0.3g A',xtick_A(i)));
    end;
    text(f(2)./2,0.3,sprintf('  %0.3g A >>>',v_03));
    text(f(2)./2,0.5,sprintf('  %0.3g A >>>',v_05));
    title(['Fourier Shell Correlation, Nyquist @ ' num2str(ny) ' Ang. FSC: 0.5 @ ' sprintf('%0.3g',v_05) ' A, 0.3 @ ' sprintf('%0.3g',v_03) ' A.']);
end