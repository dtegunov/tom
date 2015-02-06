function [bfactor decay_restore corrected]=tom_fit_bfactor(in, objectpixelsize, fit_range, apply_range, mass_Da,FSC,info_flag)

% [bfactor decay_restore corrected]=tom_fit_bfactor(in, objectpixelsize, fit_range, apply_range, mass_Da,FSC,info_flag)
%
% Input:
%           in:                         input 3D volume.
%           objectpixelsize:            objectpixelsize in Angstrom
%           fit_range:                  fit range for bfactor determination in Angstrom.
%           apply_range:                apply bfactor in this resolution
%                                       range. Accepts 'Inf' as upper
%                                       boundary. Lower boundary is
%                                       cut-off, higher frequencies are set
%                                       to zero.
%           mass_Da:                    molecular mass, not used at the moment.
%           FSC:                        Fourier Shell Curve, not used at the moment.
%           info_flag:                  verbose flag, show Guinier graph and fit.
% Output:
%           bfactor:                    determined Bfactor in 1/(A^2).
%           decay_restore:              restored decay function.
%           corrected:                  input volume corrected by determined Bfactor.
% 
% Example:
% [bfactor decay_restore corrected]=tom_fit_bfactor(vol_1.Value,2,[5 20],[5 Inf],2500000,1,1);
%
% Theoretical background:
%
% The Debye-Waller factor (DWF), named after Peter Debye and Ivar Waller,
% is used in condensed matter physics to describe the attenuation of x-ray
% scattering or neutron scattering caused by thermal motion or quenched disorder.
% It has also been called the B-factor or the temperature factor.
% The DWF depends on q, the absolute value of the scattering vector q.
% For a given q, DWF(q) gives the fraction of elastic scattering; 1-DWF(q)
% correspondingly the fraction of inelastic scattering. In diffraction studies,
% only the elastic scattering is useful; in crystals, it gives rise to distinct
% Bragg peaks. Inelastic scattering events are undesirable as they cause a diffuse
% background.
% 1. Debye, Peter (1913). "Interferenz von Roentgenstrahlen und Waermebewegung"
% (in German). Ann. d. Phys. 348 (1): 49-92.
% 2. Waller, Ivar (1923). "Zur Frage der Einwirkung der Waermebewegung auf
% die Interferenz von Roentgenstrahlen" (in German).
% Zeitschrift fuer Physik A Hadrons and Nuclei (Berlin/Heidelberg: Springer)
% 17: 398-408. 
%
% Identification and interpretation of structural features of EM
% density maps are hampered by the loss of high resolution contrast.
% Degradation of the image contrast is caused, by factors related to the experimental
% imaging process (e.g. specimen movement and charging, radiation damage, partial
% microscope coherence, etc.) (Henderson, 1992; Wade, 1992)., as
% well as by factors related to the computational procedures for
% structure determination (e.g. inaccurate determination of the orientation
% parameters, etc.) (Conway and Steven, 1999; Rosenthal
% and Henderson, 2003; Henderson, 2004).
% The combined effect of all these factors has been traditionally modeled by a Gaussian
% amplitude decay of structure factors, given by B_overall=exp(B/(4*d^2)) with an
% overall temperature factor B_overall (also called B-factor) and where
% d denotes resolution (Glaeser and Downing, 1992; Conway
% and Steven, 1999; Rosenthal and Henderson, 2003). This amplitude
% fall-off significantly affects the high resolution components, thereby
% making the density map look apparently smooth.
% For theoretical cosniderations see 3.), for a more recent implementation
% see 4.)
% 3. "Optimal Determination of Particle Orientation, Absolute Hand, and Contrast Loss in
% Single-particle Electron Cryomicroscopy." P.B. Rosenthal and R. Henderson
% Journal of Molecular Biology, 333:721-745, 2003.
% 4. "Sharpening high resolution information in single particle electron
% cryomicroscopy.", J.J. Fernandez, D. Luque, J.R. Caston, J.L. Carrascosa
% Journal of Structural Biology, 164:170-175, 2008.
%
%
%REFERENCES
% 1. Debye, Peter (1913). "Interferenz von Roentgenstrahlen und Waermebewegung"
%   (in German). Ann. d. Phys. 348 (1): 49-92.
% 2. Waller, Ivar (1923). "Zur Frage der Einwirkung der Waermebewegung auf
%   die Interferenz von Roentgenstrahlen" (in German).
%   Zeitschrift fuer Physik A Hadrons and Nuclei (Berlin/Heidelberg: Springer)
%   17: 398-408.
% 3. "Optimal Determination of Particle Orientation, Absolute Hand, and Contrast Loss in
%   Single-particle Electron Cryomicroscopy." P.B. Rosenthal and R. Henderson
%   Journal of Molecular Biology, 333:721-745, 2003.
% 4. "Sharpening high resolution information in single particle electron
%   cryomicroscopy.", J.J. Fernandez, D. Luque, J.R. Caston, J.L. Carrascosa
%   Journal of Structural Biology, 164:170-175, 2008.
%
%SEE ALSO
%   tom_apply_bfactor
%
%   created by SN 05/06/10
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

if nargin<7
    info_flag=0;
end;
if nargin<6
    FSC=1.0;
end;
zero_angle_scattering=0.28.*(mass_Da./12); % 12 for carbon, 700 kDa, 700000/12=N_atoms
wilson_regime=sqrt(mass_Da./12); 
if info_flag==1
disp(['------------------------------------------------------------------------------------']);
disp(['EM density with:']);
disp(['Objectpixelsize: ' num2str(objectpixelsize) ' Angstrom.']);
disp(['Molecular Mass: ' num2str(mass_Da) ' Da.']);
disp(['Bfactor Fit Range: ' num2str(fit_range(1)) ' nm to ' num2str(fit_range(2)) ' nm.']);
disp(['Bandpass Filter: ' num2str(apply_range(1)) ' nm to ' num2str(apply_range(2)) ' nm.']);
if nargin<6
    disp(['Fourier Shell Correlation Curve constant 1.0 .']);
    FSC=1.0;
else
   disp(['Fourier Shell Correlation Curve provided and applied.']);
    if length(FSC)==1
        disp(['FSC: ' num2str(FSC)]);
    end;
end;
disp(['------------------------------------------------------------------------------------']);
disp(['Contrast correction:']);
disp(['Wilson regime: ' num2str(wilson_regime) ', Zero Angle Scattering: ' num2str(zero_angle_scattering)]);
disp('calculate Guinier graph and fit line ...');
end;
% Experimental scattering amplitudes may be placed on an absolute scale by setting
% the zero angle scattering equal to 0.28 x N_atoms and the average
% amplitude in the high resolution region (Wilson regime) to
% sqrt(N_atoms). N_atoms denotes the number of carbon atom equivalents
% corresponding to the molecular mass of the protein, and 0.28 is the solvent contrast factor.
% Ref 4.), page 171, top, right.

bfactor=0;
corrected=0;
fit_range_pixel=(2.*objectpixelsize./fit_range).*size(in,1)./2;
sz=size(in,1);
sz_2=size(in,1)./2;

lln=calc_fourier_shell(sqrt(tom_ps(in)));

% normalize decay to 0.28*N_atoms, zero_angle_scattering
lln=lln./max(lln);
lln=lln.*zero_angle_scattering;

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
if info_flag==1
    plot_decay_curve(decay,objectpixelsize,decay_area,x_area_plot,fitline,q_square_area,q_square,zero_angle_scattering,wilson_regime,FSC); drawnow;
    disp(['Line Fit: y=a+b*x, a: ' num2str(var(2)) ', b: ' num2str(var(1))]);
    disp(['B-factor: ' num2str(bfactor)]);
    disp('correct contrast ...');
end;

[corrected decay_restore decay_restore_3d]=tom_apply_bfactor(in,objectpixelsize,bfactor,FSC,apply_range);

if info_flag==1
    disp('done.');
end;

function [l nr]=calc_fourier_shell(in)

l=0;
sz=size(in,1);
sz_2=sz./2+1;
% for r=1:sz_2
%     mask_1=tom_spheremask(ones(sz,sz,sz),r,0,[sz_2+1 sz_2+1 sz_2+1]);
%     mask_2=tom_spheremask(ones(sz,sz,sz),r+1,0,[sz_2+1 sz_2+1 sz_2+1]);
%     mask=mask_2-mask_1;
%     s=sum(sum(sum(in.*mask)));
%     p=sum(sum(sum(mask)));
%     l(r)=s./p;
%     %tom_dspcub(mask);drawnow;
% end;
idx=1;
[x,y,z]=ndgrid(0:size(in,1)-1,0:size(in,2)-1,0:size(in,3)-1);
v = sqrt((x+1-sz_2).^2+(y+1-sz_2).^2+(z+1-sz_2).^2);
for r=0:sz_2-1
    ind = find(round(v)==r);    
    l(idx)=mean(in(ind));
    nr(idx)=length(ind);
    idx=idx+1;
end

function plot_decay_curve(decay,objectpixelsize,decay_area,x_area_plot,fitline,q_square_area,q_square,zero_angle_scattering,wilson_regime,FSC)
x=([0:q_square./64:q_square]);
fig=figure;set(fig,'Position',[680   487   802   613]);
drawnow;
FSC_weight=sqrt((2.*FSC)./(1+FSC));
p1=plot(q_square,decay,'Linewidth',2); hold on;
p2=plot(q_square_area,decay_area,'r-','Linewidth',2);
p3=plot(q_square_area,fitline,'g-','Linewidth',2);
p4=plot(q_square,FSC_weight.*max(decay),'m-','Linewidth',2); hold on;
p=get(p1,'Parent');
%set(fig,'Position',[80 12 60 20]);
%set(p,'Ylim',[(min(decay)-1) (max(decay)+1)]);
%set(t,'Color','white');
xlabel('Resolution, 1/q^2 [1/(A^2)]');
set(p,'Ytick',[round((min(decay)-1)):1:round((max(decay)+1))]);
%set(p,'Xtick',[0:x(end)./10:x(end)]);
ylabel('log(spherically averaged structure factor curve)');
grid on;
p=get(p1,'Parent');
%xtick=get(p,'Xtick');
%xtick_nm=zeros(size(xtick));
%xtick_nm(2:end)=2.*objectpixelsize.*max(xtick).*(1./xtick(2:end));

%for i=2:size(xtick,2)
%    t=text(xtick(i),-5,sprintf('%0.3g A',xtick_nm(i)));    
%end;

plot(0,log(zero_angle_scattering),'rx','Linewidth',8);
tn=text(0,log(zero_angle_scattering),'0.28*N_a_t_o_m_s');
set(tn,'Color','red');set(tn,'Fontsize',14);
%l=line([0 q_square(end)],[log(wilson_regime) log(wilson_regime)]);
wilson_x=[0 q_square(end)];
wilson_y=[log(wilson_regime) log(wilson_regime)];
l=plot(wilson_x,wilson_y,'c-');
set(l,'Color','cyan','Linewidth',2);
legend('Guinier graph','Fit area','Line fit','FSC weight scaled','Zero Scattering','Wilson regime');
title(['Guinier Graph, Spherically Averaged Structure Factor Curve, Nyquist @ ' num2str(objectpixelsize.*2) ' A.']);
hold off;
