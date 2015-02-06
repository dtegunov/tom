function y = tom_ctf(Dz, pix_size, voltage, pixs, Cs, alpha, Cc, deltaE, infoflag)
%TOM_CTF calculates and plots CTF (pure phase contrast)
%
%   ctf_out = tom_ctf(Dz, pix_size, voltage, pixs, Cs, alpha, Cc, deltaE)
%
%   The one-dimensional CTF is calculated and plotted. In case no figure is
%   displayed open one (see also example). The function can be used to get
%   a quick overview of the CTF for cterian imaging conditions.
%
%PARAMETERS
%  INPUT
%   Dz       : Defocus (<0 underfocus, >0 overfocus) (in \mu m);
%   pix_size : pixel size (in nm) (default: 0.72 nm)
%   voltage  : accelerating Voltage (in keV) (default: 300 kV)
%   pixs     : Number of pixels used (default: 2048)
%   Cs       : sperical aberration in mm(default: 2 mm)
%   alpha    : illumination aperture in mrad (default 0.02 - reasonable for FEG)
%   Cc       : chrmatic aberration in mm -(default 2.2 mm)
%   deltaE   : energy width in eV (default 0.8 eV)
%   infoflag : display axes information (default 1)
%
%  OUTPUT
%   ctf_out  : vector of dim pixs/2 containig the plotted ctf (if pixs is 
%                                    not even: dim = ceil(pixs/2))
%              the ouput is beeing plotted where the absyssis are inverse
%              pixels.
%EXAMPLE
%   figure;
%   ctf = tom_ctf(-4,1, 120,1024,2);
%
%REFERENCES
%
% SEE ALSO
%    TOM_CREATE_CTF TOM_FOURIER TOM_IFOURIER 
%
%   created by FF 09/10/02
%   updated by ...
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom%


error(nargchk(1,9,nargin));

if nargin<9
    infoflag = 1;
end

if nargin<8
  deltaE=0.8; 
end;
if nargin<7
  Cc=2.2/1000; 
else
   Cc=Cc/1000; 
end;
if nargin<6
  alpha=0.02/1000; 
else
   alpha=alpha/1000; 
end;
if nargin<5
  Cs=2*10^(-3); 
else
   Cs=Cs*10^(-3); 
end;
if nargin<4
    pixs=2048;
end;
if nargin<3
    voltage=300000; 
else
    voltage = voltage * 1000;
end;
if nargin <2
    pix_size=0.72*10^(-9);
else
    pix_size=pix_size*10^(-9);
end;
Dz=Dz*10^(-6); %from micrometer to meter
Dzn=Dz*1000000; %for display
Csn=Cs*1000;%for display
voltagen=voltage/1000;%for display
voltagest=voltage*(1+voltage/1022000); %for relativistic calc
lambda=sqrt(150.4/voltagest)*10^-10;
q=0:1/(pixs*pix_size):1/(2*pix_size);% von, Increment, Nyqvist
nyqvist = 2*pix_size*10^9;
y = sin( pi/2* (Cs*lambda.^3.*q.^4 - 2*Dz*lambda*q.^2) );
% new function: include also envelope function:
% 1) spatial coherence
%  alpha: illumination aperture - assume 0.02mrad
Ks = exp(- ((pi*(Cs.*lambda.^2.*q.^3 - Dz.*q)*alpha).^2) / log(2) );
% 2) temporal coherence
delta = Cc*deltaE/voltage;
%Kt = exp(- (pi*lambda*delta*q.^2/2));% alter according to Reimer
Kt = exp(- (pi*lambda*delta*q.^2/(4*log(2))).^2);
K = Kt.*Ks;
%1st zero (approximately for high defocus)
ctf_zero = sqrt(lambda*abs(Dz)*10^18);
plot(0:floor(pixs/2),K.*y,'LineWidth',1.5);
hold on; plot(0:floor(pixs/2),K,'r','LineWidth',1.5);hold off;
grid on;
axis([0 floor(pixs/2) -1.1 1.1] );

if infoflag == 1
    title(['CTF for: Defocus:',num2str(Dzn),'\mum, Voltage:',num2str(voltagen),'kV, C_s:',num2str(Csn),...
            'mm, Nyqvist:', num2str(nyqvist),'nm, 1st zero:', num2str(ctf_zero,2),'nm']);
    ylabel('CTF');
    xlabel('Frequency');
end
