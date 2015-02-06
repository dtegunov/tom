function [vol, amplitude]= tom_create_ctf(Dz, vol, pix_size, voltage, Cs, sigma)
%TOM_CREATE_CTF calculates 2D or 3D CTF (pure phase contrast)
%
%   Note: only tested for even dimensions!
%
%   [ctf_out amplitude] = tom_create_ctf(Dz, vol, pix_size, voltage, Cs, sigma)
%
%PARAMETERS
%  INPUT
%   Dz       : Defocus (<0 underfocus, >0 overfocus) (in \mu m);
%   vol      : Volume (or image)
%   pix_size : pixel size (in nm) (default: 0.72 nm)
%   voltage  : accelerating Voltage (in keV) (default: 300 kV)
%   Cs       : sperical aberration
%   Sigma    : envelope of ctf (optional). If chosen decay 
%              ctf ~exp(-(freq/sigma)^2) is assumed. SIGMA is in units of Nyqvist.
%              => ctf(sigma) = 1/e 
%
%  OUTPUT
%   ctf_out  : output containing the centrosymmetric ctf. It can be used as
%              a fourier filter.
%
%EXAMPLE
%   im = tom_emread('proteasome.em');
%   ctf = tom_create_ctf(-4.4,im.Value,1.0240, 120,2);
%   tom_imagesc(ctf)
%
%REFERENCES
%
%SEE ALSO
%   TOM_CTF TOM_FOURIER TOM_IFOURIER 
%
%
%   created by FF 09/14/04
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

error(nargchk(2,6,nargin));

if nargin<5
  Cs=2*10^(-3); 
else
   Cs=Cs*10^(-3); 
end;
if nargin<4
    voltage=300000; 
else
    voltage = voltage * 1000;
end;

if nargin <3
    pix_size=0.72*10^(-9);
else
    pix_size=pix_size*10^(-9);
end;
Dz=Dz*10^(-6); 
Dzn=Dz*1000000; %for display
Csn=Cs*1000;%for display
voltagen=voltage/1000;%for display
voltagest=voltage*(1+voltage/1022000); %for relativistic calc

lambda=sqrt(150.4/voltagest)*10^-10;


Ny = 1/(2*pix_size);
nyqvist = 2*pix_size*10^9;
disp(['CTF is calculated for: Defocus' num2str(Dzn) '\mum Voltage = ' num2str(voltagen) 'kV, Nyqvist = ' num2str(nyqvist) 'nm'  ]);
vol = 0*vol; %empty
[r y z] = meshgrid(-Ny:2*Ny/size(vol,1):Ny-Ny/size(vol,1), -Ny:2*Ny/size(vol,2):Ny-Ny/size(vol,2), -Ny:2*Ny/size(vol,3):Ny-Ny/size(vol,3) );
if size(vol,3)>1
    r = sqrt(r.^2+y.^2+z.^2);clear y; clear z;
else
    r = sqrt(r.^2+y.^2);clear y; clear z;
end;
vol = sin(pi/2*(Cs*lambda^3*r.^4 - 2*Dz*lambda*r.^2));
amplitude = cos(pi/2*(Cs*lambda^3*r.^4 - 2*Dz*lambda*r.^2));
if nargin==6
    vol = vol.*exp(-(r./(sigma*Ny)).^2);
    amplitude = amplitude.*exp(-(r./(sigma*Ny)).^2);
end;
clear r;

