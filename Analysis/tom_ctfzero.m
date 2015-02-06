function qzero= tom_ctfzero(Dz, pix_size, voltage, pixs, Cs)
%TOM_CTFZERO calculates zeros of CTF (pure phase contrast)
%
%   qzero = tom_ctfzero(Dz, pix_size, voltage, pixs, Cs)
%
%PARAMETERS
%  INPUT
%   Dz       : Defocus (<0 underfocus, >0 overfocus) (in \mu m);
%   pix_size : pixel size (in nm) (default: 0.72 nm)
%   voltage  : accelerating Voltage (in keV) (default: 300 kV)
%   pixs     : Number of pixels used (default: 2048)
%   Cs       : sperical aberration
%
%  OUTPUT
%   qzero     : vector of dim 20 containing the first 20 zeros of the CTF.
%
%EXAMPLE
%   qzero = tom_ctfzero(-4,1, 120,1024,2);
%
%REFERENCES
%
% SEE ALSO
%    TOM_CREATE_CTF TOM_FOURIER TOM_IFOURIER TOM_CTFFIT
%
%   created by FF 04/22/03
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

error(nargchk(1,5,nargin));

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
Dz=Dz*10^(-6); 
Dzn=Dz*1000000; %for display
Csn=Cs*1000;%for display
voltagen=voltage/1000;%for display
voltagest=voltage*(1+voltage/1022000); %for relativistic calc
lambda=sqrt(150.4/voltagest)*10^-10;
n=1:20;
% calculate 1st 20 zeros
qzero = sqrt(sqrt(2.*n./(Cs.*lambda^3)+(Dz/(Cs*lambda^2))^2) + Dz/(Cs*lambda^2));
% scale to sampling
qzero = (pixs*pix_size).*qzero;
