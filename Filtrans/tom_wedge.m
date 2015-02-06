function wedge=tom_wedge(image,angle,rot_angle)
%TOM_WEDGE produces a wedge shaped array
%
%   wedge=tom_wedge(image,angle)
%
%PARAMETERS
%
%  INPUT
%   image               input array - 3d
%   angle               semi angle of missing wedge in deg
%   rot_angle           (optional) angle the wedge is rotated
%  
%  OUTPUT
%   wedge               output - filter
%
%   TOM_WEDGE produces a wedge shaped array
%   This array can be used as a window filter in Fourier space...
%
%EXAMPLE
%   yyy = zeros(64,64,64);
%   wedge=tom_wedge(yyy,30);
%   tom_dspcub(wedge,1);
%   yyy(1,1,1) =1;
%   psf = real(tom_ifourier(ifftshift(fftshift(tom_fourier(yyy)).*wedge)));
%   figure;tom_dspcub(psf); % creates PSF of missing wedge
%
%   dual axis:
%   yyy = zeros(32,32,32);
%   wedge=tom_wedge(yyy,30);
%   wedge_rot_90=tom_rotate(wedge,90);
%   dual_wedge=(wedge_rot_90+wedge)==0;
%   dual_wedge=abs(dual_wedge-1);
%   yyy(1,1,1) =1;
%   psf_dual_wedge = real(tom_ifourier(ifftshift(fftshift(tom_fourier(yyy)).*dual_wedge)));
%   figure;tom_dspcub(psf_dual_wedge); % creates PSF of missing wedge
%
%
%   from particles (test):
%   ps=0;for i=1:50;in=tom_emread(['particle_' num2str(i) '.em']);ps=tom_ps(in.Value)+ps;end;
%   yyy = zeros(40,40,40);
%   yyy(1,1,1) =1;
%   psf = real(tom_ifourier(ifftshift(fftshift(tom_fourier(yyy)).*wedge)));
%
%REFERENCES
%
%SEE ALSO
%   TOM_FILTER, TOM_BANDPASS
%
%   created by FF 07/20/03
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


if nargin<3
    rot_angle=[0 0 0];
end;

warning off MATLAB:divideByZero;
angle = angle*pi/180;
[dimx, dimy, dimz] = size(image);
[x,y,z] = ndgrid(-floor(dimx/2):-floor(dimx/2)+dimx-1,-floor(dimy/2):-floor(dimy/2)+dimy-1,-floor(dimz/2):-floor(dimz/2)+dimz-1);
wedge = ones(dimx, dimy, dimz);
ind = find(tan(angle) > abs(x)./abs(z));
wedge(ind)=0;


if ((sum(rot_angle==[0 0 0])==3)==0 )
    wedge=tom_rotate(wedge,rot_angle);
end;

