function pol = tom_cart2sph(I)
%TOM_CART2SPH transforms 3D-volumes from cartesian to spherical coordinates
%
%   pol = tom_cart2sph(I)
%
%   A 3D volume I given in cartesian coordinates is sampled in spherical coordinates 
%   R, PHI, THETA using trilinear interpolation. The input image I is assumed to be of 
%   equal dimensions NX, NY, and NZ. The ouput POL(R,PHI,THETA) has dimensions NR=NX/2, 
%   NPHI=4*NR, NTHETA=2*NR.
%
%PARAMETERS
%
%  INPUT
%   I                   3 dim array
%  
%  OUTPUT
%   pol                 3 dim array in spherical coordinates
%
%EXAMPLE
%   cart = zeros(32,32,32);
%   cart(8,8,17) = 1;
%   cart = tom_symref(cart,4);
%   sph = tom_cart2sph(cart);
%   tom_dspcub(sph);
%
%REFERENCES
%
%SEE ALSO
%   TOM_POLAR2CART, TOM_CART2POLAR, TOM_SPH2CART
%
%   created by FF 08/17/03
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


nx = size(I,1);ny = size(I,2);nz = size(I,3);
nradius = max(max(nx,ny),nz)/2;
ntheta = 2*nradius;
nphi = 2*ntheta;
[r phi theta] = ndgrid(0:nradius-1,0:2*pi/nphi:2*pi-2*pi/nphi,0:pi/(ntheta-1):pi);
% polar coordinates in cartesian space
%eps = 10^(-12);%added due to numerical trouble with floor
eps = 0;
px = r.*cos(phi).*sin(theta)+nradius+1+eps;
py = r.*sin(phi).*sin(theta)+nradius+1+eps;
pz = r.*cos(theta)+nradius+1;
clear r phi theta; 
%%%%%%%%%%%%%% trilinear interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate levers
tx = px-floor(px);
ty = py-floor(py);
tz = pz-floor(pz);
%perform interpolation
pol = (1-tx).*(1-ty).*(1-tz).*I(floor(px)+nx*(floor(py)-1)+ny*nx*(floor(pz)-1)) + ...
    (tx).*(1-ty).*(1-tz).*I(ceil(px)+nx*(floor(py)-1)+ny*nx*(floor(pz)-1)) + ...
    (1-tx).*(ty).*(1-tz).*I(floor(px)+nx*(ceil(py)-1)+ny*nx*(floor(pz)-1)) + ...
    (1-tx).*(1-ty).*(tz).*I(floor(px)+nx*(floor(py)-1)+ny*nx*(ceil(pz)-1)) + ...
    (tx).*(ty).*(1-tz).*I(ceil(px)+nx*(ceil(py)-1)+ny*nx*(floor(pz)-1)) + ...
    (tx).*(1-ty).*(tz).*I(ceil(px)+nx*(floor(py)-1)+ny*nx*(ceil(pz)-1)) + ...
    (1-tx).*(ty).*(tz).*I(floor(px)+nx*(ceil(py)-1)+ny*nx*(ceil(pz)-1)) + ...
    (tx).*(ty).*(tz).*I(ceil(px)+nx*(ceil(py)-1)+ny*nx*(ceil(pz)-1));