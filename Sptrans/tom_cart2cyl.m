function pol = tom_cart2cyl(I)
%TOM_CART2CYL transforms 3D-volumes from cartesian to cylinder coordinates
%
%   pol = tom_cart2cyl(I)
%
%   A 3D image I given in cartesian coordinates is sampled in cylinder coordinates 
%   R, PHI, Z using trilinear interpolation. The input image I is assumed to be of 
%   equal dimensions NX and NY. The ouput POL(R,PHI,Z) has dimensions NR=NX/2, 
%   NPHI=4*NR, NZ=NZ. Function works only for 3D, for 2D data use
%   TOM_CART2POLAR, which transforms to POLAR coordinates
%
%PARAMETERS
%
%  INPUT
%   I                   3 dim array
%  
%  OUTPUT
%   pol                 3 dim array in cylinder coordinates
%
%EXAMPLE
%   cart = zeros(32,32,32);
%   cart(8,8,:) = 1;
%   cart = tom_symref(cart,4);
%   cyl = tom_cart2cyl(cart);
%   tom_dspcub(cyl);
%
%REFERENCES
%
%SEE ALSO
%   TOM_POLAR2CART, TOM_CART2SPH, TOM_SPH2CART
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
%if nz > 1
%    error(' use for tom_cart2sph for 3D arrays!')
%end;
nradius = max(nx,ny)/2;
nphi = 4*nradius;
[r phi pz] = ndgrid(0:nradius-1,0:2*pi/nphi:2*pi-2*pi/nphi,1:nz);
% polar coordinates in cartesian space
%eps = 10^(-12);%added due to numerical trouble with floor
eps = 0;
px = r.*cos(phi)+nradius+1+eps;
py = r.*sin(phi)+nradius+1+eps;
clear r phi ; 
%%%%%%%%%%%%%% bilinear interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate levers
tx = px-floor(px);
ty = py-floor(py);
%perform interpolation
%   check for undefined indexes
%mxx = max(max(floor(px+1)));
%if  mxx > nx
%    I(mxx,ny) =0;
%    nx = mxx;
%end;
%mxy = max(max(floor(py+1)));
%if mxy > ny
%    I(nx,mxy)=0;
%    ny = mxy;
%end;    
%pol = (1-tx).*(1-ty).*I(floor(px)+nx*(floor(py)-1)) + ...
%    (tx).*(1-ty).*I(floor(px+1)+nx*(floor(py)-1)) + ...
%    (1-tx).*(ty).*I(floor(px)+nx*(floor(py+1)-1)) + ...
%    (tx).*(ty).*I(floor(px+1)+nx*(floor(py+1)-1));
pol = (1-tx).*(1-ty).*I(floor(px)+nx*(floor(py)-1)+nx*ny*(pz-1)) + ...
    (tx).*(1-ty).*I(ceil(px)+nx*(floor(py)-1)+nx*ny*(pz-1)) + ...
    (1-tx).*(ty).*I(floor(px)+nx*(ceil(py)-1)+nx*ny*(pz-1)) + ...
    (tx).*(ty).*I(ceil(px)+nx*(ceil(py)-1)+nx*ny*(pz-1));
