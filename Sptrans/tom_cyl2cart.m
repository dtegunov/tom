function I = tom_cyl2cart(pol)
%TOM_CYL2CART transforms 3D-images from polar to cartesian coordinates
%
%   I = tom_cyl2cart(pol)
%
%   A 3D-image I is sampled in cylinder coordinates RHO, PHI, Z. TOM_CYL2CART 
%   samples the data on cartesian coordinates using bilinear interpolation.
%   The dimensions of the polar data POL are the radius R, the polar angle
%   PHI, and Z. The dimensions are assumed to be (NDIMR, 4*NDIMR, NZ). The 
%   dimensions of I are (2*NDIMR,2*NDIMR,NZ)
%
%PARAMETERS
%
%  INPUT
%   pol                 3 dim array in cylinder coordinates
%  
%  OUTPUT
%   I                   3 dim array
%
%EXAMPLE
%   cyl = zeros(16,64,32);
%   cyl(8,:,:) = 1;
%   cart = tom_cyl2cart(cyl);
%   tom_dspcub(cart);
%
%REFERENCES
%
%SEE ALSO
%   TOM_CART2POLAR, TOM_CART2SPH, TOM_SPH2CART, TOM_CART2CYL, TOM_POLAR2CART
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


nphi = size(pol,2);nradius = size(pol,1);nz=size(pol,3);
nx = nphi/2;ny = nphi/2;
%if ntheta > 1
%    error(' use for tom_cart2sph for 3D arrays!')
%end;
[x y z] = ndgrid(-nradius:nradius-1,-nradius:nradius-1,1:nz);
% cartesian coordinates in polar space
%eps = 10^(-12);%added due to numerical trouble with floor
eps = 0;
cr = sqrt(x.^2+y.^2)+1;
cphi = atan2(y,x);
indx = find( cphi < 0);
cphi(indx) = 2*pi-abs(cphi(indx));
cphi=cphi/(2*pi)*nphi+1;
clear x y; 
%360 deg == 0 deg
nphi = nphi+1;
pol(:,nphi,:) = pol(:,1,:);
%extend in r dimension ... cheap solution...
nradius = ceil(sqrt(2)*nradius)+1;
pol(nradius,nphi,nz) = 0;

%%%%%%%%%%%%%% bilinear interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate levers
tr = cr-floor(cr);
tphi = cphi-floor(cphi);
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
I = (1-tr).*(1-tphi).*pol(floor(cr)+nradius*(floor(cphi)-1)+nradius*nphi*(z-1)) + ...
    (tr).*(1-tphi).*pol(ceil(cr)+nradius*(floor(cphi)-1)+nradius*nphi*(z-1)) + ...
    (1-tr).*(tphi).*pol(floor(cr)+nradius*(ceil(cphi)-1)+nradius*nphi*(z-1)) + ...
    (tr).*(tphi).*pol(ceil(cr)+nradius*(ceil(cphi)-1)+nradius*nphi*(z-1));