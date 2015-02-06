function tri=tom_triangle(sz)

%TOM_TRIANGLE creates an image containing a triangle
%
%    vol=tom_triangle(sz)
%
%PARAMETERS
%
%  INPUT
%   sz              image size
%  
%  OUTPUT
%   tri             triangle image
%
%EXAMPLE
%   triangle=tom_triangle([128]);
%   tom_imagesc(triangle)
%
%REFERENCES
%
%SEE ALSO
%   TOM_SPHERE, TOM_CIRCLE, TOM_TORUS
%
%   created by SN 10/29/09
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

tri=zeros(sz,sz);
ix=0;
for iy=1:1:sz
        ix=ix+1;
        tri(ix:end-ix+1,iy)=1;
end;
tri=tri(:,1:round(sz/2));

