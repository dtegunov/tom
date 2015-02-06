function im = tom_shift(im, delta)
%TOM_SHIFT shifts image by a vector
%
%   im = tom_shift(im, delta)
%
%PARAMETERS
%
%  INPUT
%   im                  1D, 2D or 3D array
%   delta               1D, 2D or 3D vector for shift
%  
%  OUTPUT
%   im                  shifted 1D, 2D or 3D array
%
%   The shift is performed in Fourier space, thus periodic 
%   boundaries are implicit. The shift vectors do not have to be
%   integer.
%   currently restricted to even dimensions of IMAGE
%
%EXAMPLE
%   yyy = tom_sphere([64 64 64],10,1,[16 16 16]);
%   yyy = tom_shift(yyy, [1,2,3]);
%   tom_dspcub(yyy);
%
%REFERENCES
%
%SEE ALSO
%   TOM_MOVE, TOM_SHIFT_FFT
%
%   created by FF 08/01/02
%   updated by FF 01/16/04
%   profiled and optimized by AK 09/08/07
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

if size(delta,1) > 1
    delta = delta';
end

[dimx,dimy,dimz]=size(im);
%c1 = double(int32(dimx/2));
%c2 = double(int32(dimy/2));
%MeshGrid with the sampling points of the image
[x,y,z]=ndgrid( -floor(dimx/2):-floor(dimx/2)+(dimx-1),...
    -floor(dimy/2):-floor(dimy/2)+dimy-1, ...
    -floor(dimz/2):-floor(dimz/2)+dimz-1);
indx = find([dimx,dimy,dimz] == 1);
delta(indx)=0;
delta = delta./[dimx dimy dimz];
x = delta(1)*x + delta(2)*y + delta(3)*z; clear y; clear z;

im = (fftn(im));
fim = im.*exp(-2*pi*i*ifftshift(x));
im = real(ifftn(fim));
