function im = tom_shift_fixed(im, delta)
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

indx = find([dimx,dimy,dimz] == 1);
delta(indx)=0;

x_shift = exp(i * 2 * pi * delta(1) * [0:ceil(dimx/2-1) ceil(-dimx/2):-1] / dimx);
y_shift = exp(-i * 2 * pi * delta(2) * [0:ceil(dimy/2-1) ceil(-dimy/2):-1] / dimy);
z_shift = exp(-i * 2 * pi * delta(3) * [0:ceil(dimz/2-1) ceil(-dimz/2):-1] / dimz);

if mod(dimx, 2) == 0
	x_shift(dimx/2+1) = real(x_shift(dimx/2+1));
end 
if mod(dimy, 2) == 0
	y_shift(dimy/2+1) = real(y_shift(dimy/2+1));
end
if mod(dimz, 2) == 0
	z_shift(dimz/2+1) = real(z_shift(dimz/2+1));
end


im = fftn(im);
shift = reshape(reshape(x_shift' * y_shift, dimx*dimy, 1) * z_shift, dimx, dimy, dimz);
fim = im .* shift;


im = ifftn(fim);