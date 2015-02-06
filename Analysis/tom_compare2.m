function out = tom_compare(vol1, vol2, num)
% TOM_COMPARE comparison of 2 fourier transforms
%
%   out = tom_compare(vol1, vol2, num)
%
%PARAMETERS
%
%  INPUT
%   vol1        VOLUME 1 - in real space!
%   vol2        VOLUME 2 - in real space!
%   num         Number of shells in reciprocal space
%
%  OUTPUT       Array of dimension 10xNUM
%
%EXAMPLE
%   results = tom_compare(vol1, vol2, num)
%
%   The Fourier space is subdivided into num rings  of  spherical-shell
%   zones.  For  each  zone  values  are  calculated  which  serve  as
%   measures for the agreement of Fourier coefficients, e.g., the mean
%   square deviation.
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by M. Riedlberger 10/18/03
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
%   http://www.biochem.mpg.de/tom

error(nargchk(0, 3, nargin, 'struct'))

if (ndims(vol1) ~= ndims(vol2));
    error('Number of dimensions of input volumes do not equal!');
    return;int
end;
if (size(vol1) ~= size(vol2));
    error('Dimensions of input volumes do not equal!');
end;
if (size(vol1,1) ~= size(vol1,2) || size(vol1,3) ~= size(vol1,2));
    error('x, y and z dimensions must be equal!');
end;


s = size(vol1,1);
f1 = fftshift(fftn(double(vol1)));
f2 = fftshift(fftn(double(vol2)));
f12 = f1-f2;
f1 = single((f1.*conj(f1))/s^3);
f2 = single((f2.*conj(f2))/s^3);
f12 = single((f12.*conj(f12))/s^3);
out = single(zeros(num,11));
tom_comparec2(f1, f2, f12, num, 1, num, out);


