function c = tom_fourier(a)
%TOM_FOURIER calculates the Fourier Transformation of Matrix A
%
%   c = tom_fourier(a)
%
%PARAMETERS
%
%  INPUT
%   a                   Matrix
%  
%  OUTPUT
%   c           		Matrix
%
%   C=TOM_FOURIER(A) This function calculates the Fourier Transformation of the input matrix A.
%   Matrix A can be a 1-D, 2-D or 3-D matrix. the result is stored in the
%   matrix C.
%
%EXAMPLE
%          11    3   -6  -14             16      3+19i   -6      3+19i
%       a = 1  -13    9    6        c = -24-2i   50+9i    8-64i  1+21i
%           0    8   13   -3             8       5-37i    54     5+37i
%          -8   17  -15    7            -24+2i   1-21i    8+64i  59-9i
%
%REFERENCES
%
%SEE ALSO
%   TOM_FILTER, TOM_ORCD, TOM_RTSYM
%
%   created by AL 08/07/02
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


[s1 s2 s3]=size(a);

if s1 == 1 | s2 == 1
    c=fft(a);
elseif s3 == 1
    c=fft2(single(a));
else
    c=fftn(a);
end

