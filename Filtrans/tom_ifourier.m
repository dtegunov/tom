function c = tom_ifourier(a)
%TOM_IFOURIER calculates the inverse Fourier Transformation of Matrix A
%
%   c = tom_ifourier(a)
%
%PARAMETERS
%
%  INPUT
%   a                   ...
%  
%  OUTPUT
%   c           		...
%
%EXAMPLE
%   ... = tom_ifourier(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_FOURIER, TOM_FILTER
%
%   created by AL 09/11/02
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
    c=ifft(a);
elseif s3 == 1
    c=ifft2(a);
else
    c=ifftn(a);
end

