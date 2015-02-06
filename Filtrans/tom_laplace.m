function image = tom_laplace(image)
%TOM_LAPLACE performs Laplace filtering of image or volume
%
%   image = tom_laplace(image)
%
%PARAMETERS
%
%  INPUT
%   image               image or volume to be filtered
%  
%  OUTPUT
%   image       		filtered image or volume
%
%   Laplacian filtering is essentially the 2nd derivative of the input
%   image. The filtering is performed in Fourierspace, i.e. the Fourier
%   transformed of the image is multipied by |k|^2. Works for 1D, 2D, 3D
%   data.
%
%EXAMPLE
%   yyy = tom_sphere([64 64 64],10,1,[16 16 16]);
%   yyy = tom_laplace(yyy);
%   tom_dspcub(yyy);
%
%REFERENCES
%
%SEE ALSO
%   TOM_FILTER TOM_BANDPASS
%
%   created by FF 04/17/03
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


[x,y,z]=ndgrid( -floor(size(image,1)/2):-floor(size(image,1)/2)+(size(image,1)-1),...
    -floor(size(image,2)/2):-floor(size(image,2)/2)+size(image,2)-1, ...
    -floor(size(image,3)/2):-floor(size(image,3)/2)+size(image,3)-1);
x = x.^2+y.^2+z.^2;
image = fftshift(tom_fourier(image));
image = real(tom_ifourier(ifftshift(x.*image)));
