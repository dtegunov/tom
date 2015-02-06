function fim = tom_shift_fft(fim, delta)
%TOM_SHIFT_FFT shifts an image by a vector
%
%   fim = tom_shift_fft(fim, delta)
%
%   An array FIM is shifted in Fourier space - the array FIM has to be
%   complex (~ in Fourier space. currently restricted to even dimensions of
%   IMAGE
%PARAMETERS
%
%  INPUT
%   fim                 1D, 2D array in Fourier space
%   delta               1D, 2D vector for shift
%  
%  OUTPUT
%   fim                 shifted 1D, 2D array in Fourier space
%
%EXAMPLE
%   xxx= zeros(128,128);xxx(1,1,1)=1;
%   tom_imagesc(xxx);
%   xxx=tom_fourier(xxx);xxx=tom_shift_fft(xxx,[1.5,1.5]);
%   xxx=real(tom_ifourier(xxx));figure;imagesc(xxx')
%
%REFERENCES
%
%SEE ALSO
%   TOM_MOVE, TOM_SHIFT
%
%   created by FF 02/19/04
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


[dimx,dimy,dimz]=size(fim);
%c1 = double(int32(dimx/2));
%c2 = double(int32(dimy/2));
%MeshGrid with the sampling points of the image
[x,y,z]=ndgrid( -floor(size(fim,1)/2):-floor(size(fim,1)/2)+(size(fim,1)-1),...
    -floor(size(fim,2)/2):-floor(size(fim,2)/2)+size(fim,2)-1,...
    -floor(size(fim,3)/2):-floor(size(fim,3)/2)+(size(fim,3)-1));
indx = find([dimx,dimy,dimz] == 1);
delta(indx)=0;
delta = delta./[dimx dimy dimz];
x = delta(1)*x + delta(2)*y + delta(3)*z; clear y; clear z;
fim = fftshift(fim);
fim = ifftshift(fim.*exp(-2*pi*i*x));