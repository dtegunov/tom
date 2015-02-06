function powers = tom_ps(im)
%TOM_PS calculates and displays the power spectrum.
%
%   powers = tom_ps(im)
%
%PARAMETERS
%
%  INPUT
%   im                  input data (2D or 3D)
%  
%  OUTPUT
%   powers              power spectrum
%
%   Calculates the powerspectrum (squared amplitudes of the 
%   Fourier transform) of a multidimansional (2D or 3D) image
%   Zero frequency located in the middle!
%
%EXAMPLE
%   im = tom_emread('proteasome.em');
%   ps = tom_ps(im.Value);
%   tom_imagesc(log(ps));
%
%REFERENCES
%
%SEE ALSO
%   TOM_FOURIER, TOM_IFOURIER, FFTSHIFT, TOM_CTFFIT
%
%   created by FF
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

[s1,s2,s3]=size(im);

if isequal(s3,1)
    powers = abs(fftshift(fft2(im)));
    powers = powers.^2;
elseif s3>1
    powers = abs(fftshift(fftn(im)));
    powers = powers.^2;
    
else
    error('Dimensions do not match');
end
