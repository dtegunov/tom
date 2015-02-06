function out= tom_binc(im, binning)
%TOM_BINC performs binning of 1D, 2D or 3D images
%
%   out= tom_binc(im, binning)
%
%PARAMETERS
%
%  INPUT
%   im                  ...
%   binning             ...
%  
%  OUTPUT
%   out                 ...
%
%EXAMPLE
%   E.g. 1D: pixels 1 and 2 are averaged and stored as pixel one of BINNED, pixels 3 and 4 
%   of IM are averaged and stored as pixel 2 of BINNED. The dimensions of
%   BINNED are half of IM. For multiple binning TOM_BINC nbin >1 can be
%   used
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 03/05/06
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


if nargin == 1
    binning = 1;
end

binning = round(binning);

if nargin > 2
    error('Two input arguments required.');
end

if size(im,1) > 1 && size(im,2) > 1
    if binning > 0
       im = single(im);
       binning = 2^binning;
       out = tom_bininc(im, [binning binning binning]);
    else
        out = im;
    end
else
    error('cannot bin 1D data, use tom_bin instead!');
end