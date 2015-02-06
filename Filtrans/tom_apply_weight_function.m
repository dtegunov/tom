function out=tom_apply_weight_function(in,w_func,mask,ang)
%TOM_APPLY_WEIGHT_FUNCTION applies weighting function to an image or volume
%
%   out=tom_apply_weight_function(in,w_func,mask,ang)
%
%PARAMETERS
%
%  INPUT
%   in                  image or volume
%   w_func              weighting function
%   mask                mask for imag or volume
%   ang                 angle for rotation
%  
%  OUTPUT
%   out                 weighted function
%                       rotation angles
%
%EXAMPLE
%  sphere=tom_spheremask(ones(64,64,64),18);
%  wedge=tom_wedge(ones(64,64,64),18);
%  wedged_sphere=tom_apply_weight_function(sphere,wedge);
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB 2005
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

if (nargin==2)
    mask=ones(size(in));
    ang=0;
end;

if (nargin==3)
   ang=0;
end;


if (ang~=0)
    in=single(tom_rotate(single(in),ang,'linear'));
end;

%check fourier space
fin = fftshift(fftn(double(in)));

%apply weighting
out=double(fin).*double(w_func).*double(mask);

%back to real space
out = real(ifftn(ifftshift(out)));

%rotate back
if (ang~=0)
    out=single(tom_rotate(single(out),-ang,'linear'));
end;

