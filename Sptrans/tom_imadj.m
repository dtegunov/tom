function I_out=tom_imadj(I,fmin,fmax)
%TOM_IMADJ  scales image between 0 and 1, range can be constrained by fmin and fmax
%
%   I_out=tom_imadj(I,fmin,fmax)
%
%PARAMETERS
%
%  INPUT
%   I                   ...
%   fmin                ...
%   fmax                ...
%  
%  OUTPUT
%   I_out               ...
%
%EXAMPLE
%   ... = tom_imadj(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FF 07/26/02 changed function of AF
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


error(nargchk(1,3,nargin))
%Use of the hole dynamic range per default
fmaxt = max(max(max(I)));
fmint = min(min(min(I)));

if (nargin < 3)
    fmax = fmaxt;
end;
if (nargin < 2)
    fmin = fmint;
end;
if (fmin < fmint)
    disp(['given minimum (' num2str(fmin) ') smaller than minimum value in image (' num2str(fmint) ')' ])
    fmin = fmint;
end;
if (fmax > fmaxt)
    disp(['given maximum (' num2str(fmax) ') bigger than maximum value in image (' num2str(fmaxt) ')' ])
    fmax = fmaxt;
end;

%Justify Image between 0 and 1
I=tom_limit(I,fmin, fmax);
I_out = (I-fmin)/((fmax-fmin)+0.00000001);

