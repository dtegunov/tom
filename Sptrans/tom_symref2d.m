function  out=tom_symref2d(in,factor)
%TOM_SYMREF2D creates ...
%
%   out=tom_symref2d(in,factor)
%
%PARAMETERS
%
%  INPUT
%   in                  ...
%   factor              ...
%  
%  OUTPUT
%   out                 ...
%
%EXAMPLE
%   ... = tom_symref2d(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   last change 
%   10/20/10 FF - corrected normalization bug
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

out=ones(size(in)).*mean2(in);
for angle=0:(360./factor):359    
    out=(imrotate(in,angle,'bilinear','crop')+out);
end
out = out ./ factor;

