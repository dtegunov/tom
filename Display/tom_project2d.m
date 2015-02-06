function m=tom_project2d(volume,axis)
%TOM_PROJECT2D is a function thaat projects the 3-D Volume in a axis-oriented projection. In
% other words one can aquire the sum of the volume in the xy,xz or yz
% direction.
%
%   m=tom_project2d(volume,axis)
%
%PARAMETERS
%
%  INPUT
%   volume              ...
%   axis                ...
%  
%  OUTPUT
%   m           		...
%
%EXAMPLE
%   x=randn(3,3,3);
%   prj=tom_project2d(x,'xy');
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AL 11/14/05
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

error(nargchk(0, 2, nargin, 'struct'))

[dim1 dim2 dim3]=size(volume);  
switch (axis)
 
case 'xy'

    m=sum(volume,3);
case 'yz'
    volume=permute(volume,[3 2 1]);
    m=sum(volume,3);
case 'xz'
    volume=permute(volume,[3 1 2]);
    m=sum(volume,3);
otherwise
       error('Unknown orientation'); 
end


