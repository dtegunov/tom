function cyl=tom_cylinder(radius_in, radius_out, overall_size)
%TOM_CYLINDER creates an 3D array with a cylinder inside
%
%   cyl=tom_cylinder(radius_in, radius_out, overall_size)
%
%PARAMETERS
%
%  INPUT
%   radius_in           a scalar defining the inner radius
%   radius_out          a scalar defining the outer radius
%   overall_size        size of 3D array
%  
%  OUTPUT
%   cyl                 a 3D array with the cylinder
%
%EXAMPLE
%   cyl=tom_cylinder(8, 10, [32 32 32]);
%   creates a 32 cube cyl with a cylinder of an inner
%   radius of 8 and an outer radius of 10.
%
%REFERENCES
%
%SEE ALSO
%   TOM_CIRCLE
%
%   created by SN 04/04/04
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

% check input arguments
error(nargchk(0, 3, nargin, 'struct'))

cyl=zeros(overall_size);
cyli=zeros(overall_size);
c=tom_circle(radius_out);
for lauf=1:overall_size(3)
    cyl=tom_paste(cyl,c,[overall_size(1)./2-radius_out+1 overall_size(2)./2-radius_out+1 lauf]);
end;
c=-tom_circle(radius_in);
for lauf=1:overall_size(3)
    cyli=tom_paste(cyli,c,[overall_size(1)./2-radius_in+1 overall_size(2)./2-radius_in+1 lauf]);
end;
cyl=cyl+cyli;
