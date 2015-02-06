function [proj] = tom_vol2proj(vol,angles)
%TOM_VOL2PROJ calculates a 2D projection of a 3D volume. C implementation.
%
%   proj=tom_vol2proj(vol,angles)
%
%PARAMETERS
%
%  INPUT
%   vol            3D volume (single or double, in C code handled as float), 
%   angeles    a 2-vector [phi theta] (single or double,  in C code handled as double) 
%                    in degrees, which is describing the orientation of the projection in 3d-space.
%                    For example, phi is the projection direction and theta the tiltangle.
%  
%  OUTPUT
%   proj          2D projection (single, in C code handled as float)
%
%EXAMPLE
%   cyl = tom_cylinder(8, 10, [32 32 32]); % Creates volume
%   [cylproj] = tom_vol2proj(cyl,[45 45]); % Calculates projection 
%
%REFERENCES
%   EM package, R. Hegerl
%
%SEE ALSO
%   TOM_BACKPROJ3D
%
%   created by ME 01/02/08
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

% C implementation
[proj] = tom_vol2projc(single(vol),double(angles));
