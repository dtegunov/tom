function [PROJ] = tom_proj3d2c(VOL,ANGLES)
%TOM_PROJ3D2C calculates a 2D projection of a 3D volume. C implementation.
%
%   proj=tom_proj3d2c(vol,angles)
%
%PARAMETERS
%
%  INPUT
%   VOL          3D volume
%   ANGLES    a 2-vector [PHI THETA] in degrees, which is describing the 
%                    orientation of the projection in 3d-space. 
%                    For example, PHI is the projection direction and THETA the tiltangle.
%  
%  OUTPUT
%   PROJ         2D projection
%
%EXAMPLE
%   cyl = tom_cylinder(8, 10, [32 32 32]); % Creates volume
%   [cylproj] = tom_proj3d2c(cyl,[45 45]); % Calculates projection 
%
%REFERENCES
%   EM package, R. Hegerl
%
%SEE ALSO
%   TOM_BACKPROJ3D
%
%   created by ME 13/12/06
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

error(['No compiled mex-function found with name ' mfilename()]);
