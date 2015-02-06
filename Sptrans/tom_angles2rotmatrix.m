function rotmatrix = tom_angles2rotmatrix(euler_in)
%TOM_ANGLES2ROTMATRIX converts euler angles to the corresponding rotation matrix
%
%   rotmatrix = tom_angles2rotmatrix(euler_in)
%
% For more information about 3d Euler rotation visit:
% wiki -> Definitions -> Overview 3d Euler rotation -> 
% Rotation matrices, zxz and zyz conventions
% 
% For more information about sequential translation and Euler rotation visit:
% wiki -> Definitions -> Overview 3d Euler rotation -> 
% Sequential translation and rotation
%
%PARAMETERS
%
%  INPUT
%   euler_in            3x1 vector with Euler angles (phi psi theta)s
%  
%  OUTPUT
%   rotmatrix           3x3 zxz rotation matrix
%
%EXAMPLE
%   [rotmatrix]=tom_angles2rotmatrix([30 60 90]);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 10/24/05
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

[e s rotmatrix] = tom_sum_rotation(euler_in,[0 0 0]);