function tom_euler2matrix(ang,not)
%TOM_EULERCONVERT_XMIPP converts euler angles from xmipp to the TOM
%convention
%
%   [rotmatrix,euler_out] = tom_eulerconvert_xmipp(rot, tilt, psi)
%   
%
%PARAMETERS
%
%  INPUT
%   rot                 input xmipp angle rot
%   tilt                input xmipp angle tilt
%   psi                 input xmipp angle psi
%  
%  OUTPUT
%
%   rotmatrix           resulting rotation matrix
%   euler_out           resulting Euler angles
%
%
%EXAMPLE
%   [rotmatrix,euler_out] = tom_eulerconvert_xmipp(10,20,30)
%
%   For more information about 3d Euler rotation visit:
%   wiki -> Definitions -> Overview 3d Euler rotation -> 
%   Rotation matrices, zxz and zyz conventions
% 
%
%REFERENCES
%
%SEE ALSO
%   tom_sum_rotation
%
%   created by FB 04/05/07
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

if (strcmp(not,'zyz'))
    
    
end;

disp('NOT Implemented!')