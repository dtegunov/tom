function [euler_out] = tom_rotmatrix2angles(rott)
%TOM_ROTMATRIX2ANGLES converts a rotation matrix to the corresponding ...
%
%   [euler_out] = tom_rotmatrix2angles(rott)
%
%   TOM_ROTMATRIX2ANGLES converts a rotation matrix to the corresponding
%   euler angles. The rotation matrix has to be given in the zxz form
%
%PARAMETERS
%
%  INPUT
%   rott                 3x3 zxz rotation matrix
%  
%  OUTPUT
%   euler_out            resulting Euler angles
%
%EXAMPLE
%   [euler_out]=tom_sum_rotmatrix2angles([])
%
%   For more information about 3d Euler rotation visit:
%   wiki -> Definitions -> Overview 3d Euler rotation -> 
%   Rotation matrices, zxz and zyz conventions
% 
%   For more information about sequential translation and Euler rotation visit:
%   wiki -> Definitions -> Overview 3d Euler rotation -> 
%   Sequential translation and rotation
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

%extract euler angles
euler_out(1)=atan2(rott(3,1),rott(3,2));
euler_out(2)=atan2(rott(1,3),-rott(2,3));
euler_out(3)=acos(rott(3,3));
euler_out=euler_out.*180./pi;

if -(rott(3,3)-1)<10e-8
    euler_out(3)=0;
    euler_out(2)=0;
    euler_out(1)=atan2(rott(2,1),rott(1,1)).*180./pi;
end