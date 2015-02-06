function  [euler_out shift_out rott]=tom_sum_rotation_zyz(rots)
%TOM_SUM_ROTATION_ZYZ sums up N 3-tupel of Euler angles in zyz
%
%   [euler_out shift_out rott]=tom_sum_rotation_zyz(rots)
%
%   TOM_SUM_ROTATION sums up N translations and N 3-tupel of Euler angles
%   to only one rotation and translation
%
%PARAMETERS
%
%  INPUT
%   rots                N x 3 matrix with Euler angles, like this:
%                        rots = [phi[1] psi[1] theta[1]; ... phi[N] psi[N]
%
%  OUTPUT
%   euler_out           resulting Euler angles
%   rott                resulting rotation matrix
%
%EXAMPLE
%   [euler_out shift_out rott]=tom_sum_rotation([10 20 30; -20 -10 -30],[5 5 5; 5 5 5]);
% 
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
%   created by ... (author date)
%   updated by ...
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


for i=1:size(rots,1)
    
    phi=rots(i,1).*pi./180;
    psi=rots(i,2).*pi./180;
    theta=rots(i,3).*pi./180;
    
    % use zyz matrix
%     rotM(:,:,i)=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1]*...
%         [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]*...
%         [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    

   rotM(:,:,i)=[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]*...
        [cos(psi) 0 sin(psi); 0 1 0; -sin(psi) 0 cos(psi)]*...
        [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];


end;

rott=eye(size(rots,2));
z=size(rots,1);
for ii=1:size(rots,1)
    rott=rott*rotM(:,:,ii);
    z=z-1;
end;

%extract euler angles
% euler_out(1)=atan2(rott(2,3),rott(1,3));
% euler_out(2)=atan2(rott(2,3),-rott(1,3));
% euler_out(3)=acos(rott(3,3));
% euler_out=euler_out.*180./pi;



euler_out(1)=atan2(rott(2,3),rott(1,3));
euler_out(2)=acos(rott(3,3));
euler_out(3)=atan2(rott(3,2),-rott(3,1));

euler_out=euler_out.*180./pi;

if -(rott(3,3)-1)<10e-8
    euler_out(1)=0;
    euler_out(2)=0;
    euler_out(3)=atan2(rott(2,1),rott(1,1)).*180./pi;
end



    
    
    