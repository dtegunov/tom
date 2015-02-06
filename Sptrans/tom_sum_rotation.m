function  [euler_out shift_out rott]=tom_sum_rotation(rots,shifts,order)
%TOM_SUM_ROTATION sums up N translations and N 3-tupel of Euler angles ...
%
%   [euler_out shift_out rott]=tom_sum_rotation(rots,shifts,order)
%
%   TOM_SUM_ROTATION sums up N translations and N 3-tupel of Euler angles
%   to only one rotation and translation
%
%PARAMETERS
%
%  INPUT
%   rots                N x 3 matrix with Euler angles, like this:
%                        rots = [phi[1] psi[1] theta[1]; ... phi[N] psi[N]
%                        theta[N]];
%   shifts              N x 3 matrix with translations, like this:
%                        shifts = [x[1] y[1] z[1]; ... x[N] y[N] z[N]];
%   order               Either translate - rotate (trans_rot) or rotate -
%                       translate (rot_trans)
%
%  OUTPUT
%   euler_out           resulting Euler angles
%   shift_out           resulting translation vector
%   rott                resulting rotation matrix
%
%EXAMPLE
%   [euler_out shift_out rott]=tom_sum_rotation([10 20 30; -20 -10 -30],[5 5 5; 5 5 5]);
% 
%   [euler_out shift_out rott]=tom_sum_rotation([10 20 30; -20 -10 -30],[5 5 5; -2.9506 -3.7517  -7.2263]);
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

if nargin==2
    order='trans_rot';
end;



if strcmp(order,'rot_trans')
    
    %build up zxz rotation Matrix
    
    for i=1:size(rots,1)
        
        phi=rots(i,1).*pi./180;
        psi=rots(i,2).*pi./180;
        theta=rots(i,3).*pi./180;
        
        % use zxz matrix
        rotM(:,:,i)=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1]*...
            [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)]*...
            [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
        
    end;
    
    
    %sum up the shifts
    shift_out=zeros(size(shifts(1,:)))';
    for i=1:size(rots,1)
        z=size(rots,1);
        rott=eye(size(rots,2));
        for ii=(size(rots,1)):-1:i+1
            rott=rott*rotM(:,:,z);
            z=z-1;
        end;
        shift_out=rott*shifts(i,:)'+shift_out;
    end;
    
    %sum up the rotations
    rott=eye(size(rots,2));
    z=size(rots,1);
    for ii=1:size(rots,1)
        rott=rott*rotM(:,:,z);
        z=z-1;
    end;
    
    
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
    
else
%%
%   order == trans_rot
    %build up zxz rotation Matrix
    
    for i=1:size(rots,1)
        
        phi=rots(i,1).*pi./180;
        psi=rots(i,2).*pi./180;
        theta=rots(i,3).*pi./180;
        
        % use zxz matrix
        rotM(:,:,i)=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1]*...
            [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)]*...
            [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
        
    end;
    
    
    %sum up the shifts
    shift_out=zeros(size(shifts(1,:)))';
    for i=1:size(rots,1)
        z=size(rots,1);
        rott=eye(size(rots,2));
        for ii=1:(size(rots,1)-(i-1))
            rott=rott*rotM(:,:,z);
            z=z-1;
        end;
        shift_out=rott*shifts(i,:)'+shift_out;
    end;
    
    %sum up the rotations
    rott=eye(size(rots,2));
    z=size(rots,1);
    for ii=1:size(rots,1)
        rott=rott*rotM(:,:,z);
        z=z-1;
    end;
    
    
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
    
end;
