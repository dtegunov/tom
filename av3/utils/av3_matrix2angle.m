function [phi, psi, theta] = av3_matrix2angle(matr)
% AV3_MATRIX2ANGLE converts rotation matrix to angles
%
%   [phi psi theta] = av3_matrix2angle(matr)
%
%   routine determines Euler angles corresponding to 
%   a given 3x3 rotation matrix.
%
zz = matr(3,3);
theta = acos(zz);
yz = matr(3,2);
zy = matr(2,3);
if sin(theta) ~=0
    cosphi = yz/sin(theta);
    phi = acos(cosphi);
    cospsi = -zy/sin(theta);
    psi = acos(cospsi);
else
    disp('theta = 0 remains to be implemented ...');
end;
psi = psi*180/pi;phi = phi*180/pi;theta = theta*180/pi;
