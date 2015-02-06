function [angle_out] = av3_index2angle(index, angle_start, angular_incr,angle_end)
% AV3_INDEX2ANGLE creates Euler angles from molmatch indices
%
%   [angle_out] = av3_index2angle(index, angle_start,
%   angular_incr,angle_end)
%
%   if phi_start, psi_start, and theta_start are not given they are set to zero.
%   angles mod(angular_incr) has to be zero!
%
%   routine determines Euler angles phi, psi, and theta from a given
%   index. index is assumed to be generated by loop, from inner to outer: phi, psi, theta,
%   lowest index is zero.
%   
% PARAMETERS
%  INPUT
%   index        loop number to be converted into Euler angles
%   angle_end   [phi_end psi_end theta_end] phi psi and theta choosen as 
%                         start value as vector
%   angular_incr  [phi_increment psi_increment theta_increment] as vector
%   angle_start     [phi_start psi_start theta_start]      angle chosen as
%                            start value as vector 
%
%  OUTPUT
%   phi          [phi psi theta] angle corresponding to INDEX as vector
%
%   
%  SEE ALSO
%   molmatch.exe (MPI c-program), omnimatch.exe, AV3_CREAMOTL
%
%   08/09/02 FF
%  last change 03/31/05 FF - docu updated


nangle=(angle_end-angle_start)./(angular_incr)+1;

angle_out(3) = floor(index /(nangle(1).*nangle(2)));
rest = index - angle_out(3).*(nangle(1).*nangle(2));
angle_out(2) = floor(rest/nangle(1));
angle_out(1) = rest - angle_out(2).*nangle(1);

angle_out = angle_out .*angular_incr + angle_start; 





