function out=rot3d(in,phi,psi,theta,ip);
%  rot3d performs 3d rotation
%
%  Syntax: out = rot3d(in,phi,psi,theta,ip);
%          out = 3D image output volume
%          in = 3D image input volume
%          phi, psi, theta = Euler angles 
%	   ip = interpolation type: 'linear', 'splines'

%  Last changes: Oct. 21, 2003
%  M. Riedlberger


out = single(zeros(size(in)));
rot3dc (single(in),out,phi,psi,theta,ip);
out = double(out);
