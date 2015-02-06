function motl = normalvec(motl, cent)
% NORMALVEC determines Euler angles of particles on a sphere
%
%   motl = normalvec(motl, cent)
%
%   assume particles are on sphere with center = cent. Euler angles psi and
%   theta are calculated and stored in motl.
%
% PARAMETERS
%  INPUT
%   motl        : motivelist
%   cent        : center (3D vector)
%
%  OUTPUT
%   MOTL          motivelist with new angles
%
% last change 03/31/05 FF - docu updated
disp(['x= ' num2str(cent(1)) 'y= ' num2str(cent(2)) 'z= ' num2str(cent(3))])
npart = size(motl,2);
tmp(1:npart) = cent(1);
x = motl(8,:)-tmp;
tmp(1:npart) = cent(2);
y = motl(9,:)-tmp;
tmp(1:npart) = cent(3);
z = motl(10,:)-tmp;
tmp(:) = 90;
theta= 180/pi*atan2(sqrt(x.^2+y.^2),z); %theta
psi  = 90+180/pi*atan2(y,x); %psi
motl(18,:)=psi;
motl(19,:)=theta;