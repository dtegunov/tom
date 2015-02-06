function [az,elev,rho] = av3_cart2sph(r)
%CART2SPH Transform Cartesian to spherical coordinates.
%   [TH,PHI,RHO] = CART2SPH(R) transforms corresponding elements of
%   data stored in Cartesian coordinates R=[X,Y,Z] to spherical
%   coordinates (azimuth TH, elevation PHI, and radius R).  The arrays
%   X,Y, and Z must be the same size (or any of them can be scalar).
%   TH and PHI are returned in radians.
%
%   TH is the counterclockwise angle in the xy plane measured from the
%   positive x axis.  PHI is the elevation angle from the xy plane.
%
%   EXTENDED version: can handle arrays of coordinates:
%   INPUT
%   array of dim 3 x number of coordinates
%
%   See also CART2POL, SPH2CART, POL2CART.

%   L. Shure, 4-20-92.
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.12 $  $Date: 2002/04/09 00:29:46 $

rho = sqrt(sum(r.*r,1));
elev = atan2(r(3,:),sqrt(r(1,:).^2+r(2,:).^2));
az = atan2(r(2),r(1));
