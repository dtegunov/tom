function [R C z_0] = ellipsedet(R_1,z_1,R_2,z_2,R_3,z_3)
%ELLIPSEDET           determine ellipse parameters from parameters
%
%   [R C z_0] = ellipsedet(R_1,z_1,R_2,z_2,R_3,z_3)
%
% routine determines ellipse parameters z_0, C, and R
%
%
%   ellipse is assumed to have the form:
%
%  ( (x - x_0)/R )^2 + ( (y - y_0)/R )^2 + ( (z - z_0)/C )^2 = 1
%

%   3 x-y Slices:
%   1)      R_1 + (R/C)^2 * (z_1-z_0)^2 = R^2
%   2)      R_2 + (R/C)^2 * (z_2-z_0)^2 = R^2
%   3)      R_3 + (R/C)^2 * (z_3-z_0)^2 = R^2
%
%   R_1, R_2 and R_3: squares of interactively determined radii in slices
%   take care - in input assumed as real radii, therefore squared for
%
% PARAMETERS
%  INPUT
%   R_1, R_2, R_3   radii of circles in 3 different xy slices
%   z_1, z_2, z_3   corresponding z's
%
%  OUTPUT
%   R               xy-radius of ellipse (see eq.)
%   C               z-radius of ellipse (see eq.)
%   z_0             z-height of ellipse (see eq.)
%
% SEE ALSO
%   ELLIPSEPOS2ANGLE
%
%   last change 03/31/05 FF

%   calculations
R_1 = R_1^2;
R_2 = R_2^2;
R_3 = R_3^2;
%   eps_sq = (R/C)^2
eps_sq = (R_2 - R_3)/( (z_2-z_3)*(z_1+z_2) ) + (R_2 - R_1)/( (z_1-z_2)*(z_1+z_2) );
z_0 = 1/2* ( 1/eps_sq *(R_1-R_2)/(z_1-z_2) + z_1 + z_2 );
R = sqrt(R_1 + eps_sq*(z_1-z_0)^2);
C = R/sqrt(eps_sq);
