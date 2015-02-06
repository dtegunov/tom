function vol = tom_rotate_tform(vol, angs, interpolation, vol_center)
%TOM_ROTATE_TFORM rotates a 3d volume
%
%   vol = tom_rotate_tform(vol, angs, interpolation, vol_center)
%
%PARAMETERS
%
%  INPUT
%   vol                 3d input volume to rotate
%   angs                [phi psi theta] rotation angles according to zxz convention
%   interpolation       (optional) 'nearest','linear','cubic'
%   vol_center          (optional) [x y z], default is (size(vol)./2) + 1
%  
%  OUTPUT
%   vol                 rotated volume
%
% For more information about 3d Euler rotation visit:
% wiki -> Definitions -> Overview 3d Euler rotation -> 
% Rotation matrices, zxz and zyz conventions
% 
% For more information about sequential translation and Euler rotation visit:
% wiki -> Definitions -> Overview 3d Euler rotation -> 
% Sequential translation and rotation
%
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 08/18/06
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


if nargin < 3 || isempty(interpolation)
    interpolation = 'cubic';
end
if nargin < 4 || isempty(vol_center)
    vol_center = (size(vol) ./ 2) + 1;
end

%convert angles from degree to radians
phi=angs(1).*pi./180;
psi=angs(2).*pi./180;
theta=angs(3).*pi./180;

% use zxz matrix
T2 = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1]*...
     [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)]*...
     [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
T2(:,4) = 0;
T2(4,:) = 0;
T2(4,4) = 1;

 
%shift volume center
T1 = [1 0 0 0; 0 1 0 0; 0 0 1 0; -vol_center 1];
T3 = [1 0 0 0; 0 1 0 0; 0 0 1 0; vol_center 1];

T = T1 * T2 * T3;

tform = maketform('affine', T);

if ~(tformfwd(vol_center, tform) == vol_center); error('Sanity check failed.'); end;

R = makeresampler(interpolation, 'fill');

TDIMS_A = [1 2 3];
TDIMS_B = [1 2 3];
TSIZE_B = size(vol);
TMAP_B = [];
F = 0;
vol = tformarray(vol, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);