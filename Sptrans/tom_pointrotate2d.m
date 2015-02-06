function trans_p=tom_pointrotate2d(point,angle,center)
%TOM_POINTROTATE rotates point
%
%   r = tom_pointrotate(point,angle,center)
%
%   A vector in 3D is rotated around the origin = [0 0 0]. The puropose is
%   for example to predict the location of a point in a volume after
%   rotating it with tom_rotate3d. Take care that the coordinates are with
%   respect to the origin!
%
%PARAMETERS
%
%  INPUT
%   point                   point
%   angle                 angle in deg 
%   center                 (optional) center default (0 0)
%
%  
%  OUTPUT
%   trans_p                 rotated coordinates
%
%EXAMPLE
%   p = [12 13]
%   r = tom_pointrotate2d(p,30,[11 15])
%
%REFERENCES
%
%SEE ALSO
%   tom_rotate
%
%   created by fb(eckster)
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



if nargin < 3
    center=[0 0];
end;

m=[cosd(angle) -sind(angle); sind(angle) cosd(angle)];

%shift da coordinat system!
trans_p=point-center;

%rotate point
try
    trans_p=m*trans_p;
catch
    trans_p=m*trans_p';
    trans_p=trans_p';
end;


%shift coordinate system back!
trans_p=trans_p+center;