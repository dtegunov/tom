function cm =tom_cm(image)
%TOM_CM Calculates the Center of Mass for 1D, 2D or 3D data
%   CM=TOM_CM(A) Calculates the coordinates of the center of mass
%   for the input image. Note that the first element of the result
%   is the horizontal coordinate (or x-coordinate) of the center of mass,
%   and the second element is the vertical coordinate (or y-coordinate).
%   The result does not have to be an integer.
%
%   cm =tom_cm(image)
%
%PARAMETERS
%
%  INPUT
%   image               ...
%  
%  OUTPUT
%   cm      			...
%
%   Example
%
%           0  0  0  0  0  0
%           0  0  0  0  0  0
%           0  0  0  0  0  0
%       a = 0  0  0  1  1  1
%           0  0  0  1  1  1
%           0  0  0  1  1  1
%
%       [x y] = tom_cm(a)
%
%       x = 
%           5
%
%       y =
%           5
%
%REFERENCES
%
%   See also TOM_LIMIT, TOM_MIRROR, TOM_PEAK, TOM_SPC
%
%   08/16/02   AL
%   01/19/03   revised version FF
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


image=image+abs(min(min(min(image))));
if nargin<1
    error('Not enough Input Arguments');
elseif nargin == 1
    [s1 s2 s3] = size(image);
    norm = sum(sum(sum(image)));
    [x y z]=ndgrid(1:1:s1,1:1:s2,1:1:s3);
    cm(1) = sum(sum(sum(image.*x)))/norm;
    cm(2) = sum(sum(sum(image.*y)))/norm;
    cm(3) = sum(sum(sum(image.*z)))/norm;
else
    error('Too many Input Arguments');
end

