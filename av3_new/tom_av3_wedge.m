function wedge=av3_wedge(image,minangle, maxangle)
%TOM_AV3_WEDGE produces a wedge shaped array.
%
%   wedge=av3_wedge(image,minangle, maxangle)
%
%PARAMETERS
%
%  INPUT
%   image               input array - 3d
%   minangle            semi angle of missing wedge in deg
%   maxangle            semi angle of missing wedge in deg
%  
%  OUTPUT
%   wedge               output - filter
%
%EXAMPLE
%   ... = av3_wedge(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FF 07/20/03
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


warning off MATLAB:divideByZero;
maxangle = maxangle*pi/180;
minangle = minangle*pi/180;
[dimx, dimy, dimz] = size(image);
[x,y,z] = ndgrid(-floor(dimx/2):-floor(dimx/2)+dimx-1,-floor(dimy/2):-floor(dimy/2)+dimy-1,-floor(dimz/2):-floor(dimz/2)+dimz-1);
wedge = ones(dimx, dimy, dimz);
% 1st angles
ind = find(tan(pi/2-maxangle) > (x)./(z));
% merge with negative
ind2=find(tan(-pi/2-minangle) < (x)./(z));
ind = intersect(ind, ind2);
wedge(ind)=0;

