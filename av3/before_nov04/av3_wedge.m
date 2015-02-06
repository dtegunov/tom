function wedge=av3_wedge(image,minangle, maxangle)
%   wedge=av3_wedge(image,minangle, maxangle)
%
%   TOM_WEDGE produces a wedge shaped array. 
%   This array can be used as a window filter in Fourier space...
%
%   image   input array - 3d 
%   angle   semi angle of missing wedge in deg
%   wedge   output - filter
%
%   FF 07/20/03
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

