function [x y] = av3_locateonproj(coord, dims, tiltangle)
% AV3_LOCATEONPROJ computes 3D coordinates to 2D coordinates in projection
%
%  [x y] = av3_locateonproj(coord, dims, tiltangle)
%
% PARAMETERS
%  INPUT
%   COORD       coordinates in volume (3dim array)
%   DIMS        dimensions of tomgram (e.g. [512 512 128])
%   TILTANGLE   tiltangle (in deg)
%
%  OUTPUT
%   X           X in micrograph - respective to center
%   Y           Y in micrograph - respective to center
%
%  last change 03/31/05 FF - docu updated
%
diff_coord = coord - floor(dims/2)+1 ;%calculate respective to center
tiltangle = tiltangle/180*pi;
x = cos(tiltangle)*coord(1) -sin(tiltangle)*coord(3);
y = coord(2);