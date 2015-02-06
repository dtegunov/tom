function [x y] = av3_locateonproj(coord, dims, tiltangle)
%
%
%
%   COORD       coordinates in volume
%   DIMS        dimensions of tomgram
%   TILTANGLE   tiltangle (in deg)
%
%   X           X in micrograph - respective to center
%   Y           Y in micrograph - respective to center
diff_coord = coord - floor(dims/2)+1 ;%calculate respective to center
tiltangle = tiltangle/180*pi;
x = cos(tiltangle)*coord(1) -sin(tiltangle)*coord(3);
y = coord(2);