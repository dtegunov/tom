function [ m ] = tom_mat_roty( angle )

c = cos(angle);
s = sin(angle);

m = [c, 0, -s;  0, 1, 0;  s, 0, c]';

end

