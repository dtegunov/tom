function [ m ] = tom_mat_rotz( angle )

c = cos(angle);
s = sin(angle);

m = [c, s, 0;  -s, c, 0;  0, 0, 1]';

end

