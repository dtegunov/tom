function [ m4 ] = tom_mat_euler4( angles )

m = tom_mat_rotz(angles(1)) * tom_mat_roty(angles(2)) * tom_mat_rotz(angles(3));
m4 = eye(4);
m4(1:3, 1:3) = m;

end