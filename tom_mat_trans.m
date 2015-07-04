function [ m ] = tom_mat_trans( vec )

m = eye(4);
m(1, 4) = vec(1);
m(2, 4) = vec(2);
m(3, 4) = vec(3);

end

