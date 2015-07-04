function [ m ] = tom_mat_euler( angles )

m = tom_mat_rotz(angles(1)) * tom_mat_roty(angles(2)) * tom_mat_rotz(angles(3));

% alpha = angles(1);
% beta  = angles(2);
% gamma = angles(3);
% 
% ca = cos(alpha);
% cb = cos(beta);
% cg = cos(gamma);
% sa = sin(alpha);
% sb = sin(beta);
% sg = sin(gamma);
% cc = cb * ca;
% cs = cb * sa;
% sc = sb * ca;
% ss = sb * sa;
% 
% A = eye(3);
% A(1, 1) =  cg * cc - sg * sa;
% A(1, 2) =  cg * cs + sg * ca;
% A(1, 3) = -cg * sb;
% A(2, 1) = -sg * cc - cg * sa;
% A(2, 2) = -sg * cs + cg * ca;
% A(2, 3) = sg * sb;
% A(3, 1) =  sc;
% A(3, 2) =  ss;
% A(3, 3) = cb;

end