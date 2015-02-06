function F = tom_mark_cv_get_Faffine_x(x)
% Returns the fundamental matrix F for an affine camera (parallel projection)
% minimizing the reprojection error. 
% At least n>=4 points (in general position) are needed.
% x must be a 2x2xn matrix in case of inhomogenous coordinates
% and 3x2xn in case of homogenous.
% Only finite homogen points are allowed.
% Fa satisfies x2' * Fa * x1 * 0 
% Input arguments are not checked for Inf, Nan, size, type, etc.
% From Hartley/Zissermann; ISBN 0-521-54051-8; Algorithm 14.1, pg 351
% Returns an empty matrix in case of degeneracy
msize = size(x);
if (msize(1) == 3)
    x = x(1:2, :, :) ./ repmat(x(3, :, :), [2, 1, 1]);
end;


% The 4point algorithm tom_mark_cv_get_Faffine_x4 gives the same result, as
% they both minimize the geometric distance.
% But it is slower in the current implementation :)
%
%if (msize(3) == 4)
%    F = tom_mark_cv_get_Faffine_x4(x);
%else
%    F = tom_mark_cv_get_Faffine_xn(x, msize);
%end;
%return;


Xi = squeeze(cat(1, x(:, 2, :), x(:,1,:)));
Xc = sum(Xi, 2) / msize(3);
[U, D, V] = svd((Xi - repmat(Xc, [1, msize(3)]))', 0);

N = V(:,end);

F = [[   0,    0,    N(1)];
     [   0,    0,    N(2)];
     [N(3), N(4), - N'*Xc];];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = tom_mark_cv_get_Faffine_x4(x)


A = [ [x(1,1,1), x(2,1,1), 1,        0,        0, 0, -x(1,2,1)]; ...
      [       0,        0, 0, x(1,1,1), x(2,1,1), 1, -x(2,2,1)]; ...
      [x(1,1,2), x(2,1,2), 1,        0,        0, 0, -x(1,2,2)]; ...
      [       0,        0, 0, x(1,1,2), x(2,1,2), 1, -x(2,2,2)]; ...
      [x(1,1,3), x(2,1,3), 1,        0,        0, 0, -x(1,2,3)]; ...
      [       0,        0, 0, x(1,1,3), x(2,1,3), 1, -x(2,2,3)]; ...
    ];
Ha = null(A);
Ha = [reshape(Ha(1:6), 3, 2)'; [0, 0, Ha(7)]];

l = cross(Ha * [x(:, 1, 4); 1], [x(:, 2, 4); 1]);
F = crossm([-l(2), l(1), 0]) * Ha;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = tom_mark_cv_get_Faffine_xn(x, msize)


Xi = squeeze(cat(1, x(:, 2, :), x(:,1,:)));
Xc = sum(Xi, 2) / msize(3);
[U, D, V] = svd((Xi - repmat(Xc, [1, msize(3)]))', 0);
N = V(:,end);

F = [[   0,    0,    N(1)];
     [   0,    0,    N(2)];
     [N(3), N(4), - N'*Xc];];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix of cross product cross(a,b) == crossm(a)*b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function am = crossm(a)

am = [[    0, -a(3),  a(2)]; ...
      [ a(3),     0, -a(1)]; ...
      [-a(2),  a(1),     0]];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test routine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doTest;
n = 4;

% Generate random affine cameras.
P1 = rand(3,4); P1(3, 1:3) = 0;
P2 = rand(3,4); P2(3, 1:3) = 0;


C1 = null(P1);
C2 = null(P2);

% generate random points
X = (rand(4,n));

x1 = P1 * X; x1 = x1./repmat(x1(3,:), 3, 1);
x2 = P2 * X; x2 = x2./repmat(x2(3,:), 3, 1);
x = permute(cat(3, x1(1:2, : ,:), x2(1:2, : ,:)), [1 3 2]);
msize = size(x);



Fo = crossm(P2*C1) * P2 * pinv(P1); Fo = Fo/norm(Fo);
F4 = tom_mark_cv_get_Faffine_x4(x(:,:,1:4)); F4 = F4/norm(F4);
Fn = tom_mark_cv_get_Faffine_xn(x, msize); Fn = Fn/norm(Fn);

diag(x2' * Fo * x1);
diag(x2' * F4 * x1);
diag(x2' * Fn * x1);
