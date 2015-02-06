function [H] = tom_mark_cvaf_metricX(X, P, idxpro);

Xo = X;
X = X(:, idxpro);

nm = size(X, 2);
x = P(1:3,1:4) * X(1:4, :);


H = eye(4);


[U, D, V] = svd(X', 0);
dD = diag(D);

y = ((x(1,:) * U)') ./ dD;
y(abs(dD) < 1e3*eps) = 0;
H(1, 1:4) = V * y;

y = ((x(2,:) * U)') ./ dD;
y(abs(dD) < 1e3*eps) = 0;
H(2, 1:4) = V * y;

a = H(1,1);
b = H(1,2);
c = H(1,3);
d = H(1,4);
e = H(2,1);
f = H(2,2);
g = H(2,3);
h = H(2,1);
a_ = P(1,1);
b_ = P(1,2);
c_ = P(1,3);
d_ = P(1,4);
e_ = P(2,1);
f_ = P(2,2);
g_ = P(2,3);
h_ = P(2,4);

A = ...
[[ - c*f + b*g - b_*g + c_*f, - a*g + e*c - c_*e + a_*g, -a_*f + b_*e + a*f - e*b, 0 ]; ...
 [ - c_*b + b_*c, + c_*a - a_*c, + a_*b - b_*a, 0 ]; ...
 [ - d_*b*g - c_*d*f - b_*c*h + b_*d*g + c_*b*h + d_*c*f, + d_*a*g - d_*e*c + a_*c*h - c_*a*h - a_*d*g + c_*e*d, - d_*a*f + b_*h*a - a_*h*b + a_*d*f + d_*e*b - b_*e*d, + b_*e*c - b_*g*a + a_*g*b - a_*c*f + c_*a*f - c_*e*b ]; ...
 [ - f_*g + g_*f, + e_*g - g_*e, - e_*f + f_*e, 0 ]; ...
 [ + f_*c - g_*b - c*f + b*g, + g_*a - a*g - e_*c + e*c, + e_*b - f_*a + a*f - e*b, 0 ]; ...
 [ - f_*c*h + f_*d*g - g_*d*f + g_*b*h + h_*c*f - h_*b*g, + h_*a*g - g_*a*h + g_*e*d - h_*e*c + e_*c*h - e_*d*g, - h_*a*f - e_*h*b + e_*d*f + f_*h*a - f_*e*d + h_*e*b, + f_*e*c - f_*g*a + g_*a*f - g_*e*b + e_*g*b - e_*c*f ]];

A = (A(:,1:3));
[U, D, V] = svd(A, 0);
dD = diag(D);

y = ([0 0 1 0 0 1] * U)' ./ dD;
y(abs(dD) < 1e3*eps) = 0;
H(3,1:4) = [(V * y)', 1];


return;

X = Xo;
l = max(X(3,:)) - min(X(3,:));
if (l > 0)
    H(3, 3) = 0.3*mean(max(X(1:2,:), [], 2) - min(X(1:2,:), [], 2)) / l;
    H(3, 4) = - min(X(3,:)) * H(3, 3);
end;

