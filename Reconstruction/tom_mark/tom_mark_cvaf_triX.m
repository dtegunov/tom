function X = tom_mark_cvaf_triX(P, x)
% Triangulates the projections x of inhomogen points 
% (2xMxN) from M affine cameras P (3x4xM)
% Returns a homogen 4xN vector with last coordinate 1.
% Minimizes reprojection error per point (linear)
%
% Projections x(:,:,i) with nan coordinates are ignored and the corresponding
% 3D point X(1:4,i) is triangulated by the remaining. At least two known 
% Positions are needed. If less then 2 positions are given, X(:,i) is set
% to nan. 

Pexist = squeeze(all(all(isfinite(P), 1), 2));
if (~all(Pexist))
    P = P(:,:,Pexist);
    x = x(:,Pexist,:);
end;

m = size(P,3);
n = size(x, 3);

if (m < 2)
    X = nan(4, n);
    return;
end;

A = nan(2*m, 3);
T = nan(2*m, 1);
for (i=1:m)
    A(2*(i-1)+(1:2), 1:3) = P(1:2,1:3,i) / P(3, 4, i);
    T(2*(i-1)+(1:2), 1) = -P(1:2, 4, i) / P(3, 4, i);
end;

X = nan(4, n);
U = [];
xmask = all(isfinite(x(1:2, :, :)), 1);
xmask2 = repmat(xmask, [2,1,1]);
nxmask = squeeze(sum(xmask, 2))';
for (i=1:n)
    nxmaski = nxmask(i);
    if (nxmaski < 2)
        %X(1:3, i) = nan;
    elseif (nxmaski == m)
        if (isempty(U))
            [U, D, V] = svd(A);
            D = diag(D);
        end;
        TT_ = U' * (T + reshape(x(1:2, :, i), [2*m, 1]));
        yi = TT_(1:3) ./ D;
        X(1:4, i) = [V*yi; 1];
    else
        xmask2i = reshape(xmask2(:,:,i), [2*m, 1]);
        [U__, D__, V__] = svd(A(xmask2i, 1:3));
        D__ = diag(D__);
        TT_ = U__' * (T(xmask2i) + reshape(x(1:2, xmask(1,:,i), i), [2*nxmaski, 1]));
        yi = TT_(1:3) ./ D__;
        yi(abs(D__) < 1e3*eps) = 0;
        X(1:4, i) = [V__*yi; 1];
    end
end;
    


