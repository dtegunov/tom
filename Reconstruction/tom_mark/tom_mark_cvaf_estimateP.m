function P = tom_mark_cvaf_estimateP(X, x)


m = size(x, 2);
n = size(X, 2);

% Drop non existing 3D points.
finite_m = all(isfinite(X(1:3,:)), 1);
X = X(1:3, finite_m);
x = x(1:2, :, finite_m);

m = size(x, 2);
n = size(X, 2);

P = zeros(3, 4, m);
P(3,4,:) = 1;

A = zeros(2*n, 8);
for (i=1:n)
    A(i*2+(-1:0), 1:8) = [[X(:,i)', [0 0 0 1 0]]; ...
                          [[0 0 0], X(:,i)', [0 1]]];
end;


D = [];

for (i=1:m)
    finite_m = squeeze(all(isfinite(x(1:2,i,:)), 1))';
    nfinite_m = sum(finite_m);
    
    if (nfinite_m < 4)
        P(:,:,i) = nan;
    else
        if (nfinite_m == n)
            if (isempty(D))
                [U, D, V] = svd(A);
                D = diag(D);
            end;
            TT_ = U' * reshape(x(1:2, i, :), [2*n, 1]);
            yi = TT_(1:8) ./ D;
            yi(abs(D) < 1e2*eps) = 0;
            PP = V*yi;
        else
            finite_m2 = reshape([finite_m; finite_m], [2*n,1]);
            [U__, D__, V__] = svd(A(finite_m2,:));
            D__ = diag(D__);
            TT_ = U__' * reshape(x(1:2, i, finite_m), [2*nfinite_m, 1]);
            yi = TT_(1:8) ./ D__;
            yi(abs(D__) < 1e2*eps) = 0;
            PP = V__*yi;
        end;        
        
        P(1:3, 1:4, i) = [[PP(1:3)', PP(7)]; ...
             [PP(4:6)', PP(8)]; ...
             [0, 0, 0,    1,]];
    end;    
end;

 

    