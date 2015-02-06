function H = tom_markh_H_from_X_3d_lin(X1, X2)
% Computes a projective transformation between
% 3D points X1 to X2. X2 = H * X1.
% Currently only finite points reguarded.
% Uses the DLT algorithm, which minimizes
% algebraical error!


[r,c] = size(X1);

if ([r,c] ~= size(X2))
   error ('Input point sets are different sizes!')
elseif (c < 1)
   error('No input points given');
end

if (r == 3)
   C1 = getCondition(X1);
   C2 = getCondition(X2);

   X1 = [X1; ones(1, c)];
   X2 = [X2; ones(1, c)];
else
   X1_not_ideal = abs(X1(4,:)) > 1e4*eps;
   X2_not_ideal = abs(X2(4,:)) > 1e4*eps;
   sX1_not_ideal = sum(X1_not_ideal);
   sX2_not_ideal = sum(X2_not_ideal);
   doInverse = sX2_not_ideal < sX1_not_ideal;
   if (doInverse)
       tmp = X1; X1 = X2; X2 = tmp;
       tmp = X1_not_ideal; X1_not_ideal = X2_not_ideal; X2_not_ideal = tmp;
       sX2_not_ideal = sX1_not_ideal;
   end;
   if (sX2_not_ideal < c)
       X1 = X1(:, X2_not_ideal);
       X2 = X2(:, X2_not_ideal);
       if (isempty(X1))
           error('Only ideal points');
       end;
       warning([num2str(c-sX2_not_ideal) ' ideal points ignored.']);
       c = size(X1, 2);
   end;

   X1c = X1(1:4,abs(X1(4,:))>3e4*eps);
   X1c = X1c(1:3,:) ./ repmat(X1c(4,:),3,1);
   X2 = X2 ./ repmat(X2(4,:),4,1);


   C1 = getCondition(X1c);
   C2 = getCondition(X2(1:3,:));
end

if (c < 5)
   error('At least 5 points needed');
end;


X1 = C1 * X1;
X2 = C2 * X2;

A = zeros(3*c, 16);
OOOO  = zeros(1,4);



for (i=1:c)
   pX1 = X1(:,i)';
   pX2 = X2(:,i)';
   A(3*i+(-2:0),1:16) = [
           [pX1,  OOOO, OOOO, -pX2(1)*pX1] ;
           [OOOO,  pX1, OOOO, -pX2(2)*pX1] ;
           [OOOO, OOOO, pX1,  -pX2(3)*pX1]
       ];

end

% Extract nullspace
[u,s,v] = svd(A, 0);

s = diag(s);
nullspace_dimension = sum(s < eps * s(1) * 1e3);
if (nullspace_dimension > 1)
   warning('Nullspace is a bit roomy...');
end;

H = inv(C2) * reshape(v(:,end), 4, 4)' * C1;

if (doInverse)
   H = inv(H);
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X must be inhomogenous 3xN vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = getCondition(X, isotropic)

if (isempty(X))
   C = eye(4);
   return;
end;
mX = mean(X, 2);
sX = std(X, 0, 2);

sX = sX + (sX==0);

if (nargin==1 || ~isotropic)
   sC = diag(sqrt(2)./sX);
else
   sC = eye(3) * sqrt(2)/mean(sX);
end;

C = [[sC, -sC*mX]; [0 0 0 1]];

