function [P, X] = tom_mark_cvaf_cam_from_x(x0)
%TOM_MARK_CVAF_CAM_FROM_X computes the affine cameras from a set of
%projections x
%
%   [P, X] = tom_mark_cvaf_cam_from_x(x)
%
%PARAMETERS
%  INPUT
%   x : a 2xMxN matrix of image points. N is the number of points
%       which are projected onto the M views. All points must be
%       finite (i.e. ~nan and ~inf)
%
%  OUTPUT
%   P : a 3x4xM matrix containing the affine cameras for the M views.
%       the n-th 3D point is projected onto the image i-th by the equation
%       x(1:2, i, n) ~ P(1:3, 1:4, i) * X(1:4, n)
%       P are affine cameras, that means, the last row P(3,1:4,i) is [0 0 0 1]
%       for all views i.
%   X : a 4xN matrix with the corresponding 3D points
%
% Notice that the camera matrices P (and the 3D Points X) are up to an 
% affine ambiguity. That means, that the REAL cameras and the REAL 3D points
% are undetremined up to an (unknown) affine transformation H which  holds
% Preal = P * inv(H) and Xreal = H*X
%
% EXAMPLE
%   execute TOM_MARK_CVAF_CAM_FROM_X without parameters...
%
%REFERENCES
%   Implements algorithm 18.1, page 437 from Hartley/Zissermann, Multiple
%   View Geometry in Computer Vision. ISBN 0-521-54051-8
%
%SEE ALSO
%   -
%
%   created by Thomas Haller, 20. july 2007
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom
%%


if (nargin < 1)
    warning('Called with no parameter. DO TEST');
    doTest;
    return;
end;


% First prepare the input points x0....

if (size(x0,1) == 2)
    idxundefm = squeeze(any(any(isnan(x0) | isinf(x0), 1), 2));
    x = x0(:, :, ~idxundefm); 
elseif (size(x0,1) == 3)
    idxundefm = squeeze(any(any(isnan(x0(1:2,:,:)) | isinf(x0(1:2,:,:)), 1) | (abs(x0(3,:,:))<1e3*eps), 2));
    x = x0(:, :, ~idxundefm);
    x = x(1:2,:,:) ./ repmat(x(3,:,:), [2, 1, 1]);
else
    error('x must be a 2xIxN (inhomogenous) or a 3xIxN (homogenous) coordinates.');
end;


t = mean(x, 3);

x = x - repmat(t, [1, 1, size(x,3)]);

W = reshape(x, 2*size(x,2), size(x, 3));

if (nargout >= 2)
    [U,D,V] = svd(W);

    M = [D(1,1)*U(:,1), D(2,2)*U(:,2), D(3,3)*U(:,3)];
    X = nan(4, size(x0, 3));
    X(1:4, ~idxundefm) = [V(:,1), V(:,2), V(:,3), ones(size(x,3),1)]';
else
    [U,D] = svd(W, 0);

    M = [D(1,1)*U(:,1), D(2,2)*U(:,2), D(3,3)*U(:,3)];
end;

P = zeros(3, 4, size(x, 2));
for (i=1:size(x,2))
    P(1:2, 1:3, i) = M(2*(i-1)+(1:2), 1:3);
    P(1:2, 4, i) = t(1:2, i);
end;
P(3, 4, :) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test routine with random values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doTest



% Size of Testdata
M = 100; N = 100;
% Random test points and random matrices.
X = rand(3, N); X(4, :) = 1;
P = rand(3,4,M); P(3,1:3,:) = 0; P(3,4,:) = 1;
% Project points
x = nan(3, M, N);
for (i=1:M)
    x(1:3, i, 1:N) = P(1:3, 1:4, i) * X(1:4,1:N);
end;
% Add gausian noise...
xnoise = x(1:2,:,:) + randn(2, M, N) .* repmat(max(sum((x - repmat(mean(x, 3), [1, 1, N])) .^ 2, 1), [], 3), [2, 1, N]) / 10;
xnoise(3,:) = 1;

% Compute affine cameras.
[P2, X2] = tom_mark_cvaf_cam_from_x(xnoise(1:2,:,:));

% Reproject 
xnoisereproj = nan(3, M, N);
for (i=1:M)
    xnoisereproj(1:3, i, 1:N) = P2(1:3, 1:4, i) * X2(1:4,1:N);
end;

% Deviation of distance between noisy data and reprojected data.
tom_dev(sqrt(sum((xnoise - x) .^ 2, 1)));
tom_dev(sqrt(sum((xnoisereproj - x) .^ 2, 1)));

% Plot points for one view
i=1;pN = randperm(N); pN = pN(1:100);
figure(1); cla, hold on; axis equal;
plot(squeeze([x(1,i,pN),xnoise(1,i,pN),xnoisereproj(1,i,pN)]), squeeze([x(2,i,pN),xnoise(2,i,pN),xnoisereproj(2,i,pN)]), 'm-');
plot(squeeze(x(1,i,pN)), squeeze(x(2,i,pN)),'.b', squeeze(xnoise(1,i,pN)), squeeze(xnoise(2,i,pN)),'.g', squeeze(xnoisereproj(1,i,pN)), squeeze(xnoisereproj(2,i,pN)),'.r');






