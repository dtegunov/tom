function [inliers] = tom_mark_cvaf_inliersRANSAC(x, threshold, maxTrials, nocheckx)
%TOM_MARK_CVAF_INLIERSRANSAC returns the index of the inliers using ransac.
%
%     [inliers] = tom_mark_cvaf_inliersRANSAC(x, threshold, maxTrials, nocheckx)
%
% Computes the set of inliers, of point correspondences in M views using
% RANSAC
% 
%PARAMETERS
%  INPUT
%    x: a 2xMxN matrix of N inhomogen (finite) image points or a 3xMxN matrix
%       with homogen image points. if "nocheck" is true, all points must be
%       finite and (in case of homogenous points) the last row must
%       be normalized to 1.
%       Those points which contain at least one, non finite point are
%       ignored.
%    threshold: Value to decide whether point is inlier or not.
%    maxTrails: maximum number of tries in RANSAC.
%    nocheckx: Optional parameter for checking the parameter x. if false or
%       omitted, the pointset x is checkt, so that
%       point, which are not finite in all views are rejected and in case of homogenous 
%       coordinates, the last row can differ from 1. Useful to save
%       unnecessary work.
%      
%  OUTPUT
%    inliers: indices of inliers (1...N). If not enough points are given, 
%       [] is returned.
%
% Uses RANSAC algorithm: Calles TOM_MARK_CVAF_CAM_FROM_X with a random sample of
% 4 points to obtain a camera configuration. The obtained cameras are used
% to triangulate all points and reproject them. A point is an inlier, if 
% the maximum!!! from the geometric distance between its positions in the M
% views and the reprojected points is smaller then the input argument
% threshold.
%
%REFERENCES
%   Hartley/Zissermann, Multiple View Geometry in Computer Vision. ISBN 0-521-54051-8
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


msize = size(x);
if (length(msize) == 2)
    msize(3) = double(~isempty(x));
end;

if (~exist('nocheckx', 'var') || ~nocheckx)
    
    idxdefined = squeeze(all(all(isfinite(x), 1), 2));
    if (msize(1) == 3)
        idxdefined = idxdefined & squeeze(all(abs(x(3,:,:)) > eps*1e3, 2));
        x = x ./ repmat(x(3,:,:), 3, 1);
    end;
    
    x = x(:,:,idxdefined);
    msize(3) = size(x, 3);
end;

if (msize(3) < 4)
    inliers = [];    
    return;
end;

% x1 and x2 are 'stacked' to create a 6xN array for ransac
[P, inliers] = tom_mark_kovesi_ransac(reshape(x(1:2, :, :), [2*msize(2), msize(3)]), @fittingfn, @distfn, @degenfn, 4, threshold, 100, maxTrials);

if (exist('idxdefined', 'var'))
    i = find(idxdefined);
    inliers = i(inliers);
end;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fittingfn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = fittingfn(x)

x = reshape(x, [2, size(x,1)/2, size(x,2)]);
P = tom_mark_cvaf_cam_from_x(x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inliers, P] = distfn(P, x, t);

x = reshape(x, [2, size(x,1)/2, size(x,2)]);

dist = getDist(P, x);

inliers = find(max(dist, [], 2) < t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = getDist(P, x)

X = tom_mark_cvaf_triX(P, x);

m = size(P, 3);
n = size(X, 2);
xr = nan(3, m, n);
for (i=1:m)
    xr(1:3, i, :) = P(:,:,i) * X;
end;

dist = sqrt(sum((xr(1:2, :, :) - x) .^ 2, 1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (Degenerate!) function to determine if a set of matched points will result
% in a degeneracy in the calculation of a fundamental matrix as needed by
% RANSAC.  This function assumes this cannot happen...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = degenfn(x)
r = false;
    