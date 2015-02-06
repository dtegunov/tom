function [P, x_proj, X_trian, x_reproj] = tom_mark_cvaf_alignment3d_reproj(Tiltangles, psi, tx, ty, Imdim, X, x, ret_inhomogen)
%TOM_MARK_CVAF_ALIGNMENT3D_REPROJ
% Converts the output from tom_alignment3d (tom_mark_cvaf_alignment3d) into projection
% matrices P. P has dimension 3x4xN where N is the number of
% images. P is a set of homogenous affine projection matrices.
% Thus the last row of every projection matrix is P(3,1:4,i) = [0 0 0 1]
%
% Notice, that if some parameters (e.g. tx or ty) are nan, the returned
% projection matrix P contains some unknown entries nan. Set them as you
% wish:
%     P(1:3,1:4,~all(all(isfinite(P), 1), 2)) = nan;
%     P(~isfinite(P)) = 0;
%
%INPUTPARAMETERS:
% Tiltangles is the (known) angle for every image in DEG!!!
%   It is exactly what was handled to tom_alignment3d in Matrixmark(1,1:N,1)
% psi is the 2nd output parameter from tom_alignment 3D. It is in RAD!!!!
% tx and ty are vectors of length N with the shifts for every image. They
%   are returned from tom_alignment3d as Matrixmark(7,1:N,1) and
%   Matrixmark(8,1:N,1) respectively.
% Imdim is the last parameter handled to tom_alignment3d.
%
%EXAMPLE: reprojecting the computed 3d coordinate from tom_alignment3d
%   irefmark = 1;        % the reference marker.
%   imintilt = 10;       % Reference projection
%   irefX = rand(3,1) * imsize; % 3D coordinate of reference marker. 
%   [Matrixmark, psi, rms, x, y, z]  = ...
%       tom_alignment3d(Matrixmark, irefmark, imintilt, irefX, imsize);
%   X = cat(1,x,y,z) + (imsize/2 + 1); % The 3D coordinates. Notice:
%                                      % X(1:3,irefmark) == irefX
%   P = tom_mark_cvaf_alignment3d_reproj(Matrixmark(1,:,1), psi, 
%                          Matrixmark(7,:,1), Matrixmark(8,:,1), imsize);
%   X = cat(X, ones(1,size(X,2));
%   for (i=1:size(Matrixmark,2))
%       x(1:3,i,1:length(x)) = P(1:3,1:4,i) * X;
%   end;
%   % x are the reprojected points. In case of exact marker points
%   % it holds x(1:2,:,:) == Matrixmark(2:3,:,:)
%
%SEE ALSO
%   TOM_ALIGNMENT3D
%   
%   created by Thomas Haller, 27. july 2007
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

ni = length(Tiltangles);
tx = reshape(tx, [1, ni]);
ty = reshape(ty, [1, ni]);

Tiltangles = Tiltangles * pi() / 180;

P = zeros(3,4,ni);

spsi = sin(psi);
cpsi = cos(psi);

cTiltangles = cos(Tiltangles);
sTiltangles = sin(Tiltangles);

for (i=1:ni)
    P(1:2,1:3,i) = ...
         [ spsi^2 * cTiltangles(i)                                + cpsi^2, ...
           spsi   * cpsi               * (1-cTiltangles(i)), ...
           spsi   * sTiltangles(i); ...
           spsi   * cpsi               * (1-cTiltangles(i)), ...
           cpsi^2 * cTiltangles(i)                                + spsi^2, ...
          -cpsi   * sTiltangles(i)];
end;


P(1:2, 4, 1:ni) = [tx; ty] - squeeze(sum(P(1:2,1:3,:), 2) * (Imdim/2+1)) + (Imdim/2+1);
P(3,4,:) = 1;
 

if (~exist('X', 'var'))
    return;
end;
if (~exist('ret_inhomogen','var'))
    ret_inhomogen = true;
end;

if (isempty(X))
    x_proj = [];
else
    if (size(X,1) == 3)
        X(4,:) = 1;
    elseif (size(X,1)~=4)
        error('X must be a vector of 3D points with size 3xN (or 4xN, last row equal 1).');
    end;
    x_proj = nan(3,ni,size(X,2));
    for (i=1:ni)
        x_proj(1:3, i, :) = P(1:3,1:4,i) * X;
    end;
    if (ret_inhomogen)
        x_proj = x_proj(1:2,:,:);    
    end;
end;

if (~exist('x', 'var'))
    return;
end;
if (isempty(x))
    X_trian = [];
    x_reproj = [];
else
    if (size(x,1)~=2 || size(x,2)~=ni)
        error(['x must be a inhomogenous vector of 2D points with size 2x' num2str(ni) 'xN.']);
    end
    X_trian = tom_mark_cvaf_triX(P, x);
    if (nargout > 3)
        x_reproj = nan(3,ni,size(X_trian,2));
        for (i=1:ni)
            x_reproj(1:3, i, :) = P(1:3,1:4,i) * X_trian;
        end;
        if (ret_inhomogen)
            x_reproj = x_reproj(1:2,:,:);    
        end;
    end;
    if (ret_inhomogen)
        X_trian = X_trian(1:3,:,:);
    end;
end;





