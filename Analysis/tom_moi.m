function [toi,eigvec,eigval] = tom_moi(image,flag)
% TOM_MOI calculates the moments of inertia for 3D volumes
%
%   [toi,eigvec,eigval] = tom_moi(image,flag) 
%
%   Calculates the the moments of inertia  for the input volume. The
%   function can be used as crude alignment criterion for particle or to
%   align particle along symmetry axis. 
%
%PARAMETERS
%
%  INPUT
%   image      3D array
%   flag[char] can be chosen as 'cm' or 'abs' - optional, default 'cm'
%               cm: moments of inertia relative to center of mass (default)
%               abs: absolute
%
%  OUTPUT
%   toi        tensor of inertia 3x3 matrix sum_{x,y,z}(y^2+z^2)*I  sum_{x,y,z}x*y*I   sum_{x,y,z}x*z*I  
%                                           sum_{x,y,z}x*y*I  sum_{x,y,z}(x^2+z^2)*I   sum_{x,y,z}y*z*I 
%                                           sum_{x,y,z}x*z*I  sum_{x,y,z}y*z*I   sum_{x,y,z}(x^2+y^2)*I 
%   eigvec     eigenvectors of toi - 1st eigenvector eig
%   eigval     eigenvalues of toi
%
%EXAMPLE
%   vol = tom_emread('20S_core_1.6nm.em');
%   [toi,eigvec,eigval] = tom_moi(tom_limit(-vol.Value,0,-min(min(min(vol.Value)))));
%   % moments of inertia are calculated for proteasome; density has to be
%   % inverted (dense = positive) and density origin is set. 
%
%REFERENCES
%
%SEE ALSO
%   TOM_CM
%
%   created by FF 01/19/03
%   updated by FF 05/04/05
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


error(nargchk(1,2,nargin));
if (nargin < 2)
    flag = 'cm';
end;
[s1 s2 s3] = size(image);
if strmatch(flag,'cm')
    cm = tom_cm(image);
elseif strmatch(flag,'abs')
    cm = [0, 0, 0];
else
    cm = [0, 0, 0];
end;
[x y z]=ndgrid(1-cm(1):1:s1-cm(1),1-cm(2):1:s2-cm(2),1-cm(3):1:s3-cm(3));
toi(1,1) = sum(sum(sum(image.*y.*y+image.*z.*z)));% fixed FF 04/05/05
toi(1,2) = sum(sum(sum(image.*x.*y)));
toi(1,3) = sum(sum(sum(image.*x.*z)));
toi(2,1) = toi(1,2);
toi(2,2) = sum(sum(sum(image.*x.*x+image.*z.*z)));% fixed FF 04/05/05
toi(2,3) = sum(sum(sum(image.*y.*z)));
toi(1,3) = sum(sum(sum(image.*y.*z)));
toi(3,1) = toi(1,3);
toi(3,2) = toi(2,3);
toi(3,3) = sum(sum(sum(image.*y.*y+image.*x.*x)));% fixed FF 04/05/05
[eigvec,eigval]=eig(toi);
%normalize eigvec
%for ind=1:3
%    eigvec(ind,:)= eigvec(ind,:)/sqrt(sum(eigvec(ind,:).*eigvec(ind,:)))
%end;
