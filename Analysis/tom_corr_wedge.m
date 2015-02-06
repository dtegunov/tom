function ccf=tom_corr_wedge(a,b,wedge_a,wedge_b,flag)
%TOM_CORR_WEDGE creates ...
%
%ccf=tom_corr_wedge(a,b,wedge_a,wedge_b,flag)
%
%PARAMETERS
%
%  INPUT
%   a                   ...
%   b                   ...
%   wedge_a             ...
%   wedge_b             ...
%   flag                ...
%  
%  OUTPUT
%   ccf         		...
%
%EXAMPLE
%   ... = tom_corr_wedge(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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

    error(nargchk(0, 5, nargin, 'struct'))

    mask=tom_spheremask(ones(size(wedge_b)),30);

    wedge_all=wedge_a.*wedge_b.*mask;

    n_all = size(find(wedge_all),1); % number of voxel


% index_a=find(a.*wedge_a);
% mn = sum(sum(sum(a(index_a))))/size(index_a,1); % mean of a
% stdv = sqrt(sum(sum(sum(a(index_a).*a(index_a))))/size(index_a,1) - mn^2); % standard deviation of a
% if (stdv~=0)
%     a = (a -mn)/stdv;
%     a = a .*wedge_a;
% end;
% 
% index_b=find(b.*wedge_b);
% mn = sum(sum(sum(b(index_b))))/size(index_b,1); % mean of b
% stdv = sqrt(sum(sum(sum(b(index_b).*b(index_b))))/size(index_b,1) - mn^2); % standard deviation of a
% if (stdv~=0)
%     b = (b -mn)/stdv;
%     b= b.*wedge_b;
% end;



    n = (size(a,1)*size(a,2)*size(a,3)); % number of voxel
    mn = sum(sum(sum(a)))/n; % mean of a 
    stdv = sqrt(sum(sum(sum(a.*a)))/n - mn^2); % standard deviation of a
    if (stdv~=0)
        a = (a -mn)/stdv;
    end;

    mn = sum(sum(sum(b)))/n;
    stdv = sqrt(sum(sum(sum(b.*b)))/n - mn^2);
    if  (stdv~=0)
        b = (b - mn)/stdv;
    end;


ccf = real(ifftshift(tom_ifourier(tom_fourier(b).*conj(tom_fourier(a)))))/n_all;%perform correlation