function [a,b,c,d,e]=tom_dev(A,suppress,mask)
%TOM_DEV Calculates the mean, max, min, standard-deviation, variance
%   of an input image
%
%   [a,b,c,d,e]=tom_dev(A,suppress,mask)
%
%PARAMETERS
%  INPUT
%   A           data.
%   suppress    (opt) verbose flag (meaning suppress == 'noinfo' -> will suppress output)
%   mask        (opt) statistic is only calculated inside mask  
%
%
%
%  OUTPUT
%   a           Mean value
%   b           maximum value
%   c           minimum value
%   d           standard deviation
%   e           Variance of image
%
% EXAMPLE
%   [mean, max, min, std, variance] = tom_dev(IN,flag)
%   im = tom_emread('proteasome.em');
%   tom_dev(im.Value);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AF 03/18/02
%   updated by WDN 07/11/03
%   updated by 05/23/06
%   updated by fb(eckster) statistics inside Mask 05/23/06
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

%error(nargchk(1,2,nargin))
if nargin < 3
    mask=ones(size(A));
end;

idx=find(mask > 0);
A=double(A); % SN, different mean value for single!!! -> cast to double
%A=A.*mask;

[s1,s2,s3,s4]=size(A(idx));
a=sum(sum(sum(sum(A(idx)))))./(s1.*s2.*s3.*s4);
b=max(max(max(max(A(idx)))));
c=min(min(min(min(A(idx)))));
d=std(A(idx),1);
e=d.^2;
if nargin <2
    f=sprintf('Mean= %g,  Max= %g,  Min= %g,  Std= %g,  Variance= %g', a,b,c,d,e);disp(f);
elseif nargin < 4
    switch suppress
        case 'noinfo'
            %no action
        otherwise
            f=sprintf('Mean= %g,  Max= %g,  Min= %g,  Std= %g,  Variance= %g', a,b,c,d,e);disp(f);
    end
end
