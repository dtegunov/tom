function ccc = tom_ccc(a,b,flag)
% TOM_CCC calculates normalized 3d cross correlation coefficient
%
%   tom_ccc(a,b,flag)
%
% PARAMETERS
%   INPUT
%   A       array - 1d, 2d or 3d
%   B       array - dimsions as A
%   FLAG    canbe set to 'norm' for normalized CCC
%
%   OUTPUT
%   CCC     cross correlation COEFFICIENT
%
%   The cross correlation coefficient is calculated in real space, thus the
%   translation is NOT determined. If FLAG is 'norm', then the
%   normalized CCC is computed.
%
% EXAMPLE
%   im = tom_emread('proteasome.em');
%   ccc = tom_ccc(im.Value,im.Value,'norm');
%
%REFERENCES
%
% SEE ALSO
%   TOM_CORR, TOM_ORCD
%
%   last change
%   03/30/03 FF
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


if (nargin > 2)
    %    if isequal(flag,'norm')
    a= a-mean(mean(mean(a)));
    b= b-mean(mean(mean(b)));
    if (sqrt(sum(sum(sum(a.*a)))*sum(sum(sum(b.*b)))) ==0)
        ccc = sum(sum(sum(a.*b)));
    else
        ccc = sum(sum(sum(a.*b)))/sqrt(sum(sum(sum(a.*a)))*sum(sum(sum(b.*b))));
    end;
else
    ccc=(sum(sum(sum(a.*b))))./(size(a,1).*size(a,2).*size(a,3)-1) ;
end;
