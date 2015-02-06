function [out] = tom_isdm3file(em_name)
% TOM_ISDM3FILE tests for a file in Digital Micrograph 3-Format
%
%[out] = tom_isdm3file(em_name)
%
%
%SEE ALSO
%     TOM_ISEMFILE
%
%   created by SN 03/24/10
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

try
    [a b c]=tom_dm3read(em_name);
    out=1;
catch
    out=0;
end;