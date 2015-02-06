function [out] = tom_ishpcllogfile(em_name)
% TOM_ISHPCLLOGFILE tests for a log-file of the HPCL Cluster system
%
%[out] = tom_ishpcllogfile(em_name)
%
%
%SEE ALSO
%     TOM_ISEMFILE
%
%   created by SN 09/15/10
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

[n m e]=fileparts(em_name);
if isequal([m e],'logfile.txt')
    out=1;
else
    out=0;
end;
