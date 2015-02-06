function [out]=tom_path2cell(path,basename,ext,interval)
%TOM_PATH2CELL creates ...
%
%   [out]=tom_path2cell(path,basename,ext,interval)
%
%PARAMETERS
%
%  INPUT
%   path                ...
%   basename            ...
%   ext                 ...
%   colintervalor       ...
%  
%  OUTPUT
%   out                 ...
%
%EXAMPLE
%   ... = tom_path2cell(...);
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

for i=interval(1):interval(2)
    out{i}=[path '/' basename num2str(i) ext];
end;