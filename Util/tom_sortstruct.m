function structure = tom_sortstruct(structure, fieldname, order)
%TOM_SORTSTRUCT creates ...
%
%   structure = tom_sortstruct(structure, fieldname, order)
%
%PARAMETERS
%
%  INPUT
%   structure           ...
%   fieldname           ...
%   order               ...
%  
%  OUTPUT
%   structure           ...
%
%EXAMPLE
%   ... = tom_sortstruct(...);
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


if nargin < 3
    order = 'ascend';
end

sortmatrix = zeros(1,size(structure,2));

for i=1:size(structure,2)
    sortmatrix(i) = structure(1,i).(fieldname);
end

[sorted,idx] = sort(sortmatrix,2,order);

for i=1:size(structure,2)
    outstruct(1,i) = structure(1,idx(i)); 
end

structure = outstruct;
