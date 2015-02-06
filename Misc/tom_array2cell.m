function [cell]=tom_arry2cell(array)
%TOM_ARRAY2CELL creates ...
%
%   [cell]=tom_arry2cell(array)
%
%PARAMETERS
%
%  INPUT
%   array               ...
%  
%  OUTPUT
%   cell            	...
%
%EXAMPLE
%   .. = tom_arry2cell(...);
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


for i=1:size(array,2)
    cell{i}=array(i);
end;