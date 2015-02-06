function v_out=tom_cell2num(v_in)
%TOM_CELL2NUM creates ...
%
%   v_out=tom_cell2num(v_in)
%
%PARAMETERS
%
%  INPUT
%   v_in                ...
%  
%  OUTPUT
%   v_out               ...
%
%EXAMPLE
%   .. = tom_cell2num(...);
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


for i=1:size(v_in,1)
    v_out=str2num(v_in);
end;