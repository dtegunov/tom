function tom_make_motl_list(emname, num)
%TOM_MAKE_MOTL_LIST creates ...
%
%   tom_make_motl_list(emname, num)
%
%PARAMETERS
%
%  INPUT
%   emname              ...
%   num                 ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_make_motl_list(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN (author date)
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

l=zeros(11,num);
for i=1:num
    l(4,i)=i;
end;
tom_emwrite(emname,l);