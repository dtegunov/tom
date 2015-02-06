function motl_out=tom_make_clean_motl(motl_in)
%TOM_MAKE_CLEAN_MOTL throws out particles
%
%   motl_out=tom_make_clean_motl(motl_in)
%
%PARAMETERS
%
%  INPUT
%   motl_in             ...
%  
%  OUTPUT
%   motl_out            ...
%
%EXAMPLE
%   ... = tom_make_clean_motl(...);
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

motl_out_i=find(motl_in(20,:)==1);

for i=1:size(motl_out_i,2)
    motl_out(:,i)=motl_in(:,motl_out_i(i));
    motl_out(4,i)=i;
end;

