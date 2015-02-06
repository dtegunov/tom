function align_out=tom_av2_change_align2d_basedir(align2d,new_basedir)
%TOM_AV2_CHANGE_ALIGN2D_BASEDIR creates ...
%
%   align_out=tom_av2_change_align2d_basedir(align2d,new_basedir)
%
%PARAMETERS
%
%  INPUT
%   align2d             ...
%   new_basedir         ...
%  
%  OUTPUT
%   align_out           ...
%
%EXAMPLE
%   ... = tom_av2_change_align2d_basedir(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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


align_out=align2d;
for i=1:size(align2d,2)
     [old_path filen ext]=fileparts(align2d(1,i).filename);
     
     align_out(1,i).filename=[new_basedir '/' filen ext];
     align_out(1,i).isaligned=0; 
end;