function tom_av3_particleemlist_2_particles(particle_list,cut_size,vol,outputpath)
%TOM_AV3_PARTICLEEMLIST_2_PARTICLES creates ...
%
%   tom_av3_particleemlist_2_particles(particle_list,cut_size,vol,outputpath)
%
%PARAMETERS
%
%  INPUT
%   particle_list       ...
%   cut_size            ...
%   vol                 ...
%   outputpath          ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av3_particleemlist_2_particles(...);
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


for i=1:size(particle_list,2)
    c=particle_list(6:8,i)'-16;
    part=tom_cut_out(vol,c,cut_size,'no-fill');
    
   tom_emwrite([outputpath num2str(i) '.em'],part);
    
end;

