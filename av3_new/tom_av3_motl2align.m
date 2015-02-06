function Align=tom_av3_motl2align(emfile);
%TOM_AV3_MOTL2ALIGN creates ...
%
%   Align=tom_av3_motl2align(emfile)
%
%PARAMETERS
%
%  INPUT
%   emfile              ...
%  
%  OUTPUT
%   Align               ...
%
%EXAMPLE
%   ... = tom_av3_motl2align(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
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


motl=tom_emread(emfile);
motl=motl.Value;

partname='/fs/montreal/pub/26S/alignment_tomography/data_Oana_1/PART_';

for i=1:size(motl,2)

    Align(i).ID=[partname num2str(motl(1,i))];
    Align(i).Source='Oana_1';
    Align(i).CCF=motl(2,i);
    Align(i).Angle.Phi=motl(3,i);
    Align(i).Angle.Psi=motl(4,i);
    Align(i).Angle.Theta=motl(5,i);
    Align(i).Shift.X=motl(6,i);
    Align(i).Shift.Y=motl(7,i);
    Align(i).Shift.Z=motl(8,i);
    Align(i).Class=0;
    Align(i).ProjectionClass=0;
    Align(i).Tiltseries.Max=65;
    Align(i).Tiltseries.Min=-65;
    Align(i).Tiltseries.Angles=[-65:5:65];
    Align(i).Size.X=64;
    Align(i).Size.Y=64;
    Align(i).Size.Z=64;
end;

