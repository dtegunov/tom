function tom_makefun
%TOM_MAKEFUN compiles C files
%
%   tom_makefun
%
%   Compiles C file. Create .dll on windows, .mexglx under linux
%
%PARAMETERS
%
%  INPUT
%  
%  OUTPUT
%
%EXAMPLE
%   tom_makefun
%
%REFERENCES
%
%SEE ALSO
%   MEX
%
%   created by FF 05/12/04
%   updated by AK 03/08/06
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

cd IOfun/
mex tom_emwriteinc.c  
%the next file does not work correctly with win64 anymore!!!
%mex tom_emreadinc.c

if strcmp(computer,'PCWIN')
    mex -DWIN32 tom_emreadinc_resample2.c
    mex -DWIN32 tom_emreadinc_resample3.c
    mex -DWIN32 tom_emwriteinc2.c
elseif strcmp(computer,'PCWIN64')
    mex -DWIN64 tom_emreadinc_resample2.c
    mex -DWIN64 tom_emreadinc_resample3.c
    mex -DWIN64 tom_emwriteinc2.c
elseif strcmp(computer,'MACI')
    mex tom_emreadinc_resample2.c
    mex tom_emreadinc_resample3.c
    mex -DMACOSX tom_emwriteinc2.c
else
    mex tom_emreadinc_resample2.c
    mex tom_emreadinc_resample3.c
    mex tom_emwriteinc2.c
end
cd ..
cd Reconstruction/
mex tom_backproj3dc.c 
mex tom_dist.c
cd ..
cd  Sptrans/
mex tom_rotatec.c 
mex tom_bininc.c
mex tom_rotatec.c
cd ..
cd Analysis/
mex tom_peakinc.c
cd ..
cd Filtrans
mex aniso3.c
mex tom_calc_weight_functioninc.c
cd ..
