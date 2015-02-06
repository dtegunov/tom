function tom_av3_switchfilenames(alignfile)
%TOM_AV3_SWITCHFILENAMES recovers the original filenames that were ...
% temporarily replaced during auto refinement runs
%
%   tom_amira_createisosurface(filelabel,threshold, label, color, transformmatrix, iconposition, host)
%
%PARAMETERS
%
%  INPUT
%   alignfile           filename of the input alignment file
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av3_switchfilenames(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 01/20/06
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



try
    a = load(alignfile);
    Align = a.Align;
catch
    error('could not open alignment file');
end

j = size(Align,1);
for k = 1:size(Align,2)
    tmp = Align(j,k).Filename;
    Align(j,k).Filename = Align(j,k).Tempfilename;
    Align(j,k).Tempfilename = tmp;
end

try
    save(alignfile);
catch
    error('could not save alignment file, please check file permissions!');
end