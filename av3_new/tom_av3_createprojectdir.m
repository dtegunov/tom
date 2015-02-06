function tom_av3_createprojectdir(directory)
%TOM_AV3_CREATEPROJECTDIR creates projectfolder for particle averaging ...
%and alignment
%
%   tom_av3_createprojectdir(directory)
%
%PARAMETERS
%
%  INPUT
%   directory           the target directory to generate
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av3_createprojectdir(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV3_OSCARGUI
%
%   created by AK 10/20/05
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


mkdir(directory);
cd(directory);
mkdir 'masks';
mkdir 'vols';
mkdir 'settings';
mkdir 'scripts';
mkdir 'aligns';
mkdir 'psfs';
mkdir 'templates';
mkdir 'particles';
mkdir 'auto_refinements';
mkdir 'ccfs';
mkdir 'logs';
mkdir 'picklists';

