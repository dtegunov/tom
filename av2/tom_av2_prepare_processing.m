function tom_av2_prepare_processing(destination_directory,flag)

%TOM_AV2_PREPARE_PROCESSING sets up the folder structure for a new
%experiment.
%
%   tom_av2_prepare_processing(destination_directory,flag)
%
%PARAMETERS
%
%  INPUT
%   flag                    'pairs' or 'singles', depending on acquisition scheme, default is singles
%   destination_directory   complete path of destination folder for experiment. must exist.
%
%
%EXAMPLE
%   %for focal pairs:
%   tom_av2_prepare_processing('/fs/pool/pool-nickell/26S/em/data/eri/2d/100223_sc26s_rpn1gfp/','pairs');
%   %for single shots:
%   tom_av2_prepare_processing('/fs/pool/pool-nickell/26S/em/data/eri/2d/100223_sc26s_rpn1gfp/');
%REFERENCES
%
%SEE ALSO
%
%   created by SB 10/02/24
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

if nargin<2
    flag='';
end;

unix(['mkdir ' destination_directory '/log']);
unix(['mkdir ' destination_directory '/log/pick']);
unix(['mkdir ' destination_directory '/log/rec']);

unix(['mkdir ' destination_directory '/low']);
unix(['mkdir ' destination_directory '/lowrefit']);
unix(['mkdir ' destination_directory '/log/rec/parts_low_128']);
unix(['mkdir ' destination_directory '/log/rec/parts_low_64']);
unix(['mkdir ' destination_directory '/low_del']);
unix(['mkdir ' destination_directory '/low_corr']);

switch flag;

    case 'pairs'

        unix(['mkdir ' destination_directory '/high']);
        unix(['mkdir ' destination_directory '/highrefit']);
        unix(['mkdir ' destination_directory '/log/rec/parts_high_128']);
        unix(['mkdir ' destination_directory '/log/rec/parts_high_64']);
        unix(['mkdir ' destination_directory '/high_del']);
        unix(['mkdir ' destination_directory '/high_corr']);
        
end;