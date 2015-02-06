function tom_av2_picklist2filefilter(pickList,filefilter_name)
%TOM_AV2_PICKLIST2FILEFILTER generates a filefilter according 2 picklist
%                            only images with picked particles are in the filefilter  
%
%
%   tom_av2_picklist2filefilter(pickList,filefilter_name)
%
%PARAMETERS
%
%  INPUT
%   pickList           filename or in memory  
%   filefilter_name    output filename for filefilter                   
%  
%  OUTPUT
%      
%
%EXAMPLE
%  
%  tom_av2_picklist2filefilter('al__fs_pool_pool-engelhardt_EG2_101109_dia_corr_high_corr.mat','al__fs_pool_pool-engelhardt_EG2_101109_dia_corr_high_corr_ff.mat')
%  
%REFERENCES
%
%SEE ALSO
%   
% tom_av2_split_picklist
%
%   created by fb 01/08/07
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


if (isstruct(pickList)==0)
    load(pickList);
    pickList=align2d;
end;



for i=1:size(pickList,2)
    [a b c]=fileparts(pickList(1,i).filename);
    all_filenames{i}=[b c];
end;

u_filenames=unique(all_filenames);


d=dir([a '/*' c] );
particlepicker.filelist={d.name};

for i=1:length(particlepicker.filelist)
     if (ismember(particlepicker.filelist{i},u_filenames))
        particlepicker.filefilter{i}=[1];
    else
        particlepicker.filefilter{i}=[0];
    end;
end;


save(filefilter_name,'particlepicker');




