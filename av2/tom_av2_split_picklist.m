function pl_out=tom_av2_split_picklist(align2d_in,flag,param,outpath)
%TOM_AV2_SPLIT_PICKLIST generates a seperate picklist for each data-set in
%the input picklist.
%
%    pl_out=tom_av2_split_picklist(align2d_in,flag,param,outpath)
%
%PARAMETERS
%
%  INPUT
%   align2d             tom align2d struct
%   flag                ('root_path')flag 2 split the picklist
%   param               ('')param for splitting the picklist
%   outpath             ('out') outpufolder for splitted picklists     
%
%
%EXAMPLE
%   
%    tom_av2_split_picklist(align2d,'root_path','','out');
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN/FB 01/24/06
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

if (nargin < 2)
    flag='root_path';
end;

if (nargin < 3)
    param='';
end;

if (nargin < 4)
    outpath='out';
end;

warning off; mkdir(outpath); warning on;

if (strcmp(flag,'root_path'))
    for i=1:size(align2d_in,2)
        all_f_root{i}=fileparts(align2d_in(1,i).filename);
    end;
    all_f_root_unique=unique(all_f_root);
    
    for i=1:length(all_f_root_unique)
        idx=find(ismember(all_f_root,all_f_root_unique{i}));
        pl_out{i}=align2d_in(1,idx);
        if (isempty(outpath)==0)
            align2d=pl_out{i};
            save([outpath '/al_' strrep(all_f_root_unique{i},'/','_')],'align2d');
            tom_av2_picklist2filefilter(align2d,[outpath '/al_' strrep(all_f_root_unique{i},'/','_') '_ff.mat']);
        end;
    end;
    
end;


  
  