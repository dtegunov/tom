function pl_out=tom_av2_extr_img_pickList(align2d,img_name,flag)
%TOM_AV2_EXTR_IMG_PICKLIST reduces picklist 2 one image
%   pl_out=tom_av2_extr_img_pickList(align2d,img_name)
%
%PARAMETERS
%
%  INPUT
%   align2d             align2d struct
%   img_name            name of the image
%   flag                (full) or filename
%
%  OUTPUT
%    pl_out             reduced picklist 
%
%
%EXAMPLE
%
%  pl_red=tom_av2_extr_img_pickList(align2d,'/fs/pool/pool-engelhardt/EG2/101109_dia_corr/high/dia_46.em');
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

if (nargin < 3)
    flag='full';
end;

for i=1:size(align2d,2)
    if (strcmp(flag,'full'))
        f_names{i}=align2d(1,i).filename;
    else
        [a b c]=fileparts(align2d(1,i).filename);
        f_names{i}=[b c];
    end;
end;

idx=find(ismember(f_names,img_name));

pl_out=align2d(1,idx);