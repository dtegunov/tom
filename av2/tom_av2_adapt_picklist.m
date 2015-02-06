function align2d=tom_av2_adapt_picklist(f_org_align,f_adapt_align,fieldname,param1,param2)
%TOM_AV2_ADAPT_PICKLIST changes all fields in picklist
%
%   [angle_out shift_out ccc aligned_part_sum]=tom_av2_align(ref,im,mask_im,mask_cc_rot,mask_cc_trans,filter_st,num_of_iterations,demo)
%
%PARAMETERS
%
%  INPUT
%   f_org_align        reference image
%   f_adapt_align      image to be aligned
%   fieldname          fieldname to be changed in list
%   param1             1 parameter for changing    
%   param2             2 parameter for changing
%   
%  OUTPUT
%   pl_adapt            adapted picklist
%
%
%EXAMPLE
%
% %1 change root path of images
% tom_av2_adapt_picklist('pick_high.mat','pick_high_adapt.mat','filename','/fs/sandy02/lv04/pool/pool-titan1/CWT/data/06__10042011/high/','/fs/pool/pool-titan1/CWT/data/06__10042011_3/high/');
%
% %2 change radius of picklist
% tom_av2_adapt_picklist('pick_high.mat','pick_high_adapt.mat','radius',256);
% 
% %3 adapt low 2 low_corr
% tom_av2_adapt_picklist('fr1113final.mat','pick_high_corr.mat','filename','/low/','/low_corr/');
%
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

disp(['loading ' f_org_align]);
in=load(f_org_align);
disp('done !');

align2d=in.align2d;

if (isnumeric(param1))
    param1=num2str(param1);
else
    if (strcmp(fieldname,'filename')==0)
        param1=['''' param1 ''''];
    end;
end;


disp('adapting ');
for i=1:size(in.align2d,2)
    if (mod(i,10000)==0) 
        disp([ num2str(i) ' parts processed']);
    end;

    if (strcmp(fieldname,'filename'))
        old_name=align2d(1,i).filename;
        align2d(1,i).filename=strrep(align2d(1,i).filename,param1,param2); 
        if (strcmp(old_name,align2d(1,i).filename))
            disp(['Warning: filename: ' align2d(1,i).filename ' not changed ']);
        end;
        continue;
    end;
    eval(['align2d(1,i).' fieldname '='  param1 ';']);
    
end;
disp('done!');


if (isempty(f_adapt_align)==0)
    disp(['saving ' f_adapt_align]);
    save(f_adapt_align,'align2d');
    disp('done !');
end;


