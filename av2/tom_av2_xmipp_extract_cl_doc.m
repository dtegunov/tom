function tom_av2_xmipp_extract_cl_doc(f_doc,f_outputdoc,cl_vect,inv_cl_vect,f_output_sel)
%TOM_AV2_XMIPP_EXTRACT_CL_DOC filters xmipp doc file
%   
%
%  TOM_AV2_XMIPP_EXTRACT_CL_DOC(f_doc,f_outputdoc,cl_vect,f_output_sel)
%
%  TOM_AV2_XMIPP_EXTRACT_CL_DOC filters xmipp doc file according 2 given
%  classes
%  
%
%PARAMETERS
%
%  INPUT
%   f_doc               *.doc filename use abs filename 
%   f_outputdoc         name of output em stack
%   cl_vect             classes vector
%   inv_cl_vect         1=use all classes in cl vect 0 use all classes not in cl_vect  
%   f_output_sel        (opt) sel outputname
%  
% OUTPUT
%
%EXAMPLE
%     
% tom_av2_xmipp_extract_cl_doc();
% 
%  
%
%REFERENCES
%
%NOTE 
%  
%
%SEE ALSO
%   tom_av2_em_classify3d
%
%   created by FB 08/09/09
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
    doc=tom_xmippdocread(f_doc);
catch ME
    disp(['Error reading: ' f_doc]);
    error(ME.message);
end;


idx=zeros(length(doc),1);
%alloc memory 
for i=1:size(doc,1)
    tmp=ismember(doc(i).ref,cl_vect);
    if (tmp > 0)
        idx(i)=1;
    end;
end;

if (inv_cl_vect==1)
    filt_idx=find(idx==0);
else
    filt_idx=find(idx==1);
end;
    
disp(['filt. doc has length of: ' num2str(length(filt_idx))]);

new_doc=doc(filt_idx);
new_doc(1).header=doc(1).header;
new_doc(1).part_idx_unique=doc(1).part_idx_unique;
tom_xmippdocwrite(f_outputdoc,new_doc);

fid=fopen(f_output_sel,'wt');
for i=1:size(new_doc,1)
    fprintf(fid,'%s 1  \n',new_doc(i).name);
end;
fclose(fid);



