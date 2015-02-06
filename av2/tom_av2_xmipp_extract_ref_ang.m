function tom_av2_xmipp_extract_ref_ang(f_input_doc,f_output)
%tom_av2_xmipp_extract_ref_ang extracts ref_angles.doc from doc file
%  
% tom_av2_xmipp_extract_ref_ang(f_input_doc,f_output)
%
%  TOM_AV2_XMIPP_ALIGN_STACK extracts ref_angles.doc from doc file
%  
%
%PARAMETERS
%
%  INPUT
%   f_input_doc         *.doc filename 
%   f_output              name of output ref_angles file (need for xmipp 3d rec and debug!)
%
%  OUTPUT
%
%EXAMPLE
%     
%  
%
%
%REFERENCES
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


doc=tom_xmippdocread(f_input_doc);

refs=unique([doc(:).ref]);

angs=zeros(length(refs),3);

for i=1:length(refs)
    idx=find([doc(:).ref]==refs(i));
    angs(i,:)=[doc(idx(1)).ref doc(idx(1)).rot doc(idx(1)).tilt];
    if (mod(i,1000)==0)
        disp(num2str(i));
    end;
end;



base_out_st='%5d %d  %10.5f %10.5f %10.5f\n';
fid=fopen(f_output,'wt');
%for i=1:length(refs)
for i=1:max(angs(:,1)) 
    idx=find(angs(:,1)==i);
    if (isempty(idx)==0)
        fprintf(fid,base_out_st,i,3,angs(idx(1),2),angs(idx(1),3),0);
    else
        fprintf(fid,base_out_st,i,3,0,0,0);
    end;    
end;
fclose(fid);

