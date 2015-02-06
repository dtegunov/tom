function cc=tom_av2_xmipp_cmp_ref2class(doc_name,mask,mask2)
%TOM_AV2_XMIPP_CMP_REF2CLASS correlates ang classes angainst ref lib
%
%   cc_list=tom_av2_xmipp_cmp_ref2class(doc_name)
%
%PARAMETERS
%
%  INPUT
%   doc_name                 reconstruction.doc
%   mask                     mask applied 2 images  
%   mask2                      
%   
%EXAMPLE
%   
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB 05/14/12
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


%parse inputs


find_what='/ProjMatchClasses/proj_match_class';
repl_with='/ReferenceLibrary/ref';

doc=tom_xmippdocread(doc_name);

cc=zeros(length(doc),2);
%for i=1:length(doc)
for i=635:635
    cl=tom_spiderread(doc(i).name);
    ref=tom_spiderread(strrep(doc(i).name,find_what,repl_with));
    cc(i,1)=tom_get_max_numeric({doc(i).name});
    cc(i,2)=tom_ccc(cl.Value.*mask,ref.Value.*mask,'norm');
    cc(i,3)=tom_ccc(cl.Value.*mask2,ref.Value.*mask2,'norm');
    if (cc(i,1)==636)
        disp(' ');
    end;
    if (mod(i,100)==0)
        disp(num2str(i));
    end;
end;

cc=cc(1:i,:);











