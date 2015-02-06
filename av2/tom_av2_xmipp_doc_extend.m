function new_doc=tom_av2_xmipp_doc_extend(Inputdoc,num_of_frames,find_what,replace_with,add,f_docOut)
%  tom_av2_xmipp_doc_extend extends a doc-file for sumed images 2 frames
%  
%     doc=tom_av2_xmipp_doc_extend(Inputdoc,num_of_frames,find_what,replace_with,add,f_docOut)
%  
%  PARAMETERS
%  
%    INPUT
%     Inputdoc                 input picklist
%     num_of_frames            number of frames 4 extension   
%     find_what                part 2 be replaced 4 the frame name
%     replace_with             part 2 be replaced 4 the frame name
%     add                      added part usually extension (.spi)
%     f_docOut                 (opt.) filename 4 output of the picklist  
%     
%    OUTPUT
%     doc             extended pickList
% 
%  EXAMPLE
%  
%   new_doc=tom_av2_xmipp_doc_extend();
%
%  REFERENCES
%  
%  SEE ALSO
%     ...
%  
%     created by FB 25/04/13
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom

if (nargin < 6)
    f_docOut='';
end;

if (isstruct(Inputdoc)==0)
    Inputdoc=tom_xmippdocread(Inputdoc);
end;

%alloc memory
new_doc=[];
for ii=1:num_of_frames
    new_doc=cat(1,new_doc,Inputdoc);
end;

zz=1;
for i=1:length(Inputdoc)
    for ii=1:num_of_frames
        new_doc(zz)=Inputdoc(i);
        [a b c]=fileparts(Inputdoc(i).name);
        tmp_name=Inputdoc(i).name;
        for iii=1:length(find_what)
            tmp_name=[strrep(tmp_name,find_what{iii},replace_with{iii}) ];
        end;
        new_doc(zz).name=[tmp_name '' b  '_F_' num2str(ii) add];
        zz=zz+1;
    end;
end;

if (isempty(f_docOut)==0)
    tom_xmippdocwrite(f_docOut,new_doc);
end;




