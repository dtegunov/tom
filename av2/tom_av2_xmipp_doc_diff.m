function diff_doc=tom_av2_xmipp_doc_diff(doc1,doc2,outputdoc)
%tom_av2_xmipp_doc_diff subtracts 2 doc files
%
%   diff_doc=tom_av2_xmipp_doc_diff(doc1,doc2,outputdoc)
%
%PARAMETERS
%
%  INPUT
%   doc1                  input xmipp doc-file (filename or struct)
%   doc2                  input xmipp doc-file (filename or struct) 
%   outputdoc             (opt.) output doc-filename
%
%  
%  OUTPUT
%    diff_doc             struct in mem        
%  
%
%EXAMPLE
%   
% diff_doc=tom_av2_xmipp_doc_diff('Iter_1_current_angles.doc','Iter_2_current_angles.doc');
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


if (isstruct(doc1)==0)
    doc1=tom_xmippdocread(doc1);
end;


if (isstruct(doc2)==0)
   doc2=tom_xmippdocread(doc2);
end;

all_fields=fieldnames(doc1);
zz=1;
for i=1:length(all_fields)
    if ((isnumeric(getfield(doc1(1),all_fields{i}))) && (strcmp(all_fields{i},'part_idx_unique')==0) &&   (strcmp(all_fields{i},'cols')==0) && (strcmp(all_fields{i},'part_idx')==0) && (strcmp(all_fields{i},'run_num')==0) )
        numeric_fields{zz}=all_fields{i};
        zz=zz+1;
    end;
end;

diff_doc=doc1;

for i=1:length(numeric_fields)
    for ii=1:length(diff_doc)
        diff_doc(ii)=setfield(diff_doc(ii),numeric_fields{i},(getfield(doc1(ii),numeric_fields{i})-getfield(doc2(ii),numeric_fields{i})) );
    end;
end;



if (nargin>2)
    tom_xmippdocwrite(outputdoc,diff_doc);
end;



