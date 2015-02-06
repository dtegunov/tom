function s_doc=tom_av2_xmipp_sort_doc(input_doc,field,mode,f_output_doc)
%TOM_AV2_XMIPP_SORT_DOC sorts a xmipp .doc file
%
%   sorted_doc=tom_av2_xmipp_sort_doc(input_doc,f_output_doc,field,value)
%
%  
%
%PARAMETERS
%
%  INPUT
%   input_doc       filename of the input doc or doc struct in mem
%   field           fieldname doc file should be sorted 
%   mode            ('ascend') sorting direction ('ascend' or 'descend')
%   f_output_doc    ('') filename of sorted doc empty 2 switch off
%
%  OUTPUT
%   sorted_doc      sorted doc in memory
%
%EXAMPLE
%  
%  tom_av2_xmipp_sort_doc('../data/Iter_1_current_angles.doc','ref','ascend','../data/It1_sorted.doc');
%   
%
%REFERENCES
%
%SEE ALSO
%   
%   by fb 
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
    mode='ascend';
end;

if (nargin < 4)
    mode='ascend';
end;


if (isstruct(input_doc)==0)
    s_doc=tom_xmippdocread(input_doc);
else
    s_doc=input_doc; 
end;
head_tmp=s_doc(1).header;
p_idx=s_doc(1).part_idx_unique;


vect=eval(['[s_doc.' field ']' ]);

[s_val,idx]=sort(vect,2,mode); 

s_doc=s_doc(idx);

s_doc(1).header=head_tmp;
s_doc(1).part_idx_unique=p_idx;

if (isempty(f_output_doc)==0)
    tom_xmippdocwrite(f_output_doc,s_doc);
end;



