function doc_out=tom_av2_xmipp_doc2angcldoc(input_doc,filename_out,class_base,class_ext)
%TOM_AV2_XMIPP_DOC2ANGCLDOC counts the particles per class and returns a
%                           doc with the information
%
%   doc_out=tom_av2_xmipp_doc2angcldoc(doc)
%
%PARAMETERS
%
%  INPUT
%   input_doc           name of doc or struct
%   filename_out        (opt.) outputfilename
%   class_base          ('class_') basename for classes
%   class_ext           ('.spi')   extension for classes   
%  
%  OUTPUT
%   doc_out            class doc file in memory
%
%EXAMPLE
%
%  doc_st=tom_av2_xmipp_doc2angcldoc('Iter_12_current_angles_1_cut.doc','test.doc');
%  
%REFERENCES
%
%SEE ALSO
%   tom_av2_xmipp_angular_plot
%
%   created by FB 10/04/2012
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
   filename_out='';
end;

if (nargin < 3)
   class_base='class_';
end;

if (nargin < 4)
   class_ext='.spi';
end;

if (isstruct(input_doc)==0)
    input_doc=tom_xmippdocread(input_doc);
end;

cl_vect=[input_doc(:).ref];
cl_vect_u=unique(cl_vect);

doc_out=input_doc(1:length(cl_vect_u));
doc_out=rmfield(doc_out,'maxCC');

doc_out(1).header='; Headerinfo columns: rot (1) , tilt (2), psi (3), Xoff (4), Yoff (5), Weight (6), Flip (7)';
for i=1:length(cl_vect_u)
   part_idx=find(cl_vect==cl_vect_u(i));
   doc_out(i).name=[class_base num2str(cl_vect_u(i)) class_ext];
   doc_out(i).rot=input_doc(part_idx(1)).rot;
   doc_out(i).part_idx=cl_vect_u(i);
   doc_out(i).tilt=input_doc(part_idx(1)).tilt;
   doc_out(i).weight=length(part_idx);
   doc_out(i).cols=7;
   doc_out(i).psi=0;
   doc_out(i).yoff=0;
   doc_out(i).xoff=0;
   doc_out(i).ref=input_doc(part_idx(1)).ref;
   doc_out(i).run_num=i; 
end;


if (isempty(filename_out)==0)
    doc_out=rmfield(doc_out,'ref');
    tom_xmippdocwrite(filename_out,doc_out);
end;

