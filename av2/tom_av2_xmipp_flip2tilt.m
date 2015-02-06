function new_doc=tom_av2_xmipp_flip2tilt(doc_in,doc_out)
%TOM_AV2_XMIPP_FLIP2TILT calculates projection angles on whole Euler sphere
%from data with flip=0 and 1
%
%   new_doc=tom_av2_xmipp_flip2tilt(doc_in,doc_out)
%
%PARAMETERS
%
%  INPUT
%   doc_in              input doc file name (with flip 0 and 1)
%  
%  OUTPUT
%   doc_out             output docfile name 
%
%EXAMPLE
%
%   new_doc=tom_av2_xmipp_flip2tilt('parts.doc','parts-noflip.doc')  
%   -> this gives the angular distribution of the whole Euler sphere and
%   writes the resultion doc file parameters to parts-noflip.doc and to the
%   structure new_doc written to the memory
%
%   new_doc=tom_av2_xmipp_flip2tilt('parts.doc')
%   -> the new doc file parameters are not written to a doc file but can be
%   found in the structure new_doc written to the memory
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

if nargin<2
    doc_out='';
end;

doc_name=doc_in;
doc_name_out=doc_out;
doc=tom_xmippdocread(doc_name);

trans_angle=[0 180 0]; %zyz
[rotM trans_angle_zxz]=tom_eulerconvert_xmipp(trans_angle(1),trans_angle(2),trans_angle(3));

%alloc memory
new_doc=doc;

for i=1:length(doc)
    ang=[doc(i).rot doc(i).tilt doc(i).psi];
    if doc(i).flip==1
        [rotM ang_zxz]=tom_eulerconvert_xmipp(ang(1),ang(2),ang(3));
        [euler_sum_zxz shift_out rott]=tom_sum_rotation([ang_zxz; trans_angle_zxz],[0 0 0;0 0 0]);
        [rotM ang_sum_zyz]=tom_eulerconvert_xmipp(euler_sum_zxz(1),euler_sum_zxz(2),euler_sum_zxz(3),'tom2xmipp');
        new_doc(i).rot=ang_sum_zyz(1);
        new_doc(i).tilt=ang_sum_zyz(2);
        new_doc(i).psi=ang_sum_zyz(3);       
        new_doc(i).yoff=doc(i).yoff;
        new_doc(i).xoff=doc(i).xoff;
        new_doc(i).flip=0;
    end;
end;
if isempty(doc_out)==0
    tom_xmippdocwrite(doc_name_out,new_doc);
end;