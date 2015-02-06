function tom_xmippdocrotate(input_doc,output_doc,angle)
%tom_xmippdocrotate rotates a doc file by adding a angle 
%   to all particles in the doc file
%
%  tom_xmippdocrotate(input_doc,output_doc,angle)
%
%  TOM_XMIPPDOCROTATE rotates a doc file
%                  
%                  
%PARAMETERS
%
%  INPUT
%   inputdoc              *.doc filename 
%   output_doc            name of output em stack
%   angle                 rotation angle in zxz (tom-convention) !!!!
%
%  OUTPUT
%
%EXAMPLE
%     
%  tom_xmippdocrotate('test_compl.doc','test_compl_rot.doc',[90 90 0]);
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


doc=tom_xmippdocread(input_doc);

for i=1:length(doc)
    [rr,curr_ang_tom] = tom_eulerconvert_xmipp(doc(i).rot,doc(i).tilt,doc(i).psi);
    new_ang_tom=tom_sum_rotation([angle;curr_ang_tom],[0 0 0 ;0 0 0]);
    [rotmatrix,new_ang] = tom_eulerconvert_xmipp(new_ang_tom(1),new_ang_tom(2),new_ang_tom(3),'tom2xmipp'); 
    doc(i).rot=new_ang(1);
    doc(i).tilt=new_ang(2);
    doc(i).psi=new_ang(3);
end;

tom_xmippdocwrite(output_doc,doc);


