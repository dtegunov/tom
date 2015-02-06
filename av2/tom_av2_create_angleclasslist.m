function tom_av2_create_angleclasslist(align2d,iter_num,datname)
%TOM_ANGLECLASSLIST creates ...
%
%   tom_av2_create_angleclasslist(align2d,iter_num,datname)
%
%PARAMETERS
%
%  INPUT
%   align2d             ...
%   iter_num            ...
%   datname             ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av2_create_angleclasslist(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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


 fid = fopen(datname,'W');
 fprintf(fid,'Stackname    %s Number of Particles  %d  \n',align2d(iter_num,1).filename,size(align2d(iter_num,:),2));
 fclose(fid);


for i=1:size(align2d(iter_num,:),2)
    proj_nr=align2d(iter_num,i).angleclass.proj_nr;
    angle_rot=align2d(iter_num,1).angular_scan(2,proj_nr);
    angle_nut=align2d(iter_num,1).angular_scan(1,proj_nr);
    fid = fopen(datname,'a');
    fprintf(fid,'part nr: %d   angle1: %4.2f   angle2 %4.2f \n',i,angle_rot,angle_nut);
    fclose(fid);
end;