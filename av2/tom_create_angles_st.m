function angles=tom_create_angles_st(tiltaxis,tiltangle,scheme)
%TOM_CREATE_ANGLES_ST computes an quasi-equal angular sampling
%   equal, spider- and xmipp-style.
%
%   angles=tom_create_angles_st(tiltaxis,tiltangle,scheme)
%
%PARAMETERS
%
%  INPUT
%   tiltaxis            ...
%   tiltangle           ...
%   scheme              'equal','spider', 'xmipp'
%  
%  OUTPUT
%   angles              equally spaced angles
%
%EXAMPLE
%   angles=tom_create_angles_st([0:10:90],[0:10:90],'xmipp');
%   
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_equal_angular_spacing
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

if (strcmp(scheme,'equal'))
    z=1;
    if (tiltaxis(1)==0 && tiltaxis(2)==0 && tiltaxis(3)==0)
        tiltaxis(1)=0;
        tiltaxis(2)=1;
        tiltaxis(3)=0;
    end;
    for j=tiltaxis(1):tiltaxis(2):tiltaxis(3)
        for i=tiltangle(1):tiltangle(2):tiltangle(3)
            angles(1,z)=j;
            angles(2,z)=i;
            z=z+1;
        end
    end;
    if (exist('angles','var')==0)
        angles='';
    end;

end;


if (strcmp(scheme,'spider') || strcmp(scheme,'xmipp'))
     angles=tom_av2_equal_angular_spacing([tiltaxis(1) tiltaxis(3)],[tiltangle(1) tiltangle(3)],tiltaxis(2),scheme);
     angles=angles';
     tmp=angles; 
     tmp(2,:)=tmp(2,:).*-1;
     angles=[angles tmp(:,2:end)];
end;











