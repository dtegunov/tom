function volsym = tom_symref(vol,nfold,flag,angle)
%TOM_SYMREF does n-fold symmetrization of a 3D reference
%
%   volsym = tom_symref(vol,nfold,flag,angle)
%
%   If a volume VOL is assumed to have a n-fold symmtry axis along z it can
%   be rotationally symmetrized using TOM_SYMREF
%
%PARAMETERS
%
%  INPUT
%   vol                 3D volume to be symmterized
%   nfold               rotational symmetry along z (>= 2)
%   flag                for a special case of C2 symmetry with a rotation
%   angle               angle for the special case
%  
%  OUTPUT
%   volsym              symmetrized volume
%
%EXAMPLE
%   ... = tom_symref(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FF 07/28/03
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

if(nfold > 0)
    if (nargin==2)
        flag='normal';
    end;

    if (strcmp(flag,'y-axis')==1)

        iangle = 360/nfold;
        volsym=vol;
        nphi = 360/nfold;
        for ind = 2:nfold
            phi = nphi*(ind-1);
            volsym = volsym + double(tom_rotate(vol,[270 90 phi]));
        end;
        volsym = volsym/nfold;
    end;
    
   if (strcmp(flag,'normal')==1)

        iangle = 360/nfold;
        volsym=vol;
        nphi = 360/nfold;
        for ind = 2:nfold
            phi = nphi*(ind-1);
            volsym = volsym + double(tom_rotate(vol,[phi 0 0]));
        end;
        volsym = volsym/nfold;
    end;
    
    if (strcmp(flag,'special')==1)
        vol_rot=tom_rotate(vol,angle);
        volsym=vol+vol_rot;
    end
else
    volsym = vol;
end;


