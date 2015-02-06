function tom_build_3d_phantom_stack(filename_out,number)
%TOM_BUILD_3D_PHANTOM_STACK creates ...
%
%   tom_build_3d_phantom_stack(filename_out,number)
%
%PARAMETERS
%
%  INPUT
%   filename_out        ...
%   number              ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_build_3d_phantom_stack(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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

error(nargchk(0, 2, nargin, 'struct'))




moped=tom_av2_build_artificial26S;


ang_inkre=360./number(1);

angle=0;
zz=1;
for i=1:number(1)
    
    for ii=1:number(2)
        moped_rot=tom_rotate(moped,[0 0 angle]);
        moped_rot=moped_rot+rand(size(moped_rot));
        
        yyy=ones(size(moped_rot));
        wedge=tom_wedge(yyy,30);
        moped_rot=tom_apply_weight_function(moped_rot,wedge);
        tom_emwrite([filename_out num2str(zz) '.em'],moped_rot);
        zz=zz+1;
    end;
    angle=angle+ang_inkre;
end;

