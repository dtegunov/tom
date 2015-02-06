function vol=tom_torus(sz_vol,diameter,thickness,decay)

%TOM_TORUS creates a 3D volume containing a torus
%
%    vol=tom_torus(sz_vol,diameter,thickness,decay)
%
%PARAMETERS
%
%  INPUT
%   sz_vol              volume size
%   diameter            torus diameter 
%   thickness           torus thickness 
%   decay               torus decay, smoother 
%  
%  OUTPUT
%   vol                 torus volume
%
%EXAMPLE
%   torus=tom_torus([128 128 128],35,5,2);
%
%REFERENCES
%
%SEE ALSO
%   TOM_SPHERE, TOM_CIRCLE
%
%   created by CH, SN 10/29/09
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

sphere=tom_spheremask(ones(sz_vol),thickness,decay,[sz_vol(1)./2+1 sz_vol(2)./2+1 sz_vol(3)./2+1]);

sphere_shift=tom_shift(sphere,[diameter./2+1 0 0]);

vol=tom_symref(sphere_shift,100);