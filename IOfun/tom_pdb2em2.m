function emmap = tom_pdb2em2(pdbdata, pixelsize, dim,center)
%TOM_PDB2EM creates ...
%
%   emmap = tom_pdb2em(pdbdata, pixelsize, dim)
%
%PARAMETERS
%
%  INPUT
%   PDBDATA             expected as structure as derived from TOM_PDBREAD
%   PIXELSIZE           desired pixelsize in Angstrom
%   DIM                 dimension of desired cube
%   CENTER              (1) pre center struct  
%
%  OUTPUT
%   emmap               ...
%
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FF 2003
%   updated by RK 2006 bugfixed - now centering the protein model prior to
%   mapping
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

if (nargin < 4)
    center=1;
end;

% mod rk 060428
minx=min([pdbdata.Model.Atom.X]);
miny=min([pdbdata.Model.Atom.Y]);
minz=min([pdbdata.Model.Atom.Z]);
maxx=max([pdbdata.Model.Atom.X]);
maxy=max([pdbdata.Model.Atom.Y]);
maxz=max([pdbdata.Model.Atom.Z]);

if (center==1)
    for i=1:size(pdbdata.Model.Atom,2)
        pdbdata.Model.Atom(i).X = pdbdata.Model.Atom(i).X - ((maxx + minx) / 2);
        pdbdata.Model.Atom(i).Y = pdbdata.Model.Atom(i).Y - ((maxy + miny) / 2);
        pdbdata.Model.Atom(i).Z = pdbdata.Model.Atom(i).Z - ((maxz + minz) / 2);
    end;
end;
%---

emmap = zeros(dim,dim,dim);
x = round([pdbdata.Model.Atom.X]/pixelsize)+ floor(dim/2); 
y = round([pdbdata.Model.Atom.Y]/pixelsize)+ floor(dim/2);
z = round([pdbdata.Model.Atom.Z]/pixelsize)+ floor(dim/2);
%atom = reshape([pdbdata.Model.Atom.AtomName],4,size([pdbdata.Model.Atom.AtomName],2)/4);
for iatom = 1:size(x,2)
    if ((x(iatom) > 0) &&  (y(iatom) > 0) &&  (z(iatom) > 0))
        if (strmatch(pdbdata.Model.Atom(iatom).AtomName(1),'C') )  
            emmap(x(iatom),y(iatom),z(iatom))=emmap(x(iatom),y(iatom),z(iatom))+6;
        elseif (strmatch(pdbdata.Model.Atom(iatom).AtomName(1),'N')) 
            emmap(x(iatom),y(iatom),z(iatom))=emmap(x(iatom),y(iatom),z(iatom))+7;
        elseif (strmatch(pdbdata.Model.Atom(iatom).AtomName(1),'O')) 
            emmap(x(iatom),y(iatom),z(iatom))=emmap(x(iatom),y(iatom),z(iatom))+8;
        elseif (strmatch(pdbdata.Model.Atom(iatom).AtomName(1),'H')) 
            emmap(x(iatom),y(iatom),z(iatom))=emmap(x(iatom),y(iatom),z(iatom))+1;
        elseif (strmatch(pdbdata.Model.Atom(iatom).AtomName(1),'P')) 
            emmap(x(iatom),y(iatom),z(iatom))=emmap(x(iatom),y(iatom),z(iatom))+15;
        elseif (strmatch(pdbdata.Model.Atom(iatom).AtomName(1),'S')) 
            emmap(x(iatom),y(iatom),z(iatom))=emmap(x(iatom),y(iatom),z(iatom))+16;
        end;
    end;
end;