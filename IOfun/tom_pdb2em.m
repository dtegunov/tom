function emmap = tom_pdb2em(pdbdata, pixelsize, dim, center)
%TOM_PDB2EM creates a density map based on pdb-structure
%   from the pdb database: www.pdb.org.
%
%   emmap = tom_pdb2em(pdbdata, pixelsize, dim)
%
%PARAMETERS
%
%  INPUT
%   PDBDATA             expected as structure as derived from TOM_PDBREAD
%                       or pdbread
%   PIXELSIZE           desired pixelsize in Angstrom
%   DIM                 dimension of desired cube
%   CENTER              If 'RK', use Roland's dodgy stuff
%                       otherwise: 'cent': center according to center of
%                       mass, 'none': leave it as it is
%  
%  OUTPUT
%   emmap               volume
%
%EXAMPLE
%   read the Thermoplasma 20S proteasome from the database.
%   proteasome = pdbread('http://www.rcsb.org/pdb/files/1pma.pdb');
%   create a EM density map of 96x96x96 voxel and a objectpixelsize
%   of 2.1 A.
%   emmap = tom_pdb2em(proteasome, 2.1, 96);
%
%REFERENCES
%
%SEE ALSO
%   tom_pdb2em2, www.pdb.org
%
%   created by FF 2003
%   updated by RK 2006 bugfixed - now centering the protein model prior to
%                       mapping
%   update 11/2010 FF: dodgy RK centering now optional (center='RK')
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

% mod rk 060428
if exist('center')
    if strcmp(center, 'RK')
        minx=min([pdbdata.ATOM.X]);
        miny=min([pdbdata.ATOM.Y]);
        minz=min([pdbdata.ATOM.Z]);
        maxx=max([pdbdata.ATOM.X]);
        maxy=max([pdbdata.ATOM.Y]);
        maxz=max([pdbdata.ATOM.Z]);
        for i=1:size(pdbdata.ATOM,2)
            pdbdata.ATOM(i).X = pdbdata.ATOM(i).X - ((maxx + minx) / 2);
            pdbdata.ATOM(i).Y = pdbdata.ATOM(i).Y - ((maxy + miny) / 2);
            pdbdata.ATOM(i).Z = pdbdata.ATOM(i).Z - ((maxz + minz) / 2);
        end;
    elseif strcmp(center, 'cent')
        cx = sum([pdbdata.ATOM.X])/size([pdbdata.ATOM.X],2);
        cy = sum([pdbdata.ATOM.Y])/size([pdbdata.ATOM.Y],2);
        cz = sum([pdbdata.ATOM.Z])/size([pdbdata.ATOM.Z],2);
        disp(['Center: cx=',num2str(cx),', cy=',num2str(cy),...
            ', cz=',num2str(cz)])
        for i=1:size(pdbdata.ATOM,2)
            pdbdata.ATOM(i).X = pdbdata.ATOM(i).X - cx;
            pdbdata.ATOM(i).Y = pdbdata.ATOM(i).Y - cy;
            pdbdata.ATOM(i).Z = pdbdata.ATOM(i).Z - cz;
        end;
    else
        disp('no re-centering of atoms')
    end;
end;
%---

emmap = zeros(dim,dim,dim);
x = round([pdbdata.ATOM.X]/pixelsize)+ floor(dim/2); 
y = round([pdbdata.ATOM.Y]/pixelsize)+ floor(dim/2);
z = round([pdbdata.ATOM.Z]/pixelsize)+ floor(dim/2);
ATOM = reshape([pdbdata.ATOM.AtomName],4,size([pdbdata.ATOM.AtomName],2)/4);
for iATOM = 1:size(x,2)
    if ((x(iATOM) > 0) &  (y(iATOM) > 0) &  (z(iATOM) > 0))
        if (strmatch(ATOM(2,iATOM),'C') )  
            emmap(x(iATOM),y(iATOM),z(iATOM))=emmap(x(iATOM),y(iATOM),z(iATOM))+6;
        elseif (strmatch(ATOM(2,iATOM),'N')) 
            emmap(x(iATOM),y(iATOM),z(iATOM))=emmap(x(iATOM),y(iATOM),z(iATOM))+7;
        elseif (strmatch(ATOM(2,iATOM),'O')) 
            emmap(x(iATOM),y(iATOM),z(iATOM))=emmap(x(iATOM),y(iATOM),z(iATOM))+8;
        elseif (strmatch(ATOM(2,iATOM),'H')) 
            emmap(x(iATOM),y(iATOM),z(iATOM))=emmap(x(iATOM),y(iATOM),z(iATOM))+1;
        elseif (strmatch(ATOM(2,iATOM),'P')) 
            emmap(x(iATOM),y(iATOM),z(iATOM))=emmap(x(iATOM),y(iATOM),z(iATOM))+15;
        elseif (strmatch(ATOM(2,iATOM),'S')) 
            emmap(x(iATOM),y(iATOM),z(iATOM))=emmap(x(iATOM),y(iATOM),z(iATOM))+16;
        end;
    end;
end;