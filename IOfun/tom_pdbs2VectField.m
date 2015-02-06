function [chi_bild ptx pty ptz vx vy vz len pdb_head]=tom_pdbs2VectField(pdb_tail,pdb_head,pre_align,pixelsize,dim,verbose)
%tom_pdbs2VectField subtracts the coordinates of 2 pdbs(head - tail) and
%                   returns a matlab style vect field
%
%   [chi_bild px py pz vx vy vz]=tom_pdbs2vectVield(pdb_tail,pdb_head)
%
%PARAMETERS
%
%  INPUT
%
%   pdb_tail     pdb with tail positions (struct or filename)  
%   pdb_head     pdb with head positions (struct or filename) 
%   pre_align    (0) aligns head coordinates in respect 2 tail b4 subtraction
%   pixelsize    (1) pixelsize in Ang 
%   dim          (0) new dimension
%   verbose      (1) verbose flag
%  
%  OUTPUT
%
%  chi_bild      chimera bild file as sel variable
%  ptx           array of x positions from pdb_tail  
%  pty           array of y positions from pdb_tail  
%  ptz           array of z positions from pdb_tail  
%  vx            array of x vectors (pdb_head - pdb_tail) 
%  vy            array of y vectors (pdb_head - pdb_tail)  
%  vz            array of z vectors (pdb_head - pdb_tail)  
%  len           array of length
%  pdb_head      aligned head pdb
%
%EXAMPLE
% proteasome = pdbread('http://www.rcsb.org/pdb/files/1ryp.pdb'); 
% [chibi ptx pty ptz vx vy vz]=tom_pdbs2VectField(proteasome,proteasome);
% quiver3(ptx,pty,ptz,vx,vy,vz,0.5);
%
%
%NOTE
%
% positions must be in st.Model.Atom(1..n).X (Y,Z)
%
%REFERENCES
%
%SEE ALSO
%   quiver3,tom_vol2chimera
%
%   created by FB 01/24/06
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


if (nargin < 3)
    pre_align = 0;
end;

if (nargin < 4)
    pixelsize = 1;
end;

if (nargin < 5)
    dim = 0;
end;

if (nargin < 6)
    verbose = 1;
end;

if (isstruct(pdb_tail)==0)
    if (verbose==1)
        disp(['reading: ' pdb_tail]);
    end
    pdb_tail=pdbread(pdb_tail);
end;

if (isstruct(pdb_head)==0)
    if (verbose==1)
        disp(['reading: ' pdb_head]);
    end
    pdb_head=pdbread(pdb_head);
end;


%check size
if (length(pdb_tail.Model.Atom) ~= length(pdb_head.Model.Atom))
    error(['Number of atoms in head and tail differ !']);
end;


if (dim~=0)
    pdb_tail=pre_center_mod(pdb_tail);
    pdb_head=pre_center_mod(pdb_head);
end;

num_of_atoms=length(pdb_tail.Model.Atom);

if (verbose==1)
    disp([num2str(num_of_atoms) ' atoms found:']);
end;

mid=floor(dim/2);

phx=([pdb_head.Model.Atom(:).X]./pixelsize)+mid;
phy=([pdb_head.Model.Atom(:).Y]./pixelsize)+mid;
phz=([pdb_head.Model.Atom(:).Z]./pixelsize)+mid;

ptx=([pdb_tail.Model.Atom(:).X]./pixelsize)+mid;
pty=([pdb_tail.Model.Atom(:).Y]./pixelsize)+mid;
ptz=([pdb_tail.Model.Atom(:).Z]./pixelsize)+mid;

if (pre_align==1)
    Options.Registration='Rigid';
    Options.TolX=0.0000000001;
    Options.TolP=0.0001;
    Options.Optimizer='fminsearch';
    [ph_mov,M]=tom_cp2tform(cat(1,ptx,pty,ptz)', cat(1,phx,phy,phz)', Options);
    phx=ph_mov(:,1)';
    phy=ph_mov(:,2)';
    phz=ph_mov(:,3)';
    pdb_head_org=pdb_head;
    for i=1:length(phx)
        pdb_head.Model.Atom(i).X=phx(i);
        pdb_head.Model.Atom(i).Y=phy(i);
        pdb_head.Model.Atom(i).Z=phz(i);
    end;
end;



vx=phx-ptx;
vy=phy-pty;
vz=phz-ptz;

len=sqrt((vx.^2)+(vy.^2)+(vz.^2));

chi_bild=cell(num_of_atoms+1,1);
chi_bild{1}='.color 1 0 0';


if (verbose==1)
    disp(['mean vect length: ' num2str(mean(len))]);
    
end;

for i=1:num_of_atoms
    if ((vx(i)==0 && vy(i)==0 && vz(i)==0)==0 )
        chi_bild{i+1}=['.arrow ' num2str(ptx(i)) ' '  num2str(pty(i)) ' ' num2str(ptz(i)) ' ' num2str(phx(i)) ' ' num2str(phy(i)) ' ' num2str(phz(i)) ];
    end;
    if (mod(i,10000)==0)
        disp(num2str(i));
    end;
end;


disp(' ');



function pdbdata=pre_center_mod(pdbdata)
minx=min([pdbdata.Model.Atom.X]);
miny=min([pdbdata.Model.Atom.Y]);
minz=min([pdbdata.Model.Atom.Z]);
maxx=max([pdbdata.Model.Atom.X]);
maxy=max([pdbdata.Model.Atom.Y]);
maxz=max([pdbdata.Model.Atom.Z]);
for i=1:size(pdbdata.Model.Atom,2)
    pdbdata.Model.Atom(i).X = round(pdbdata.Model.Atom(i).X - ((maxx + minx) / 2));
    pdbdata.Model.Atom(i).Y = round(pdbdata.Model.Atom(i).Y - ((maxy + miny) / 2));
    pdbdata.Model.Atom(i).Z = round(pdbdata.Model.Atom(i).Z - ((maxz + minz) / 2));
end;



