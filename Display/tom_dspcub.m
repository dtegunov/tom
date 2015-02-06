function h = tom_dspcub(a,md,columns,range)
%TOM_DSPCUB shows the object step by step in z-dimension
%
%   h = tom_dspcub(a,md,columns,range)
%
%PARAMETERS
%
%  INPUT
%   a                   3D image
%   md                  mode of presentation, can be 0 (=default), 1 or 2
%                       0: xy-slices
%                       1: yz-slices
%                       2: xz-slices
%                       Use as in EM
%   columns             number of columns to display (optional)
%   range               contrast range of image (optional)
%  
%  OUTPUT
%   h           		handle to the image
%
%EXAMPLE
%   a=tom_emread('vol1.em');
%   tom_dspcub(a.Value);
%
%REFERENCES
%
%SEE ALSO
%   TOM_VOLXY, TOM_VOLXYZ, TOM_EMBROWSE
%
%   created by AF 25/09/02
%   updated by WDN & AL 28/10/05
%   updated by AK 13/01/06 added columns, range and output handle
%   updated by FBR 19/06/08
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

error(nargchk(0, 4, nargin, 'struct'))

warning off;

var=2;
if ~exist('md')
    md=0;
end

%FBR 080619
if isstruct(a)
    a = a.Value;
elseif ischar(a)
    a = tom_emreadc(a); a = a.Value;
end

if md==2
    a=double(permute(a,[3 2 1]));
elseif md==1
    a=double(permute(a,[3 1 2]));
else
    a=double(permute(a,[2 1 3]));
end

if nargin < 4
    a=tom_imadj(a);
end

[s1,s2,s3]=size(a);
a=shiftdim(a,var);
a=shiftdim(a,-1); a=shiftdim(a,2);

if nargin < 3
    montage(a);
else
    if nargin < 4
        hh = tom_montage(a,columns);
    else
        hh = tom_montage(a,columns,range);
    end
end

for kk=0:ceil(sqrt(s3))
    set(gca,'Ytick',[0:s1:kk*s1]);
    set(gca,'Xtick',[0:s2:kk*s2]);
end
set(gca,'XAxisLocation','top');
set(gca,'GridLineStyle','-');
set(gca,'XColor', [0.416 0.706 0.780]);
set(gca,'YColor', [0.416 0.706 0.780]);

axis on; grid on;

if nargout > 0 && exist('hh','var')==0
    hh=gca;
    h = hh;
end


if exist('hh','var')==1
    h =hh;
end;
    