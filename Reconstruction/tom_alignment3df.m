function tom_alignment3df(Markerfile,imdim)
%TOM_ALIGNMENT3DF transfers data into stream
%
%   tom_alignment3df(Markerfile,imdim)
%
%PARAMETERS
%
%  INPUT
%   Markerfile          ...
%   imdim               ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_alignment3df(...);
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


zz=1;
for iz=1:size(Markerfile,3)
    for iy=1:size(Markerfile,2)
         for ix=1:size(Markerfile,1)
           AR(zz)=Markerfile(ix,iy,iz);
           zz=zz+1;
        end;
    end;
end;

NA(1)=0;
NA(2)=size(Markerfile,1);
NA(3)=size(Markerfile,2);
NA(4)=size(Markerfile,3);


IR(2)=imdim(1);

tom_alignment3dinf(AR,NA,AR);
 
