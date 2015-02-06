function out=tom_adapthist3D(in,numtiles,distribution)
%TOM_ADAPTHIST3D creates ...
%
%   out=tom_adapthist3D(in,numtiles,distribution)
%
%PARAMETERS
%
%  INPUT
%   in                  ...
%   numtiles            ...
%   distribution        ...
%  
%  OUTPUT
%   out         		...
%
%EXAMPLE
%   ... = tom_adapthist3D(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN 10/23/04
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

for lauf=1:size(in,3)
    out(:,:,lauf)=adapthisteq(in(:,:,lauf),'NumTiles',numtiles,'Distribution',distribution);
end;


    
    