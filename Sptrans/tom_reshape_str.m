function vol=tom_reshape_str(stream,size_vol,order)
%TOM_RESHAPE_STR allocates some memory
%
%   vol=tom_reshape_str(stream,size_vol,order)
%
%PARAMETERS
%
%  INPUT
%   stream              ...
%   size_vol            ...
%   order               ...
%  
%  OUTPUT
%   vol                 ...
%
%EXAMPLE
%   ... = tom_reshape_str(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
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

vol=zeros(size_vol);


if nargin==2
    order='xyz';
end;

zz=1;
for z=1:size_vol(3)
    for y=1:size_vol(2)
        for x=1:size_vol(1)
         
            if strcmp(order,'xyz')
                vol(x,y,z)=stream(zz);
            end;
            
         zz=zz+1;
        end;
    end;
end;
