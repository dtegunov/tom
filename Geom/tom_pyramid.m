function py=tom_pyramid(sz)
%TOM_PYRAMID creates a volume with a pyramid
%
%   py=tom_pyramid(sz);
%
%PARAMETERS
%
%  INPUT
%   sz                 size of pyramid
%  
%  OUTPUT
%   py                 pyramid volume
%
%EXAMPLE
%   py = tom_pyramid(32);
%   tom_dspcub(py);
%
%REFERENCES
%
%SEE ALSO
%   tom_sphere
%
%   created by SN 03/27/09
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

py=zeros(sz,sz,sz,'single');

for i=1:sz
    r=ones(i,i);
    py(:,:,i)=tom_paste(py(:,:,i),r,[round(sz/2-i/2) round(sz/2-i/2) i]);
end;