function pic=tom_show_refinemt(min_ref, max_ref, name)
%TOM_SHOW_REFINEMT displays small 3D volumes in a 2D row image
%
%   pic=tom_show_refinemt(min_ref, max_ref, name)
%
%PARAMETERS
%
%  INPUT
%   min_ref             ...
%   max_ref             ...
%   name                ...
%  
%  OUTPUT
%   pic                 ...
%
%EXAMPLE
%   ... = tom_show_refinemt(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN 06/20/05
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


idx=1;

for i=min_ref:max_ref
    in=tom_emread([name '_' num2str(i) '.em']);
    ref(:,:,:,idx)=tom_norm(in.Value,1);
    idx=idx+1;
end;

pic=ones(size(ref,1).*size(ref,3)+size(ref,3).*2+2,size(ref,2).*max_ref-min_ref+size(ref,4).*2+2).*min(min(min(min(ref))));


iy=3;
for i=1:size(ref,4)
    ix=3;
    for iz=1:size(ref,3)

        pic=tom_paste(pic,ref(:,:,iz,i),[ix iy]);
        ix=ix+size(ref,1)+2;
    end;
    iy=iy+size(ref,2)+2;
end;