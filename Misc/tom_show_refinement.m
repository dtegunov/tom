function pic=tom_show_refinement
%TOM_SHOW_REFINEMENT displays small 3D volumes in a 2D row image
%
%   pic=tom_show_refinement
%
%PARAMETERS
%
%  INPUT
%  
%  OUTPUT
%   pic                 ...
%
%EXAMPLE
%   ... = tom_show_refinement
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

directoryname = uigetdir(pwd,'select directory');
d=dir(directoryname);

idx=1;
max_ref=size(d,1)+1;
min_ref=1;
for i=3:size(d,1)
    in=tom_emread(char(d(i).name));
    ref(:,:,:,idx)=tom_norm(tom_filter(in.Value,3),1);
    idx=idx+1;
end;

pic=ones(size(ref,1).*size(ref,3)+size(ref,3).*2+2,size(ref,2).*(max_ref-min_ref)+size(ref,4).*2+2).*min(min(min(min(ref))));


iy=3;
sx=size(ref,1);
for i=1:size(ref,4)
    ix=3;
    for iz=1:size(ref,3)
        pic(ix:ix+sx-1,iy:iy+sx-1)=squeeze(ref(:,:,iz,i));
        ix=ix+size(ref,1)+2;
    end;
    iy=iy+size(ref,2)+2;
end;