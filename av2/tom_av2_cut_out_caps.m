function tom_av2_cut_out_caps(stack,align,stack_out,align_o,cap_size,offset)
%TOM_AV2_CUT_OUT_CAPS creates ...
%
%   tom_av2_cut_out_caps(stack,align,stack_out,align_o,cap_size,offset)
%
%PARAMETERS
%
%  INPUT
%   stack               ...
%   align               ...
%   stack_out           ...
%   align_o             ...
%   cap_size            ...
%   offset              ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av2_cut_out_caps(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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


stack_h=tom_reademheader(stack);
sz_stack=stack_h.Header.Size;
mid(1)=round(sz_stack(1)./2)+1;
mid(2)=round(sz_stack(2)./2)+1;

load(align);



%allocate space for output
align_out=align2d;
align_out(1,size(align2d,2):2.*size(align2d,2)-1)=align2d;
tom_emwritec(stack_out,[cap_size(1) cap_size(2) (sz_stack(3).*2)],'new');
zz=1;

for i=1:sz_stack(3)
    image=tom_emread(stack,'subregion',[1 1 i],[sz_stack(1)-1 sz_stack(2)-1 0]);
    coord1(1)=mid(1)-offset(1)-round(cap_size(1)./2);
    coord1(2)=mid(2)-offset(2)-round(cap_size(2)./2);
    coord2(1)=mid(1)+offset(1)-round(cap_size(1)./2);
    coord2(2)=mid(2)+offset(2)-round(cap_size(2)./2);
    cap1=tom_cut_out(image.Value,coord1,cap_size,'no-fill');
    cap2=tom_cut_out(image.Value,coord2,cap_size,'no-fill');
    cap2=tom_mirror(cap2,'x');
    
    tom_emwritec(stack_out,cap1,'subregion',[1 1 zz],[cap_size(1) cap_size(2) 1]);
    align_out(1,zz)=align2d(1,i);
    align_out(1,zz).radius=max(cap_size);
    align_out(1,zz).position.x=align_out(1,i).position.x-offset(1)-round(cap_size(1)./2);
    align_out(1,zz).position.x=align_out(1,i).position.x-offset(2)-round(cap_size(2)./2);
    
    zz=zz+1;
    tom_emwritec(stack_out,double(cap2),'subregion',[1 1 zz],[cap_size(1) cap_size(2) 1]);
    align_out(1,zz)=align2d(1,i);
    align_out(1,zz).radius=max(cap_size);
    align_out(1,zz).position.x=align_out(1,i).position.x+offset(1)-round(cap_size(1)./2);
    align_out(1,zz).position.x=align_out(1,i).position.x+offset(2)-round(cap_size(2)./2);
    
    zz=zz+1;
    
    i
    
end;

save(align_o,'align_out');
