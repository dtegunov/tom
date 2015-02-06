function tom_av2_flip_stack(input_stack,output_stack)
%TOM_AV2_FLIP_STACK creates ...
%
%   tom_av2_flip_stack(input_stack,output_stack)
%
%PARAMETERS
%
%  INPUT
%   input_stack         ...
%   output_stack        ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av2_flip_stack(...);
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


stack_h=tom_reademheader(input_stack);
sz_stack=stack_h.Header.Size;
zz=1;
tom_emwritec(output_stack,[sz_stack(1) sz_stack(2) (sz_stack(3).*2)],'new');

for i=1:sz_stack(3)

    part=tom_emread(input_stack,'subregion',[1 1 i],[sz_stack(1)-1 sz_stack(2)-1 0]);
    part_f=tom_mirror(part.Value,'x');
    tom_emwritec(output_stack,part.Value,'subregion',[1 1 zz],[sz_stack(1) sz_stack(2) 1]);
    zz=zz+1;
    tom_emwritec(output_stack,part_f,'subregion',[1 1 zz],[sz_stack(1) sz_stack(2) 1]);
    zz=zz+1;
    i
end;
