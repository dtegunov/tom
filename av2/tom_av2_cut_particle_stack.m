function new_stack=tom_av2_cut_particle_stack(stack,cut_size)
%tom_av2_cut_particle_stack doubles a paticle stack by cutting left and
%right side makes only sense for C2 symmetry 
%
%   new_stack=tom_av2_cut_particle_stack(stack,cut_size)
%
%PARAMETERS
%
%  INPUT
%   stack               input particle stack
%   radius              size of the new stack  
%  
%  OUTPUT
%   new_stack           cutted stack
%
%EXAMPLE
%
%im=tom_emread('murat_stack.em');
%new_im=tom_av2_cut_particle_stack(im,[96 96]);
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by fb 31/01/07
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


new_stack_left=zeros(cut_size(1),cut_size(2),size(stack,3));
new_stack_right=zeros(cut_size(1),cut_size(2),size(stack,3));

rest=size(stack,2)-cut_size(2);
cut2_start=round(rest./2);
cut2_stop=round(size(stack,2)-(rest./2))-1;

for i=1:size(stack,3)
    new_stack_left(:,:,i)=stack(1:cut_size(1),cut2_start:cut2_stop,i);
end;

new_start=size(stack,1)-cut_size(1)+1;

for i=1:size(stack,3)
    new_stack_right(:,:,i)=stack(new_start:size(stack,1),cut2_start:cut2_stop,i);
end;


clear('stack');

%build new stack out of left and right
new_stack=zeros(cut_size(1),cut_size(2),size(new_stack_right,3)+size(new_stack_left,3));
new_stack(:,:,1:size(new_stack_left,3))=new_stack_left;
new_stack(:,:,size(new_stack_left,3)+1:size(new_stack_left,3)+size(new_stack_right,3))=new_stack_right;