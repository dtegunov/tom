function [align2d stack_uncompr]=tom_av2_uncompress_stack(al_compr,al_org,stack_org,f_al_uncmpr,f_stack_uncompr)
%tom_av2_uncompress_stack uncompresses a stack   
%
%   [align2d stack_uncompr]=tom_av2_uncompress_stack(al_compr,al_org,stack_org,f_al_uncmpr,f_stack_uncompr)
%PARAMETERS
%
%  INPUT
%   al_compr            compressed align2d struct or filenaem
%   al_org              uncompressed align2d struct or filename
%   stack_org           org stack or filename   
%   f_al_uncmpr         filename for HD output  
%   f_stack_uncompr     filename for HD output       
% 
%
%  OUTPUT
%   stack_uncompr          compressed stack
%   al_uncmpr              align2d struct for stackbrowser  
%   
%
%
%EXAMPLE
%   
%   
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_compress_stack,tom_av2_stackbrowser
%
%   created by fb (eckster)
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

if (nargin<4)
    f_al_uncmpr='';
end;

if (nargin<5)
    f_stack_uncompr='';
end;



if (ischar(al_compr))
    tmp=load(al_compr);
    al_compr=tmp.align2d;
end;

if (ischar(al_org))
    tmp=load(al_org);
    al_org=tmp.align2d;
end;

if (ischar(stack_org) && isempty(stack_org)==0)
    tmp=tom_emreadc(stack_org);
    stack_org=tmp.Value;
end;

idx=al_compr(1,1).dataset;
for i=2:size(al_compr,2)
    idx=cat(1,idx,al_compr(1,i).dataset);
end;


align2d=al_org(1,idx);

if isempty(stack_org)==0
    stack_uncompr=stack_org(:,:,idx);
end

if (isempty(f_al_uncmpr)==0)
    save(f_al_uncmpr,'align2d');
end;

if (isempty(f_stack_uncompr)==0)
    tom_emwrite(f_stack_uncompr,stack_uncompr);
end;





