function max_num=tom_get_max_numeric(input_cell,direction)
%TOM_AV2_GET_MAX_NUMERIC finds the max number in a cell array of strings
%
%   max_num=tom_get_max_numeric(input_cell)
%
%PARAMETERS
%
%  INPUT
%   input_cell          cell array with strings
%   direction           forward or backward
%
%
%  OUTPUT
%   max_num             maximum numeric
%
%EXAMPLE
% 
% tom_get_max_numeric({'55_18.em','t_2.em','t_3.em'})
% tom_get_max_numeric({'55_18.em','t_2.em','t_3.em'},'forward')
%   
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV2_STACKBROWSER
%
%   created by FB(eckster) 06/16/08
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

if nargin < 2
    direction='backward';
end;

zz=1;
for i=1:size(input_cell,2)
    tmp=input_cell{i};
   
    if strcmp(direction,'backward')
        interv=size(tmp,2):-1:1;
   else
        interv=1:size(tmp,2);
    end;
    
    
    for ii=interv
        tmpp=findstr(tmp(ii),'0 1 2 3 4 5 6 7 8 9');
        if (isempty(tmpp)==0)
            num(zz)=tmp(ii);
            zz=zz+1;
        end;
        if (zz > 1) && isempty(tmpp)
            break;
        end;
    end;
    
    if  exist('num','var')==0
        max_num=[];
        return;
    end;
    
    if strcmp(direction,'backward')
        for iii=1:size(num,2); 
            new(size(num,2)-iii+1)=num(iii);
        end;
        num=new;
        clear('new');
    end;
            
    all_idx(i)=str2double(num);
    clear('num');
    zz=1;
end;

max_num=max(all_idx);