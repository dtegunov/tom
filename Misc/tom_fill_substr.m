function main_st=tom_fill_substr(main_st,sub_st,in,index)
%TOM_FILL_SUBSTR creates ...
%
%   main_st=tom_fill_substr(main_st,sub_st,in,index)
%
%PARAMETERS
%
%  INPUT
%   main_st             ...
%   thressub_sthold     ...
%   in                  ...
%   index               ...
%  
%  OUTPUT
%   main_st             ...
%
%EXAMPLE
%   main_st=tom_fill_substr(main_st,sub_st,in,index)
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
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




if (nargin>3)
    for i=1:size(index,2)
        if (size(in,2)==1)
            main_st=setfield(main_st,{index(i)},sub_st,in);
         else
            main_st=setfield(main_st,{index(i)},sub_st,in(i));
        end;
    end;
else
    for i=1:size(main_st,2)
        if (size(in,2)==1)
            main_st=setfield(main_st,{i},sub_st,in);
        else
            main_st=setfield(main_st,{i},sub_st,in(i));
        end;
        
        end;
end;