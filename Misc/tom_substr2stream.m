function stream=tom_substr2stream(main_st,sub_st,index)
%TOM_SUBSTR2STREAM creates ...
%
%   stream=tom_substr2stream(main_st,sub_st,index)
%
%PARAMETERS
%
%  INPUT
%   main_st             ...
%   sub_st              ...
%   index               ...
%  
%  OUTPUT
%   stream              ...
%
%EXAMPLE
%   ... = tom_substr2stream(...);
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

if (nargin>2)
    for i=1:size(index,2)
        stream(i)=getfield(main_st,{index(i)},sub_st);
    end;
else
    for i=1:size(main_st,2)
        stream(i)=getfield(main_st,{i},sub_st);
    end;
end;