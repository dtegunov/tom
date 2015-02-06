function openvol(name)
%OPENVOL creates ...
%
%   openvol(name)
%
%PARAMETERS
%
%  INPUT
%   name                ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_amira_createisosurface(...);
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

in=(tom_emreadc(name));

last=findstr(name,'.vol')-1;
if isequal(computer,'PCWIN')
    first=findstr(name,'\');
else
    first=findstr(name,'/');
end;
first=first(end)+1;
name=strrep(name,'.','_');
assignin('base',name(first:last),in);



