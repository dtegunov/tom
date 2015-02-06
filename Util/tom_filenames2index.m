function index=tom_filenames2index(names)
%TOM_FILENAMES2INDEX extracts indices of filenames
%
%   index=tom_filenames2index(names)
%
%PARAMETERS
%
%  INPUT
%   names              cell,folder or xmipp sel file
%   
%  OUTPUT
%   index              inex array 
%
%
%EXAMPLE
%
%
%REFERENCES
%
%SEE ALSO
%   
%
%   created by SN 01/08/07
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


flag='xmipp_sel';

if iscell(names)
    flag='cell';
else
    flag='dir';
end;

if (strcmp(flag,'dir'))
    dd=dir(names);
    error('no implemeted!');
end;


if (strcmp(flag,'xmipp_sel'))
    error('no implemeted!');
end;

index=zeros(length(names),1);

parfor i=1:length(names)
    tmp_str=names{i};
    [a b c]=fileparts(tmp_str);
    pos=strfind(b,'_');
    index(i)=str2double(b(pos+1:end));
end;