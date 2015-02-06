function out=tom_isimagicfile(filename)
%TOM_isimagicfile creates ...
%
%   out = tom_isimagicfile (filename)
%
%PARAMETERS
%
%  INPUT
%   filename            name of the file
%  
%  OUTPUT
%   out                 0/1 
%
%EXAMPLE
%   
%   out=tom_isimagicfile(filename)
%
%REFERENCES
%
%SEE ALSO
%   tom_imagicread, tom_imagicwrite 
%
%   created by fb 28/01/06
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


[path dat_name ext]=fileparts(filename);

if (strcmp(ext,'.img')==0 & strcmp(ext,'.hed')==0 )
    out=0;
    return
end;

has_body=exist([path '/' dat_name '.img'],'file');
has_head=exist([path '/' dat_name '.hed'],'file');

if ((has_body + has_head) < 4)
    out=0;
    return;
end;

out=1;
