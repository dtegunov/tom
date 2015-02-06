function tom_bindir(pathname, filename,pathname_out)
%TOM_BINDIR bins images of a directory
%
%   tom_bindir(pathname, filename,pathname_out)
%
%PARAMETERS
%
%  INPUT
%   pathname            ...
%   filename            ...
%   pathname_out        ...
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
%   created by SN 02/01/06
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


if nargin <1 
[pathname] = uigetdir(pwd,'Select the source-directory of your images');
if isequal(pathname,0) disp('Cancelled.'); return; end;
%[pathname_out] = uigetdir(pwd,'Select path of the binned images');
[filename pathname_out] = uiputfile(pwd,'Select path of the binned images');
if isequal(pathname_out,0) disp('Cancelled.'); return; end;
pathname=[pathname '/'];
end;

d=dir(pathname);
l=0;
for lauf=3:(size(d,1))
    if (tom_isemfile([pathname d(lauf).name])==1)
        l=l+1;
        in=tom_emreadc([pathname d(lauf).name]);
        in.Value=tom_filter(in.Value,3,'quadr','real');
        [mean, max, min, std, variance] = tom_dev(in.Value,'noinfo');
        in.Value=tom_limit(in.Value,mean-3.*std,mean+3.*std)-mean;
       in.Value=-tom_bin(in.Value,1);
       % in.Value=-(in.Value);
        in.Header.Size=in.Header.Size./2;
        in.Header.Objectpixelsize=in.Header.Objectpixelsize.*2;        
        in.Header.Comment=in.Header.Comment';
        tom_emwrite([pathname_out '/' filename num2str(l) '.em'],in);
    end;   
end;
