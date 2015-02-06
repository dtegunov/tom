function tom_em2tifall
%TOM_EM2TIFALL convert all EM format file to TIF format
%
%   tom_em2tifall
%
%PARAMETERS
%
%  INPUT
%  
%  OUTPUT
%
%EXAMPLE
%   tom_em2tifall
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN 12/19/05
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
[pathname] = uigetdir(pwd,'Select the source-directory of your tiltseries');
if isequal(pathname,0) disp('Cancelled.'); return; end;
[filename pathname_out] = uiputfile({'*.*'}, 'Select path and name of the sorted tiltseries (.em is added)');
if isequal(filename,0) | isequal(pathname_out,0) disp('Cancelled.'); return; end;
newextension='.em';
newname=[pathname_out filename];
pathname=[pathname '\'];
end;

d=dir(pathname);
angles=zeros((size(d,1)-2),1);
laufx=0;
for lauf=3:(size(d,1))
    if (tom_isemfile([pathname d(lauf).name])==1)
        laufx=laufx+1;
        em=tom_reademheader([pathname d(lauf).name]);
        angles(laufx)=em.Header.Tiltangle;
    else    
        laufx=laufx+1;
        angles(laufx)=-9999; %just a dummy Value to keep the index in the right order :-)
    end;
   
end;
[y,Index]=sort(angles);
lauf2=1;
for lauf=1:(size(angles,1))
    if (tom_isemfile([pathname d(Index(lauf)+2).name]) == 1)
        newemname=[newname num2str(lauf2) newextension];
        %    dos(['copy ' path d(Index(lauf)+2).name ' ' newemname]);
        i=tom_emread([pathname d(Index(lauf)+2).name]);
        disp([pathname d(Index(lauf)+2).name ' >> ' newemname ' ']);
        o.Value=tom_xraycorrect(i.Value); 
        %    o.Value=i.Value; 
        disp('correcting x-rays with a std_dev of 10 (SN)');
        o.Header=i.Header;
        tom_emwrite(newemname,o);
        lauf2=lauf2+1;    
    end;
    
end;

