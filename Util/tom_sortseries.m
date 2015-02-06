function tom_sortseries(pathname, newname,newextension,x_ray)
%TOM_SORTSERIES sorts EM-files by tiltangle
%
%   tom_sortseries(pathname, newname,newextension)
%
%    Copies and sorts series of EM-Image Files (V-Format)
%    useful to sort tiltseries in ascending tilt-angle order.
%    TOM_SORTSERIES reads an entire directory!
%    All images are also checked for x-rays by TOM_XRAYCORRECT.
%
%PARAMETERS
%
%  INPUT
%   pathname            ...
%   newname             ...
%   newextension        ...
%   x_ray               (opt.) (1)  switch x_ray corr on/off
%
%  OUTPUT
%
%EXAMPLE
%   tom_sortseries('./tmp/','./new_tmp/new_','.em');
%   tom_sortseries         
%    opens two fileselect-boxes for source and destination
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN 08/07/02
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

if (nargin < 4)
    x_ray=1;
end;


if nargin <1
    [pathname] = uigetdir(pwd,'Select the source-directory of your tiltseries');
    if isequal(pathname,0) disp('Cancelled.'); return; end;
    [filename pathname_out] = uiputfile({'*.*'}, 'Select path and name of the sorted tiltseries (.em is added)');
    if isequal(filename,0) | isequal(pathname_out,0) disp('Cancelled.'); return; end;
    newextension='.em';
    newname=[pathname_out filename];
    if (isunix==1)
         pathname=[pathname '/'];
    else
        pathname=[pathname '\'];
    end;
end;

d=dir(pathname);
angles=zeros((size(d,1)-2),1);
laufx=0;
for lauf=3:(size(d,1))
    if (tom_isemfile([pathname '/' d(lauf).name])==1)
        laufx=laufx+1;
        em=tom_reademheader([pathname '/' d(lauf).name]);
        angles(laufx)=em.Header.Tiltangle;
    else    
        laufx=laufx+1;
        angles(laufx)=-9999; %just a dummy Value to keep the index in the right order :-)
    end;
   
end;
[y,Index]=sort(angles);
lauf2=1;
for lauf=1:(size(angles,1))
    if (tom_isemfile([pathname '/' d(Index(lauf)+2).name]) == 1)
        newemname=[newname num2str(lauf2) newextension];
        %    dos(['copy ' path d(Index(lauf)+2).name ' ' newemname]);
        i=tom_emreadc2([pathname '/' d(Index(lauf)+2).name]);
        disp([pathname  '/' d(Index(lauf)+2).name ' >> ' newemname ' ']);
        if (x_ray==1)
            o.Value=tom_xraycorrect(i.Value); 
             %    o.Value=i.Value; 
            disp('correcting values with a std_dev of 10 around the mean value (SN)');
        else
           o.Value=i.Value; 
           disp(['No x-ray correction done ' ]);
        end;
        o.Header=i.Header;
        tom_emwrite(newemname,o);
        lauf2=lauf2+1;    
    end;
    
end;

