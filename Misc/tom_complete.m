function out=tom_complete(in)
%TOM_COMPLETE creates ...
%
%   out=tom_complete(in)
%
%PARAMETERS
%
%  INPUT
%   in                  ...
%  
%  OUTPUT
%   out                 ...
%
%EXAMPLE
%   ... = tom_complete(...);
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




if (size(in,3) >1)
    new=zeros(size(in,1)*2,size(in,2),size(in,3));
    newp=tom_paste(new,in,[1 1 1]);
    clear new;
    %out=fftshift((tom_mirror(tom_mirror(newp,'x'),'z')+newp));
    %out=tom_shift(fftshift((tom_mirror(tom_mirror(newp,'x'),'z'))),[0 0 1])+fftshift(newp);
    %out=tom_move(fftshift(tom_mirror((tom_mirror(tom_mirror(newp,'x'),'z')),'y')),[0 -2 0] ) + fftshift(newp);
    out=fftshift(newp)+tom_move(tom_rotate(fftshift(newp),[-1 0 0;0 -1 0; 0 0 -1],'linear'),[-1 0 0]) ;
else
    dimx=size(in,1);
    dimy=size(in,2);
    new=zeros(2.*dimx,dimy);
    newp=tom_paste(new,tom_mirror(in,'y'),[1 1]);
    clear new;
    %out=fftshift((tom_mirror(tom_mirror(newp,'x'),'y')+newp));
    % out=tom_shift(fftshift((tom_mirror(tom_mirror(newp,'x'),'y'))),[0 1]) + fftshift(newp);
    % out=tom_move(fftshift((tom_mirror(tom_mirror(newp,'x'),'y'))),[0 0]) + fftshift(newp);
    out=fftshift(newp) + imrotate(fftshift(newp),180); 
end;



 
 