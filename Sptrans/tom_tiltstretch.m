function o=tom_tiltstretch(i,Tiltaxis,phi_current,phi_previous);
%TOM_TILTSTRETCH applies an projective stretch to an image
%
%   o=tom_tiltstretch(i,Tiltaxis,phi_previous,phi_current);
%
%PARAMETERS
%
%  INPUT
%   i                   input projection image
%   Tiltaxis            tiltaxis of the tomographic tiltseries
%   phi_previous        previous tiltangle
%   phi_current         actual projection image
%  
%  OUTPUT
%   o                   tilt corrected image
%
%EXAMPLE
%   in=tom_emread('proj_1.em');
%   o=tom_tiltstretch(in.Value,0,0,60); % stretch the image recorded at 0
%   % degree according to 60 degree tilt. The tiltaxis is 0 degree therefore
%   % parallel to the x-axis.
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN 09/03/05
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


AngleFactor=1;
Tiltaxis=Tiltaxis.*pi./180;
phi_current=phi_current.*pi./180;
phi_previous=phi_previous.*pi./180;

SF=cos(phi_previous.*AngleFactor)./cos(phi_current.*AngleFactor);

B0=1+((SF-1).*sin(Tiltaxis).*sin(Tiltaxis));
B1=1+((SF-1).*cos(Tiltaxis).*cos(Tiltaxis));
B2=(1-SF).*cos(Tiltaxis).*sin(Tiltaxis);
dimx=size(i,1);
dimy=size(i,2);
o=zeros(size(i));
o=o+mean2(i);
B3=dimx./2.*(B0+B2-1);
B4=dimy./2.*(B2+B1-1);


for laufx=1:.5:dimx+.499
    for laufy=1:.5:dimy+.499
        newx=round((laufx.*B0)+(B2.*laufy-B3));
        newy=round((laufx.*B2)+(B1.*laufy-B4));
        if(newx>0 && newx<dimx && newy>0 && newy<dimy)
            o(newx,newy)=i(round(laufx),round(laufy));
        end;
    end;
end;
