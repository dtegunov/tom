function [weight3d]=tom_weight3d_exact(dims,angles,thickness);
%TOM_WEIGHT3D_EXACT creates an 'exact' weighting function in 3d
%
%   [weight3d]=tom_weight3d_exact(dims,angles,thickness)
%
%PARAMETERS
%
%  INPUT
%   dims                dimensions of the weighting volume as a vector
%   angles              full angular set of projections (phi,theta) or Euler angles [Phi Psi Theta]
%   thickness           object thickness
%  
%  OUTPUT
%   weight3d            weighting function
%
%   apply to an unweighted tomographic reconstruction by:
%
%   weighted_volume=real(ifftn(fftn(unweighted_volume).*fftshift(weight3d)));
%
%EXAMPLE
%   [weight3d]=tom_weight3d_exact([64 64 64],[10 20 30; 40 50 60],64);
%
%   create a missing wedge:
%   tiltaxis:
%               angles(1,:)=ones(1,27).*90;
%   tiltangles:
%               angles(2,:)==[-65:5:65];
%   weighting function:            
%               [weight3d]=tom_weight3d_exact([32 32 32],angles',32);
%
%REFERENCES
%
%SEE ALSO
%   Radermacher et. al., 1986, A new 3-D reconstruction scheme
%   applied to to the 50S ribosomal subunit of E.coli, J.Microssc.
%   141:RP1-RP2
%   and Electron Tomography: Three-Dimensional Imaging with the TEM, edited
%   by J. Frank, M.Radermacher, p91-115
%
%   created by SN and FB after a Fortran routine by RH 04/11/05
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

if size(angles,2)==3 % here are eulers

    tmp=angles;
    clear angles;
    angles(:,1)=90-tmp(:,1); % phi_new=90-euler_phi; Kipprichtung !!!!
    angles(:,2)=tmp(:,3); % kw=euler_theta Kippwinkel !!!

end;
dimx=dims(1);
dimy=dims(2);
dimz=dims(3);
NDIM=dimx;
IT=1;
THICK=thickness;
if (thickness<=0) THICK=dimx; end;
WIDTH=pi; % cut off for sinc function, pi = first zero of sinc, extend to ...
BR=zeros(dimx./2.*dimy.*dimz,1);
ITLT = 0;
dimangles=size(angles,1);
for I=1:dimangles
    PHI=angles(ITLT+IT,1).*pi./180;
    THETA=angles(ITLT+IT,2).*pi./180;
    ZW=sin(THETA).*pi.*THICK;
    TILT(I,1)=ZW.*cos(PHI);
    TILT(I,2)=ZW.*sin(PHI);
    TILT(I,3)=cos(THETA).*pi.*THICK;
    ITLT=ITLT+1;
end;

IND = 1;
ZST = 0;
IT=0;
for IZ=1:dimz
    YST=0;
    for IY=1:dimy
        XST=0;
        for IX=1:dimx./2
            SUM=0;
            ITLT=IT+1;
            for I=1:dimangles
                ARG=XST.*TILT(ITLT,1)+YST.*TILT(ITLT,2)+ZST.*TILT(ITLT,3);
                if abs(ARG)<1e-6
                    ADD=1;
                elseif abs(ARG)>1e-6 & abs(ARG)<WIDTH
                    ADD=sin(ARG)./ARG;
                else
                    ADD=0;
                end;
                SUM=SUM+ADD;
                ITLT=ITLT+1;
            end; % for I=1:dimangles
            if (SUM>0)
                BR(IND)=1./SUM;
                if SUM<1
                    BR(IND)=1;
                end;
            else
                BR(IND)=0;
            end;
            IND=IND+1;
            XST = XST + 1./(NDIM);
        end; % for IX=1:dimx./2
        YST = YST + 1./(NDIM);
        if IY==dimy./2
            YST = -YST;
        end;
    end; %for IY=1:dimy
    ZST=ZST+1./dimz;
    if IZ==dimz./2
        ZST=-ZST;
    end;
end;

weight3d=reshape(BR,[dimx./2 dimy dimz]);
weight3d=swap_weight(weight3d);

% Anwendung auf Volumen:
% weighted_volume=real(ifftn(fftn(unweighted_volume).*fftshift(weight3d)));


function out=swap_weight(vol_weighted)

% do the fftshift and a complete fourier space starting from the complex
% reduced one.

new=zeros(size(vol_weighted,1)*2,size(vol_weighted,2),size(vol_weighted,3));
newp=tom_paste(new,vol_weighted,[1 1 1]);
out=fftshift((tom_mirror(tom_mirror(newp,'x'),'z')+newp));

