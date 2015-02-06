function [weight2d]=tom_weight2d_exact(dims,angles_projection,angles,thickness);
%TOM_WEIGHT2D_EXACT creates an 'exact' weighting function in 2d
%
%   [weight2d]=tom_weight2d_exact(dims,angles_projection,angles,thickness)
%
%PARAMETERS
%
%  INPUT
%   dims                dimensions of the weighting volume as a vector
%   angles_projection   angular set of projection (phi,theta) or Euler angles [Phi Psi Theta], must be subset of angles !!!!!!!!!!!!!
%   angles              full angular set of projections (phi,theta) or Euler angles [Phi Psi Theta]
%   thickness           object thickness
%  
%  OUTPUT
%   weight2d            weighting function
%
%   apply to an unweighted tomographic projection image by:
%
%   weighted_proj=real(ifft2(fft2(unweighted_proj).*fftshift(weight2d)));
%
%EXAMPLE
%           [weight2d]=tom_weight2d_exact([64 64],[10 20 30] ,[10 20 30; 40 50 60],64);
%
%REFERENCES
%
%SEE ALSO
%   ...
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

if size(angles_projection,2)==3 % here are eulers

    tmp=angles_projection;
    clear angles_projection;
    angles_projection(1)=90-tmp(1); % phi_new=90-euler_phi; Kipprichtung !!!!
    angles_projection(2)=tmp(3); % kw=euler_theta Kippwinkel !!!
    
end;

SPHI = sin(angles_projection(1).*pi./180);
CPHI = cos(angles_projection(1).*pi./180);
STHE = sin(angles_projection(2).*pi./180);
CTHE = cos(angles_projection(2).*pi./180);
A1 = CTHE*CPHI;
A2 = CTHE*SPHI;

IND = 1;
VST = 0;
ITLT = 0;
dimx=dims(1);
dimy=dims(2);
NDIM=dimx;
THICK=thickness;
IT=1;
dimangles=size(angles,1);
BR=ones(dimx./2.*dimy,1);

for I=1:dimangles
    PHI=angles(ITLT+IT,1).*pi./180;
    THETA=angles(ITLT+IT,2).*pi./180;
    ZW=sin(THETA).*pi.*THICK;
    TILT(I,1)=ZW.*cos(PHI);
    TILT(I,2)=ZW.*sin(PHI);
    TILT(I,3)=cos(THETA).*pi.*THICK;
    ITLT=ITLT+1;
end;

for IY=1:dimy
    UST=0;
    for IX=1:dimx./2
        SUM=0;
        ITLT=1;
        XST = UST*A1 - VST*SPHI;
        YST = UST*A2 + VST*CPHI;
        ZST = -UST*STHE;


        for I=1:dimangles
            ARG=XST.*TILT(ITLT,1)+YST.*TILT(ITLT,2)+ZST.*TILT(ITLT,3);
            if abs(ARG)<1e-6
                ADD=1;
            elseif abs(ARG)>1e-6 & abs(ARG)<pi
                ADD=sin(ARG)./ARG;
            else
                ADD=0;
            end;
            SUM=SUM+ADD;
            ITLT=ITLT+1;
            
        end;
%        if SUM==0
%                BR(IND)=0;
%        else
            BR(IND)=BR(IND)./SUM;
%        end;
        IND=IND+1;
        UST=UST+1./NDIM; % ????????????
    end;
    VST=VST+1./NDIM;
    if IY==dimx./2
        VST=-VST;
    end;
end;

weight2d=reshape(BR,[dimx./2 dimy]);
weight2d=swap_weight(weight2d);


function out=swap_weight(weight2d)
dimx=size(weight2d,1);
dimy=size(weight2d,2);

new=zeros(2.*dimx,dimy);

newp=tom_paste(new,tom_mirror(weight2d,'y'),[1 1]);
out=fftshift((tom_mirror(tom_mirror(newp,'x'),'y')+newp));





function out=swap_weight_old(weight2d)

dimx=size(weight2d,1);
dimy=size(weight2d,2);

out2=zeros(dimx,dimy);


ii=1;
for i=dimy./2+1:dimy
    out2(:,ii)=weight2d(:,i);
    ii=ii+1;
end;
ii=dimy./2+1;
for i=1:dimy/2
    out2(:,ii)=weight2d(:,i);
    ii=ii+1;
end;

out=zeros(dimx.*2,dimy);

ii=1;
for i=dimx+1:dimx*2
    out(i,:)=out2(ii,:);
    ii=ii+1;
end;

out2=out;
dimx2=size(out2,1);
dimy2=size(out2,2);

for ix=dimx2./2+1:dimx2
    for iy=dimy2./2+1:dimy2    
        out(dimy2-iy+1,dimx2-ix+1)=out2(ix,iy);
    end;
end;
for ix=dimx2./2+1:dimx2
    for iy=1:dimy2./2    
        out(dimy2-iy+1,dimx2-ix+1)=out(ix,iy);
    end;
end;

for ix=dimx2./2+1:dimx2
    for iy=dimy2./2+1:-1:1    
        out(iy,ix)=out2(ix,iy);
    end;
end;


    % replace by tom_mirror as in tom_weight3d_exact.m

