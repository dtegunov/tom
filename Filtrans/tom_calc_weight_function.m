function w_func=tom_calc_weight_function(dims,all_angles,thickness,angle_proj)
%TOM_CALC_WEIGHT_FUNCTION calculates 2D and yes yes 3D Weighting function
%
%   w_func=tom_calc_weight_function(dims,all_angles,thickness,angle_proj)
%
%PARAMETERS
%
%  INPUT
%   dims                Dimension of the weighting function ... works as a switch between 2d and 3d
%   all_angles          All Angles used in the backproj phi and the or euler angles or 
%                        'analytical' for analytical weighting 
%   thickness           Sample Thickness
%   angle_proj          parameter for 2d angle must be an element of ALL Angles
%  
%  OUTPUT
%   w_func              weighting function 2d or 3d
%
%EXAMPLE
%   2d:
%   w_func=tom_calc_weight_function([64 64],[10 20; 30 40; 50 10],64,[10 20]);  
%   3d:
%   w_func=tom_calc_weight_function([64 64 64],[10 20; 30 40; 50 10],64);
%
%
%NOTE:
%   All input types are double. Returned w_func is double
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
%   created by SN, FB after a Fortran routine by RH 11/27/05
%   updated by SN and FB 14/01/10
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

%
if (dims(1)~=dims(2))
    error('only equal dimensions in x and y supported.');
end;



if (isnumeric(all_angles)==1)
    %check out exact weighting
    
    all_angles=double(all_angles);
    thickness=double(thickness);
    
    if (nargin==3)
        angle_proj=[0 0];
    end;
     angle_proj=double(angle_proj);
    
    %check inputs
    if  ((size(dims,2)==2) && (nargin==3))
        error('projection angle needed for 2d weighting');
    end;
    
    if  ((size(dims,2)==3) && (nargin==4))
        error('no projection angle needed for 3d weighting');
    end;
    
    %transform input angles
    if (size(all_angles,2)==3) % here are eulers
        tmp=all_angles;
        clear all_angles;
        all_angles(:,1)=90-tmp(:,1); % phi_new=90-euler_phi; Kipprichtung !!!!
        all_angles(:,2)=tmp(:,3); % kw=euler_theta Kippwinkel !!!
    end;
    
    if (size(angle_proj,2)==3) % here are eulers
        tmp=angle_proj;
        clear angle_proj;
        angle_proj(1)=90-tmp(1); % phi_new=90-euler_phi; Kipprichtung !!!!
        angle_proj(2)=tmp(3); % kw=euler_theta Kippwinkel !!!
    end;
    
    %check if angle_proj is included in all_angles
    if (size(dims,2)==2)
        is_in=0;
        for i=1:size(all_angles,1)
            chk_sum=sum(all_angles(i,:)==angle_proj);
            if (chk_sum==2)
                is_in=1;
                break;
            else
                chk_sum=0;
            end;
        end;
        if (is_in==0)
            error('proj_angle must be element of all_angles');
        end;
    end;
    
    if (size(dims,2)==2)
        w_func=weight2d(dims,all_angles,thickness,angle_proj);
    else
        w_func=weight3d(dims,all_angles,thickness);
    end;
    
    w_func=fftshift(w_func);
    
else
    %check out  analytical weighting
    [x,y]=ndgrid(-dims(1)/2:dims(1)/2-1,-dims(2)/2:dims(2)/2-1);
    if ( size(dims,2)==3)
        w_func=zeros(dims);
        for ii=1:dims(3)
            w_func(:, :,ii)=x;
        end;
    else
        w_func=x;
    end;
    w_func=tom_norm(abs(w_func),1);
end;


function [w_func]=weight3d(dims,angles,thickness)


dimx=dims(1);
dimy=dims(2);
dimz=dims(3);
NDIM=dimx;
IT=1;
THICK=thickness;
if (thickness<=0) THICK=dimx; end;
WIDTH=pi; % cut off for sinc function, pi = first zero of sinc, extend to ...
%WIDTH=0.50.*pi;
BR=zeros(dimx.*dimy.*dimz,1);
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
        for IX=1:dimx
            SUM=0;
            ITLT=IT+1;
            for I=1:dimangles
                ARG=XST.*TILT(ITLT,1)+YST.*TILT(ITLT,2)+ZST.*TILT(ITLT,3);
                if abs(ARG)<1e-6
                    ADD=1;
                elseif abs(ARG)>1e-6 && abs(ARG)<WIDTH
                    ADD=sin(ARG)./ARG;
                else
                    ADD=0;
                end;
                SUM=SUM+ADD;
                ITLT=ITLT+1;
            end;
            BR(IND)=1./max(SUM,1.0);
            IND=IND+1;
            XST = XST + 1./(dimx);
            if IX==floor(dimx./2)
                XST = -XST;
            end;
        end; 
        YST = YST + 1./(dimy);
        if IY==floor(dimy./2)
            YST = -YST;
        end;
    end; 
    ZST=ZST+1./dimz;
    if IZ==floor(dimz./2)
        ZST=-ZST;
    end;
end;

w_func=reshape(BR,[dimx dimy dimz]);

function [w_func]=weight2d(dims,angles,thickness,angles_projection)


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
BR=ones(dimx.*dimy,1); % !!!!!!!!!!!!!!

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
    for IX=1:dimx
        SUM=0;
        ITLT=1;
        XST = UST*A1 - VST*SPHI;
        YST = UST*A2 + VST*CPHI;
        ZST = -UST*STHE;
        for I=1:dimangles
            ARG=XST.*TILT(ITLT,1)+YST.*TILT(ITLT,2)+ZST.*TILT(ITLT,3);
            if abs(ARG)<1e-6
                ADD=1;
            elseif abs(ARG)>1e-6 && abs(ARG)<pi
                ADD=sin(ARG)./ARG;
            else
                ADD=0;
            end;
            SUM=SUM+ADD;
            ITLT=ITLT+1;            
        end;
        BR(IND)=BR(IND)./max(SUM,1.0);
        IND=IND+1;
        UST=UST+1./NDIM;
        if IX==floor(dimx./2)
            UST=-UST;
        end;
    end;
    VST=VST+1./NDIM;
    if IY==floor(dimy./2)
        VST=-VST;
    end;
end;

w_func=reshape(BR,[dimx dimy]);




