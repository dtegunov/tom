function w_func=tom_calc_weight_functionc(dims,all_angles,thickness,angle_proj)
%TOM_CALC_WEIGHT_FUNCTIONC calculates 2D and yes yes 3D Weighting function
%
%   w_func=tom_calc_weight_functionc(dims,all_angles,thickness,angle_proj)
%
%PARAMETERS
%
%  INPUT
%   dims                Dimension of the weighting function ... works as a switch between 2d and 3d
%   all_angles          All Angles used in the backproj phi and the or euler angles 
%   thickness           Sample Thickness
%   angle_proj          parameter for 2d angle must be an element of ALL Angles
%  
%  OUTPUT
%   w_func              weighting function 2d or 3d
%
%EXAMPLE
%   2d:
%   w_func=tom_calc_weight_functionc([64 64],[10 20; 30 40; 50 10],64,[10 20]);  
%   3d:
%   w_func=tom_calc_weight_functionc([64 64 64],[10 20; 30 40; 50 10],64);
%
%NOTE:
%   All input types are double. Returned w_func is single.
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
%   created by SN and FB after a Fortran routine by RH 11/27/05
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

if (nargin==3)
       angle_proj=[0 0];
end;

%check inputs
if  ((size(dims,2)==2) & (nargin==3))
    error('projection angle needed for 2d weighting');
end;

if  ((size(dims,2)==3) & (nargin==4))
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


w_func=tom_calc_weight_functioninc(dims,all_angles',thickness,angle_proj);
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



















