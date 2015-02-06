function out=tom_weight3d_euler(varargin)
%TOM_WEIGHT3D_EULER: weighting function for tomography
%
%   out=tom_weight3d(type,fim,p,psi,num,theta, D)
%
%PARAMETERS
%
%  INPUT
%   proj                Input image in REAL space
%   phi                 euler 1
%   psi                 euler 2              
%   theta               euler 3 all tiltangles
%   D                   Specify the thickeness of the Object in pixels only for exact weighting) (Uniform sampling in x-y-z is assumed)

%  
%  OUTPUT
%   out                 weighted input image in REAL space
%
%   tom_weight performs exact weighting of individual
%   micrographs for subsequent 3D reconstruction
%
%NOTE:
%  Transformation 
%  1: Rotation of projeciton ,90-phi_in
%  2: Calculate anglels -psi_in-90
% 
%  Parameters num and following are only needed for exact weighting.
%
%  Theoretical background: Has to be applied on the 2D images prior to the 3D 
%  reconstruction in order to get the optimal result. Analytical weighting: Multiply
%  with a ramp filter parallel to the y-Axis all input images. Exact weighting: Multiply 
%  with the exact function (inverse of the contribution of each image on a certain    
%  sampling point in the Fourier space). Can be applied also for random conical tilting
%  or angular reconstitution. Analytical weighting is significantly faster but can not 
%  be applied on an no equidistant tilt sampling (e.g. Saxton scheme).
%
%  All angles correspond to the Euler angles.
%
%  The input image should be a square image (Fill up with zeros if necessary)
%
%EXAMPLE
%  proj=tom_weight3d_euler(proj,size(vol,1)./2-2,eu(1),eu(2),log(:,5),i,thickn);
%
%REFERENCES
%
%SEE ALSO
%   tom_reconstruction3d, tom_dist.c, tom_rec3d
%
%   created by FF 19/02/04
%   last change FB ???? messed up!
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


% switch to speed up the process by using a C-implementation for a
% time-consuming loop
USE_C_CODE=1;



if (nargin < 7)
    disp('7 parameters required for exact weighting');
    out=-1; return;
end;

proj=varargin{1};
p=varargin{2};
phi_in=varargin{3};
psi_in=varargin{4};
theta=varargin{5};
num=varargin{6};
D=varargin{7};


% transform to use the weighting function with 2 angles
psi=single(90-phi_in);
proj=single(tom_rotate(single(proj),-psi_in-90,'linear'));
fim = fftshift(fft2(double(proj)));



[s1,s2]=size(fim);
%MeshGrid with the sampling points of the image
[x,y]=ndgrid(-s1/2:s1/2-1,-s2/2:s2/2-1);


N=size(theta,1);
D=s1/D*2;

%Rotation from phi ton 90� and y=0
x_st=-x.*cos(theta(num)*pi/180).*sin(psi*pi/180)-y.*cos(psi*pi/180);
y_st= x.*cos(theta(num)*pi/180).*cos(psi*pi/180)-y.*sin(psi*pi/180);
z_st= x.*sin(theta(num)*pi/180);

w=(zeros(size(fim)));

if (USE_C_CODE==1)
    w=single(w);
    tom_dist(single(theta),single(num),single(psi),single(D),single(x_st),single(y_st),single(z_st),w);
    w=double(w);
    % the C-function tom_dist performs the commented calculation below
    % coded in C ...to speed up the process
else
    %w additional weighting due to neighbouring slices in Fourier space
    for kk=1:N
        if kk~=num
            % calculate the Distance to the plane
            d_tmp=+x_st.*sin(theta(kk)*pi/180).*sin(psi*pi/180)-y_st.*sin(theta(kk)*pi/180).*cos(psi*pi/180)+z_st.*cos(theta(kk)*pi/180);
            ind=find(abs(d_tmp)>D);
            d_tmp(ind)=D;
            w=w+tom_sinc(d_tmp./D);
        end;
    end;
    % normalize 0<=x<=1
end;

w=1./(1+w);
minw=min(min(w));
maxw=max(max(w));
w=(w-minw)./(maxw-minw);

%apply filter
if p~=0
    f=tom_circle(floor(p),s1);
    w=w.*f;
end;

%apply weighting
out=fim.*w;

%back to real space
out = real(ifft2(ifftshift(out)));

%rotate back
out=single(tom_rotate(single(out),-(-psi_in-90),'linear'));



function im = tom_circle(radius,dim);
[x,y] = ndgrid(-dim/2:dim/2-1);
im = sqrt((x).^2 + (y).^2) <= radius;
return;

function y=tom_sinc(x)
%TOM_SINC Sin(pi*x)/(pi*x) function.
%   SINC(X) returns a matrix whose elements are the sinc of the elements
%   of X, i.e.
%        y = sin(pi*x)/(pi*x)    if x ~= 0
%          = 1                   if x == 0
%   where x is an element of the input matrix and y is the resultant
%   output element.
%
%   See also SQUARE, SIN, COS, CHIRP, DIRIC, GAUSPULS, PULSTRAN, RECTPULS,
%   and TRIPULS.

%   Author(s): T. Krauss, 1-14-93
%   Copyright 1988-2001 The MathWorks, Inc.
%       $Revision: 1.6 $  $Date: 2001/04/02 20:20:54 $

y=ones(size(x));
i=find(x);
h=pi*x(i);
k=sin(h);
y(i)=k./h;

