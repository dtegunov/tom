function out=tom_weight3d(varargin)
%TOM_WEIGHT3D: Weighting function for tomography
%
%   out=tom_weight3d(varargin)
%
%PARAMETERS
%
%  INPUT
%   type                'exact' or 'analytical' (default analytical)
%   fim                 Input image in FOURIER space
%   p                   Radius for the low pass filter (in pixel)
%   psi                 tiltaxis direction (only for exact weighting)
%   num                 index of projection to be weighted (only for exact weighting)
%   theta               tilt angle as a vector (dim = ntilt) 
%   D                   Specify the thickeness of the Object in pixels  only for exact weighting) (Uniform sampling in x-y-z is assumed)
%  
%  OUTPUT
%   out                 weighted input image in FOURIER space
%
%DESCRIPTION:
%  tom_weight performs analytical (r*) or exact weighting of individual
%  micrographs for subsequent 3D reconstruction
%
%  NOTE:
%  for analytical type its possible to import the weighting funktion and the
%  circlemask ... to speed up the process
%  Sytax:out=tom_weight3d(type,p,fim,'w_func',w_func,'c_mask',c_mask)
%  w_func: weighting function
%  c_mask: circle mask for Filtering
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
%  fimage=tom_weight3d('analytical',fimage,filter);  
%  fimage=tom_weight3d('analytical',fimage,filter,'w_func',w_func,'c_mask',c_mask);  
%  fimage=tom_weight3d('exact',fimage,filter,psi,itilt,theta_ang,thickness);
%
%REFERENCES
%
%SEE ALSO
%   tom_reconstruction3d, tom_dist.c, tom_rec3d
%
%   created by AF 09/11/02
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


% parse inputs 
if (nargin < 2)
    disp('Error at least 2 parameters required: out=tom_weight3d(type,fim)'); 
    out=-1; return;
end;
if (strcmp(varargin{1},'analytical'))
    type=varargin{1};
    fim=varargin{2};
    p=0;
    w_func=0;
    c_mask=0;
    if (nargin > 2)
        p=varargin{3};
    end;
    if (nargin > 3 && strcmp(varargin{4},'w_func'))
        w_func=varargin{5};   
    end;
    if (nargin > 6 && strcmp(varargin{6},'c_mask'))
        c_mask=varargin{7};   
    end;
    if (nargin > 6 && strcmp(varargin{6},'w_fun'))
        w_func=varargin{7};   
    end; 
    if (nargin > 3 && strcmp(varargin{4},'c_mask'))
        c_mask=varargin{5};   
    end; 
end;
if (strcmp(varargin{1},'exact'))
    if (nargin < 7)
        disp('7 parameters required for exact weighting');    
        out=-1; return;
    end;
    type=varargin{1};
    fim=varargin{2};
    p=varargin{3};
    psi=varargin{4};
    num=varargin{5};
    theta=varargin{6};
    D=varargin{7};
end;

if (w_func ~= 0) || (c_mask ~= 0 || strmatch(type,'exact') )
    [s1,s2]=size(fim);
    %MeshGrid with the sampling points of the image
    [x,y]=ndgrid(-s1/2:s1/2-1,-s2/2:s2/2-1);
end;

%Analytical weighting and lowpass
if strmatch('analytical',type)
    if (w_func~=0)
        w_func=tom_norm(abs(x),1);    
    end;
    if (c_mask~=0 && p~=0)
        c_mask=tom_circle(p,s1);
    end;
    
    if ( c_mask ~=0 )
        out= (fim).*w_func.*c_mask;
    else
        out =(fim).*tom_norm(abs(x),1);
    end;
end;

if strmatch(type,'exact') 
    N=size(theta,1);
    D=s1/D*2;
    
    %Rotation from phi ton 90??? and y=0
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
end;

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

