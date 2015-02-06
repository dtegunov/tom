function Iout=affine_transform_2d_double(Iin,M,ImageSize,check_bil_intp)
% Affine transformation function (Rotation, Translation, Resize)
% This function transforms a volume with a 3x3 transformation matrix 
%
% Iout=affine_transform_2d_double(Iin,Minv,ImageSize,check_bil_intp)
%
% inputs,
%   Iin: The greyscale input image
%   Minv: The (inverse) 3x3 transformation matrix
%   ImageSize: Size of output imgae
%   check_bil_intp: If true bilinear interpolation otherwise nearest neigb.
% output,
%   Iout: The transformed image
%
% example,
%   % Read image
%   I=im2double(imread('lenag2.png'))
%   % Make a transformation matrix
%   M=make_transformation_matrix([2 3],[1.0 1.1],2);
%   % Transform the image
%   Iout=rigid_transform_2d_double(I,M,size(I))
%   % Show the image
%   figure, imshow(Iout);
%
% Function is written by D.Kroon University of Twente (February 2009)
  
% Make all x,y indices
[x,y]=ndgrid(0:ImageSize(1)-1,0:ImageSize(2)-1);

% Calculate center of the output image
mean_out=ImageSize/2;

% Calculate center of the input image
mean_in=size(Iin)/2;

% Make center of the image coordinates 0,0
xd=x-mean_out(1); 
yd=y-mean_out(2);

% Calculate the Transformed coordinates
Tlocalx = mean_in(1) + M(1,1) * xd + M(1,2) *yd + M(1,3) * 1;
Tlocaly = mean_in(2) + M(2,1) * xd + M(2,2) *yd + M(2,3) * 1;

% All the neighborh pixels involved in linear interpolation.
if(check_bil_intp)
    xBas0=floor(Tlocalx); 
    yBas0=floor(Tlocaly);
    xBas1=xBas0+1;           
    yBas1=yBas0+1;

    % Linear interpolation constants (percentages)
    xCom=Tlocalx-xBas0; 
    yCom=Tlocaly-yBas0;
    perc0=(1-xCom).*(1-yCom);
    perc1=(1-xCom).*yCom;
    perc2=xCom.*(1-yCom);
    perc3=xCom.*yCom;
else
    xBas0=round(Tlocalx); 
    yBas0=round(Tlocaly);
end

% limit indexes to boundaries
check_xBas0=(xBas0<0)|(xBas0>(size(Iin,1)-1));
check_yBas0=(yBas0<0)|(yBas0>(size(Iin,2)-1));
xBas0(check_xBas0)=0; 
yBas0(check_yBas0)=0; 

if(check_bil_intp)
    % limit indexes to boundaries
    check_xBas1=(xBas1<0)|(xBas1>(size(Iin,1)-1));
    check_yBas1=(yBas1<0)|(yBas1>(size(Iin,2)-1));
    xBas1(check_xBas1)=0; 
    yBas1(check_yBas1)=0; 
end

if(ndims(Iin)==2), 
    lo=1; 
    Iout=zeros([ImageSize 1]); 
else
    lo=3; 
    Iout=zeros([ImageSize 3]); 
end
    
for i=1:lo;    
    Iin_one=Iin(:,:,i);
    if(check_bil_intp)
        % Get the intensities
        intensity_xyz0=Iin_one(1+xBas0+yBas0*size(Iin,1));
        intensity_xyz1=Iin_one(1+xBas0+yBas1*size(Iin,1)); 
        intensity_xyz2=Iin_one(1+xBas1+yBas0*size(Iin,1));
        intensity_xyz3=Iin_one(1+xBas1+yBas1*size(Iin,1));
        % Make pixels before outside Ibuffer black
        intensity_xyz0(check_xBas0|check_yBas0)=0;
        intensity_xyz1(check_xBas0|check_yBas1)=0;
        intensity_xyz2(check_xBas1|check_yBas0)=0;
        intensity_xyz3(check_xBas1|check_yBas1)=0;
        Iout_one=intensity_xyz0.*perc0+intensity_xyz1.*perc1+intensity_xyz2.*perc2+intensity_xyz3.*perc3;
    else
        % Get the intensities
        intensity_xyz0=Iin_one(1+xBas0+yBas0*size(Iin,1));
        % Make pixels before outside Ibuffer black
        intensity_xyz0(check_xBas0|check_yBas0)=0;
        Iout_one=intensity_xyz0;
    end
    
    Iout(:,:,i)=reshape(Iout_one, ImageSize);
end