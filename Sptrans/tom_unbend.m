function [x,y,struct_out]=tom_unbend(image,size_cut,x,y)
%TOM_UNBEND creates ...
%
%   [x,y,struct_out]=tom_unbend(image,size_cut,x,y)
%
%PARAMETERS
%
%  INPUT
%   image               ...
%   size_cut            ...
%   x                   ...
%   y                   ...
%  
%  OUTPUT
%   x                   ...
%   y                   ...
%   struct_out          ...
%
%EXAMPLE
%   ... = tom_unbend(...);
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


%get the samples from user
%figure; tom_imagesc(image);
%[x,y]=ginput;
%save stack_bend_4st;
%load stack_bend_4st;
%clear size;

if (nargin==2)
    %figure; tom_imagesc(image);
    kk=1;
    set(gcf,'doublebuffer','on');
    butt=1;
    while (butt==1)
        [xh,yh,butt]=ginput(1);
        xi(kk)=xh;
        yi(kk)=yh;
        %debug staff
        hold on;    
        plot(xi,yi,'r+');
        hold off;
        kk=kk+1;
        %debug staff end
    end;
    kk=kk-1;
    x(1:kk)=xi(1:kk);
    y(1:kk)=yi(1:kk);
end;

if (nargin==3)
    disp('Hello Mr. user you made a mistake only 2 or 4 parameters are allowed')
    disp('image, size to cut out, x-coordinats, y-coordinats ... have fun!')
    return;
end;


%estimate length of the structure
length=(x(size(x,2))-x(1))*5;

%allocate some momory
struct=zeros(length,size_cut);
num_struct=1;
for i=1:(size(x,2)-1) %number of samples 
    
    %calculate dx/dy and the angle of the linear equation 
    m=-(y(i+1)-y(i))/(x(i+1)-x(i));
    angle= atan(m)*180./pi;
    angle_rad=pi./180*angle;
    
    % calculate diff v for new cosy 
    v=[(x(i+1)-x(i)) (y(i+1)-y(i))]; 
    %calculat sampling distance
    dist=round((x(i+1)-x(i)));
    
    %cut out larger area
    im_box=image((round(x(i))-2*dist)+1:(round(x(i))+2*dist),(round(y(i))-2*size_cut)+1:(round(y(i))+2*size_cut) );
    middle=[ (round(size(im_box,1)./2)) (round(size(im_box,2)./2))];
    %rotate the area
    angle;
    rot_area=imrotate(im_box,angle,'bilinear','crop');
    % calculate rotation matrix
    rot_mat=[ cos(angle_rad) sin(angle_rad) ; -sin(angle_rad) cos(angle_rad) ];
    % tranform coordinates
    p_st=middle+v;
    p_rot=(v*rot_mat)+middle;
    
    %cut out the region for unbending
    str_box=rot_area(round(middle(1))+1:round(p_rot(1)),(middle(2)-round(size_cut./2))+1:(round(middle(2))+round(size_cut./2)) );
    
    %save the staff to structure
    struct(num_struct:(num_struct+size(str_box,1)-1),1:size_cut)=str_box(1:size(str_box,1),1:size(str_box,2));
    num_struct=num_struct+size(str_box,1);
end;

% if you' re still alive joy the final box out
struct_out(1:(num_struct-1),1:size_cut)=struct(1:(num_struct-1),1:size_cut);



