function c=tom_move(varargin)
%TOM MOVE Moves an image in X,Y,Z direction and fills the rest with zeros
%
%   c=tom_move(varargin)
%
%   C=TOM_MOVE(A,[dx dy dz]) This function moves image A in X,Y and Z direction. the
%   coordinates of the first pixel is now (dx+1,dy+1,dz+1). The rest of the image
%   is filled with zeros
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   c                   ...
%
%EXAMPLE
%       c=tom_move(a,[3 3])
%
%               1   2   3   4   5   6   7               
%               8   9   10  11  12  13  14              0   0   0   0   0   0   0
%               15  16  17  18  19  20  21              0   0   0   0   0   0   0
%           a = 22  23  24  25  26  27  28          c = 0   0   0   0   0   0   0
%               29  30  31  32  33  34  35              0   0   0   1   2   3   4
%               36  37  38  39  40  41  42              0   0   0   8   9   10  11
%               43  44  45  46  47  48  49              0   0   0   15  16  17  18
%
%       c=tom_move(a,[-3 -3])
%       
%              25  26  27  28   0   0   0
%              32  33  34  35   0   0   0
%          c = 39  40  41  42   0   0   0
%              46  47  48  49   0   0   0
%              0   0   0   0    0   0   0
%              0   0   0   0    0   0   0
%
%REFERENCES
%
%SEE ALSO
%   TOM_SHIFT
%
%   created by AL 08/13/02
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


if nargin<2
    error('Not Enough Input Arguments');
    return;
end
if nargin>2
    error('Too Many Input Arguments');
    return;
end
a=varargin{1};
coord=varargin{2};
[s1 s2 s3]=size(a);
% 2-D Image Case
if s3 == 1    
    if abs(coord(1))>=s1 | abs(coord(2))>=s2 % Out of limits move
        error('Wrong Input Arguments');
        return;
    end
    b=zeros(s1,s2);
    if coord(1)>=0 & coord(2)>=0   
        atmp=a(1:(s1-coord(1)),1:(s2-coord(2)));   
        b((coord(1)+1):s1,(coord(2)+1):s2)=atmp;
        c=b;
    elseif coord(1)<0 & coord(2)>=0
        atmp=a((abs(coord(1))+1):s1,1:(s2-coord(2)));   
        satmp=size(atmp);
        b(1:satmp(1),(coord(2)+1):s2)=atmp;
        c=b;    
    elseif coord(1)>=0 & coord(2)<0
        atmp=a(1:(s1-coord(1)),(abs(coord(2))+1):s2);    
        satmp=size(atmp);
        b((coord(1)+1):s1,1:satmp(2))=atmp;
        c=b;
    else
        atmp=a((abs(coord(1))+1):s1,(abs(coord(2))+1):s2);
        satmp=size(atmp);
        b(1:satmp(1),1:satmp(2))=atmp;
        c=b;
    end
else
    % 3-D Image Case
    if abs(coord(1))>=s1 | abs(coord(2))>=s2 | abs(coord(3))>=s3 % Out of limits move
        error('Wrong Input Arguments');
        return;
    end
    c=zeros(s1,s2,s3);
    b=zeros(s1,s2);
    if coord(3)>=0
        for i=1:(s3-coord(3))
            if coord(1)>=0 & coord(2)>=0   
                atmp=a(1:(s1-coord(1)),1:(s2-coord(2)),i);   
                b((coord(1)+1):s1,(coord(2)+1):s2)=atmp;
                c(:,:,(i+coord(3)))=b;
            elseif coord(1)<0 & coord(2)>=0
                atmp=a((abs(coord(1))+1):s1,1:(s2-coord(2)),i);   
                satmp=size(atmp);
                b(1:satmp(1),(coord(2)+1):s2)=atmp;
                c(:,:,(i+coord(3)))=b;    
            elseif coord(1)>=0 & coord(2)<0
                atmp=a(1:(s1-coord(1)),(abs(coord(2))+1):s2,i);    
                satmp=size(atmp);
                b((coord(1)+1):s1,1:satmp(2))=atmp;
                c(:,:,(i+coord(3)))=b;
            else
                atmp=a((abs(coord(1))+1):s1,(abs(coord(2))+1):s2,i);
                satmp=size(atmp);
                b(1:satmp(1),1:satmp(2))=atmp;
                c(:,:,(i+coord(3)))=b;
            end
        end
    else
        for i=(abs(coord(3))+1):s3
            if coord(1)>=0 & coord(2)>=0   
                atmp=a(1:(s1-coord(1)),1:(s2-coord(2)),i);   
                b((coord(1)+1):s1,(coord(2)+1):s2)=atmp;
                c(:,:,(i-abs(coord(3))))=b;
            elseif coord(1)<0 & coord(2)>=0
                atmp=a((abs(coord(1))+1):s1,1:(s2-coord(2)),i);   
                satmp=size(atmp);
                b(1:satmp(1),(coord(2)+1):s2)=atmp;
                c(:,:,(i-abs(coord(3))))=b;    
            elseif coord(1)>=0 & coord(2)<0
                atmp=a(1:(s1-coord(1)),(abs(coord(2))+1):s2,i);    
                satmp=size(atmp);
                b((coord(1)+1):s1,1:satmp(2))=atmp;
                c(:,:,(i-abs(coord(3))))=b;
            else
                atmp=a((abs(coord(1))+1):s1,(abs(coord(2))+1):s2,i);
                satmp=size(atmp);
                b(1:satmp(1),1:satmp(2))=atmp;
                c(:,:,(i-abs(coord(3))))=b;
            end            
        end
    end
end
