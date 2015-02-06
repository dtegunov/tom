function c=tom_red(varargin)
%TOM_RED Extracts a subimage/subvolume from tha input image/volume
%
%   c=tom_red(varargin)
%
%   C=TOM_RED(A,COORD,DIMS) This functions extracts a subimage C from the input image A. The COORD
%   variable refers to the coordinates of the starting pixel (upper left pixel) of the subimage
%   and the DIMS variable refers to the dimensions of the subimage. If the dimensions of the
%   subimage exceed the dimensions of the original image A then the rest of the subimage is
%   filled with zeros. This function is also compatible with 3-D Volumes. In that case a 
%   subvolume is extracted ftom the original one. The COORD variable can also take negative
%   values
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
%            0    17   -10    -3
%      A = -10     6     0    11
%            7    -6     0   -18
%            6     4     1     5
%
%   C=TOM_RED(A,[2 2],[3 3]); A subimage is going to be created and stored in the array C.
%   Its elements will start from the (2,2) pixel of image A and its dimensions will be
%   3x3.
%                            6     0    11
%                      c  = -6     0   -18
%                            4     1     5
%   C=TOM_RED(A,[0 3],[2 2]); A subimage is going to be created and stored in the array C.
%   Its elements will start from the (0 3) pixel of image A and its dimensions will be
%   2x2
%                   c =  0   0
%                      -10  -3
%   C=TOM_RED(A,[-1 -1],[4 4]);
%   
%                        0     0     0     0
%                   c =  0     0     0     0
%                        0     0     0    17
%                        0     0   -10     6
%
%REFERENCES
%
%SEE ALSO
%   TOM_PASTE, TOM_MOVE, TOM_MIRROR
%
%   created by AL 08/17/02
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


if nargin<3
    error('Not Enough Input Arguments');
elseif nargin == 3
    a=varargin{1};
    coord=varargin{2};
    dims=varargin{3};
    for i=1:size(dims(2))
        if dims(i)<=0
            error('Dimensions must be positive Integers');
            return;
        end
    end        
    [s1 s2 s3]=size(a);
    if s3 == 1
       if  coord(1)>s1 | coord(2)>s2
           error('Inproper Selection of Starting Pixel ');
       elseif coord(1)<=0 | coord(2)<=0
               if coord(1)<=0 & coord(2)>0
                   if coord(1)+dims(1)-1<1
                       error('Inproper Selection of Dimensions');
                   elseif (coord(2)+dims(2)-1)>s2
                       xtra=coord(2)+dims(2)-1-s2;
                       b=a(1:(coord(1)+dims(1)-1),coord(2):s2);
                       dimb=size(b);
                       b=[zeros((abs(coord(1))+1),dimb(2));b];
                       dimb=size(b);
                       b=[b zeros(dimb(1),xtra)];
                       c=b;
                   else                           
                       b=a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1));
                       dimb=size(b);
                       b=[zeros((abs(coord(1))+1),dimb(2));b];
                       c=b;
                   end
               elseif coord(1)>0 & coord(2)<=0
                   if coord(2)+dims(2)-1<1
                       error('Inproper Selection of Dimensions');
                   elseif (coord(1)+dims(1)-1)>s1
                       xtra=coord(1)+dims(1)-1-s1;                                               
                       b=a(coord(1):s1,1:coord(2)+dims(2)-1);                       
                       dimb=size(b);
                       b=[zeros(dimb(1),abs(coord(2))+1) b];
                       dimb=size(b);
                       b=[b;zeros(xtra,dimb(2))]; 
                       c=b;                     
                   else
                       b=a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1));
                       dimb=size(b);
                       b=[zeros(dimb(1),(abs(coord(2))+1)) b];
                       c=b;
                   end
               elseif coord(1)<=0 & coord(2)<=0
                   if coord(1)+dims(1)-1<1 | coord(2)+dims(2)-1<1
                       error('Inproper Selection of Dimensions');
                   else
                       b=a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1));
                       dimb1=size(b);
                       b=[zeros(dimb1(1),(abs(coord(2))+1)) b];
                       dimb2=size(b);
                       b=[zeros((abs(coord(1))+1),dimb2(2));b];
                       c=b;  
                   end
               end
                       
       else                     
                if (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)<=s2
                    c=a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1));
                elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)<=s2
                    xtra=coord(1)+dims(1)-1-s1;
                    b=a(coord(1):s1,coord(2):coord(2)+dims(2)-1);
                    dimb=size(b);
                    c=[b;zeros(xtra,dimb(2))];        
                elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)>s2
                    xtra=coord(2)+dims(2)-1-s2;
                    b=a(coord(1):(coord(1)+dims(1)-1),coord(2):s2);
                    dimb=size(b);
                    b=[b zeros(dimb(1),xtra)];
                    c=b;
                elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)>s2
                    xtra1=coord(1)+dims(1)-1-s1;
                    xtra2=coord(2)+dims(2)-1-s2;
                    b=a(coord(1):s1,coord(2):s2);
                    dimb1=size(b);
                    b=[b zeros(dimb1(1),xtra2)];
                    dimb2=size(b);
                    b=[b;zeros(xtra2,dimb2(2))];
                    c=b;
                end                                
       end
       
       
       
       
       
       
       
       
       
       
       
       
       
       
    else
        
        
   %----3-D Images--------     
        
        
        
        
       
        
        
       if  coord(1)>s1 | coord(2)>s2 | coord(3)>s3
           error('Inproper Selection of Starting Pixel ');
           return;
          
       elseif coord(1)<=0 | coord(2)<=0 | coord(3)<=0
           if (coord(1)+dims(1)-1)<1 | (coord(2)+dims(2)-1)<1 | (coord(3)+dims(3)-1)<1
               error('Inproper Selection Starting Pixel');
               return;
           end
       end
       % Over the volume 
       if coord(3)<=0 & ((coord(3)+dims(3)-1)<=s3)
           ad=abs(coord(3))+1;
           for i=coord(3):(coord(3)+dims(3)-1)
               if i<=0
                   c(:,:,(i+ad))=zeros(dims(1),dims(2));
               else                 
                   if coord(1)<=0 & coord(2)>0     
                        if (coord(2)+dims(2)-1)>s2
                            xtra=coord(2)+dims(2)-1-s2;
                            b=a(1:(coord(1)+dims(1)-1),coord(2):s2,i);
                            dimb=size(b);
                            b=[zeros((abs(coord(1))+1),dimb(2));b];
                            dimb=size(b);
                            b=[b zeros(dimb(1),xtra)];
                            c(:,:,(i+ad))=b;
                        else
                            b=a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i);
                            dimb=size(b);
                            b=[zeros((abs(coord(1))+1),dimb(2));b];
                            c(:,:,(i+ad))=b;
                        end                      
                   elseif coord(1)>0 & coord(2)<=0              
                        if (coord(1)+dims(1)-1)>s1
                            xtra=coord(1)+dims(1)-1-s1;                                               
                            b=a(coord(1):s1,1:coord(2)+dims(2)-1,i);                       
                            dimb=size(b);
                            b=[zeros(dimb(1),abs(coord(2))+1) b];
                            dimb=size(b);
                            b=[b;zeros(xtra,dimb(2))]; 
                            c(:,:,(i+ad))=b;
                        else                                      
                            b=a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i);
                            dimb=size(b);
                            b=[zeros(dimb(1),(abs(coord(2))+1)) b];
                            c(:,:,(i+ad))=b;
                        end              
                  elseif coord(1)<=0 & coord(2)<=0                  
                            b=a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i);
                            dimb1=size(b);
                            b=[zeros(dimb1(1),(abs(coord(2))+1)) b];
                            dimb2=size(b);
                            b=[zeros((abs(coord(1))+1),dimb2(2));b];
                            c(:,:,(i+ad))=b;            
                  elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)<=s2              
                            c(:,:,(i+ad))=a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i);              
                  elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)<=s2               
                            xtra=coord(1)+dims(1)-1-s1;
                            b=a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i);
                            dimb=size(b);
                            c(:,:,(i+ad))=[b;zeros(xtra,dimb(2))];               
                  elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)>s2               
                            xtra=coord(2)+dims(2)-1-s2;
                            b=a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i);
                            dimb=size(b);
                            b=[b zeros(dimb(1),xtra)];
                            c(:,:,(i+ad))=b;               
                  elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)>s2              
                            xtra1=coord(1)+dims(1)-1-s1;
                            xtra2=coord(2)+dims(2)-1-s2;
                            b=a(coord(1):s1,coord(2):s2,i);
                            dimb1=size(b);
                            b=[b zeros(dimb1(1),xtra2)];
                            dimb2=size(b);
                            b=[b;zeros(xtra2,dimb2(2))];
                            c(:,:,(i+ad))=b;
              
                  end
              end
           end

          
          
       % Inside the volume            
       elseif coord(3)>=1  
           ad=coord(3)-1;
           for i=coord(3):(coord(3)+dims(3)-1)
               if i>s3
                   c(:,:,(i-ad))=zeros(dims(1),dims(2));
               else                   
                   if coord(1)<=0 & coord(2)>0               
                         if (coord(2)+dims(2)-1)>s2
                            xtra=coord(2)+dims(2)-1-s2;
                            b=a(1:(coord(1)+dims(1)-1),coord(2):s2,i);
                            dimb=size(b);
                            b=[zeros((abs(coord(1))+1),dimb(2));b];
                            dimb=size(b);
                            b=[b zeros(dimb(1),xtra)];
                            c(:,:,(i-ad))=b;
                         else
                            b=a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i);
                            dimb=size(b);
                            b=[zeros((abs(coord(1))+1),dimb(2));b];
                            c(:,:,(i-ad))=b;
                         end             
                   elseif coord(1)>0 & coord(2)<=0               
                        if (coord(1)+dims(1)-1)>s1
                             xtra=coord(1)+dims(1)-1-s1;                                               
                             b=a(coord(1):s1,1:coord(2)+dims(2)-1);                       
                             dimb=size(b);
                             b=[zeros(dimb(1),abs(coord(2))+1) b];
                             dimb=size(b);
                             b=[b;zeros(xtra,dimb(2))]; 
                            c(:,:,(i-ad))=b;
                        else                                      
                            b=a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i);
                            dimb=size(b);
                            b=[zeros(dimb(1),(abs(coord(2))+1)) b];
                            c(:,:,(i-ad))=b;
                        end               
                   elseif coord(1)<=0 & coord(2)<=0               
                        b=a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i);
                        dimb1=size(b);
                        b=[zeros(dimb1(1),(abs(coord(2))+1)) b];
                        dimb2=size(b);
                        b=[zeros((abs(coord(1))+1),dimb2(2));b];
                        c(:,:,(i-ad))=b;              
                   elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)<=s2               
                        c(:,:,(i-ad))=a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i);              
                   elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)<=s2               
                        xtra=coord(1)+dims(1)-1-s1;
                        b=a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i);
                        dimb=size(b);
                        c(:,:,(i-ad))=[b;zeros(xtra,dimb(2))];               
                   elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)>s2              
                        xtra=coord(2)+dims(2)-1-s2;
                        b=a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i);
                        dimb=size(b);
                        b=[b zeros(dimb(1),xtra)];
                        c(:,:,(i-ad))=b;              
                   elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)>s2               
                        xtra1=coord(1)+dims(1)-1-s1;
                        xtra2=coord(2)+dims(2)-1-s2;
                        b=a(coord(1):s1,coord(2):s2,i);
                        dimb1=size(b);
                        b=[b zeros(dimb1(1),xtra2)];
                        dimb2=size(b);
                        b=[b;zeros(xtra2,dimb2(2))];
                        c(:,:,(i-ad))=b;               
                   end
               end
           end
       end
       
       
    %new
    
    
    if coord(3)<=0 && ((coord(3)+dims(3)-1)>s3)
           ad=abs(coord(3))+1;
           for i=coord(3):(coord(3)+dims(3)-1)
               if i<=0
                   c(:,:,(i+ad))=zeros(dims(1),dims(2));
               end;
               
               if i <= size(a,3)  && i > 0                
                   if coord(1)<=0 & coord(2)>0     
                        if (coord(2)+dims(2)-1)>s2
                            xtra=coord(2)+dims(2)-1-s2;
                            b=a(1:(coord(1)+dims(1)-1),coord(2):s2,i);
                            dimb=size(b);
                            b=[zeros((abs(coord(1))+1),dimb(2));b];
                            dimb=size(b);
                            b=[b zeros(dimb(1),xtra)];
                            c(:,:,(i+ad))=b;
                        else
                            b=a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i);
                            dimb=size(b);
                            b=[zeros((abs(coord(1))+1),dimb(2));b];
                            c(:,:,(i+ad))=b;
                        end                      
                   elseif coord(1)>0 & coord(2)<=0              
                        if (coord(1)+dims(1)-1)>s1
                            xtra=coord(1)+dims(1)-1-s1;                                               
                            b=a(coord(1):s1,1:coord(2)+dims(2)-1,i);                       
                            dimb=size(b);
                            b=[zeros(dimb(1),abs(coord(2))+1) b];
                            dimb=size(b);
                            b=[b;zeros(xtra,dimb(2))]; 
                            c(:,:,(i+ad))=b;
                        else                                      
                            b=a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i);
                            dimb=size(b);
                            b=[zeros(dimb(1),(abs(coord(2))+1)) b];
                            c(:,:,(i+ad))=b;
                        end              
                  elseif coord(1)<=0 & coord(2)<=0                  
                            b=a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i);
                            dimb1=size(b);
                            b=[zeros(dimb1(1),(abs(coord(2))+1)) b];
                            dimb2=size(b);
                            b=[zeros((abs(coord(1))+1),dimb2(2));b];
                            c(:,:,(i+ad))=b;            
                  elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)<=s2              
                            c(:,:,(i+ad))=a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i);              
                  elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)<=s2               
                            xtra=coord(1)+dims(1)-1-s1;
                            b=a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i);
                            dimb=size(b);
                            c(:,:,(i+ad))=[b;zeros(xtra,dimb(2))];               
                  elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)>s2               
                            xtra=coord(2)+dims(2)-1-s2;
                            b=a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i);
                            dimb=size(b);
                            b=[b zeros(dimb(1),xtra)];
                            c(:,:,(i+ad))=b;               
                  elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)>s2              
                            xtra1=coord(1)+dims(1)-1-s1;
                            xtra2=coord(2)+dims(2)-1-s2;
                            b=a(coord(1):s1,coord(2):s2,i);
                            dimb1=size(b);
                            b=[b zeros(dimb1(1),xtra2)];
                            dimb2=size(b);
                            b=[b;zeros(xtra2,dimb2(2))];
                            c(:,:,(i+ad))=b;
              
                  end
               end
               
               if (i > size(a,3)) 
                    c(:,:,i)=zeros(dims(1),dims(2));
               end;
               
           end
           
           
    end;
    
    %end new   
       
       
    end
      
    
    
else
    error('Too Many Input Arguments');
end

