function a=tom_paste(varargin)
%TOM_PASTE paste one array in another one
%
%   a=tom_paste(varargin)
%   a=TOM_PASTE(A,B,COORD)
%
%   The array B will be pasted into the array A so that
%   its upper left corner will become pixel (COORD(1) COORD(2)) in the array A.
%   COORD(1), COORD(2) may be outside A. The scaling of the rows/columns is :
%   -n -(n-1) ...... -1 0 1 2 3 which means that the COORD(1), COORD(2) can also
%   take the 0 value. This function is compatible with 3-D arrays also. In that
%   case the paste command is applied in every slice of the 3-D image. The
%   user can use the extra option 'min' or 'max' with the result that the
%   new values will be pasted only if they are smaller or bigger than the
%   original data. 
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   a                   ...
%
%EXAMPLE
%           -4     0    -1   -13   -14     7   -10
%          -16     4     2     8     6    12    15             1  1  1  1  1
%            2     2    11    17    -3   -12    -8             1  1  1  1  1
%       a =  3    -1     1    -6     7     0     6         b = 1  1  1  1  1 
%          -11     8     0     9     9    -1     3             1  1  1  1  1
%           12    -5    -8    13     8   -16    -9             1  1  1  1  1
%           12    22     3   -15    13     3   -21
%
%   c=tom_paste(a,b,[3 3]);                             c=tom_paste(a,b,[0 -1])
%
%      -4     0    -1   -13   -14     7   -10         1     1     1   -13   -14     7   -10
%     -16     4     2     8     6    12    15         1     1     1     8     6    12    15
%       2     2     1     1     1     1     1         1     1     1    17    -3   -12    -8
%  c =  3    -1     1     1     1     1     1   c  =  1     1     1    -6     7     0     6
%     -11     8     1     1     1     1     1       -11     8     0     9     9    -1     3
%      12    -5     1     1     1     1     1        12    -5    -8    13     8   -16    -9
%      12    22     1     1     1     1     1        12    22     3   -15    13     3   -21
%
%   c=tom_paste(b,a,[2 2]);
% 
%         1     1     1     1     1
%         1    -4     0    -1   -13
%   c  =  1   -16     4     2     8
%         1     2     2    11    17
%         1     3    -1     1    -6
%
%REFERENCES
%
%SEE ALSO
%   TOM_MOVE, TOM_PEAK, TOM_RED
%
%   created by AL 08/15/02
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
    b=varargin{2};
    coord=varargin{3};
    [dims(1) dims(2) dims(3)]=size(b);  
    [s1 s2 s3]=size(a);
    if s3 == 1
       if  coord(1)>s1 | coord(2)>s2
           error('Inproper Selection of Starting Pixel ');
       elseif coord(1)<=0 | coord(2)<=0
               if coord(1)<=0 & coord(2)>0
                   if coord(1)+dims(1)-1<1
                       error('Inproper Selection of Starting Pixel');
                   elseif (coord(2)+dims(2)-1)>s2                                                                  
                       a(1:(coord(1)+dims(1)-1),coord(2):s2)=b((abs(coord(1))+2):dims(1),1:(s2-coord(2)+1));                     
                   else                                             
                       a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1))=b((abs(coord(1))+2):dims(1),1:dims(2));                             
                   end
               elseif coord(1)>0 & coord(2)<=0
                   if coord(2)+dims(2)-1<1
                       error('Inproper Selection of Starting Pixel');
                   elseif (coord(1)+dims(1)-1)>s1                                            
                       a(coord(1):s1,1:coord(2)+dims(2)-1)=b(1:(s1-coord(1)+1),(abs(coord(2))+2):dims(2));                                          
                   else                                             
                       a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1))=b(1:dims(1),(abs(coord(2))+2):dims(2));                       
                   end
               elseif coord(1)<=0 & coord(2)<=0
                   if coord(1)+dims(1)-1<1 | coord(2)+dims(2)-1<1
                       error('Inproper Selection of Starting Pixel');
                   else                                             
                       a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1))=b((abs(coord(1))+2):dims(1),(abs(coord(2))+2):dims(2));                        
                   end
               end
                       
       else                     
                if (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)<=s2                   
                    a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1))=b;                   
                elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)<=s2                                       
                    a(coord(1):s1,coord(2):(coord(2)+dims(2)-1))=b(1:(1+s1-coord(1)),1:dims(2));                          
                elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)>s2                   
                    a(coord(1):(coord(1)+dims(1)-1),coord(2):s2)=b(1:dims(1),1:(s2-coord(2)+1))    ;                  
                elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)>s2                   
                    a(coord(1):s1,coord(2):s2)=b(1:(s1-coord(1)+1),1:(s2-coord(2)+1));                   
                end                                
        end                 
    else     
        
        
       if  coord(1)>s1 | coord(2)>s2 | coord(3)>s3
           error('Inproper Selection of Starting Pixel ');
           return;
          
       elseif coord(1)<=0 | coord(2)<=0 | coord(3)<=0
           if (coord(1)+dims(1)-1)<1 | (coord(2)+dims(2)-1)<1 | (coord(3)+dims(3)-1)<1
               error('Inproper Selection Starting Pixel');
               return;
           end
       end
      
       if coord(3)<=0 & ((coord(3)+dims(3)-1)<=s3)
           ad=abs(coord(3))+1;
           for i=coord(3):(coord(3)+dims(3)-1)
               if i<=0                 
               else                 
                   if coord(1)<=0 & coord(2)>0     
                        if (coord(2)+dims(2)-1)>s2                    
                            a(1:(coord(1)+dims(1)-1),coord(2):s2,i)=b((abs(coord(1))+2):dims(1),1:(s2-coord(2)+1),(i+ad));                           
                        else
                            a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i)=b((abs(coord(1))+2):dims(1),1:dims(2),(i+ad));                           
                        end                      
                   elseif coord(1)>0 & coord(2)<=0              
                        if (coord(1)+dims(1)-1)>s1                                                                          
                            a(coord(1):s1,1:coord(2)+dims(2)-1,i)=b(1:(s1-coord(1)+1),(abs(coord(2))+2):dims(2),(i+ad));                            
                        else                                      
                            a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i)=b(1:dims(1),(abs(coord(2))+2):dims(2),(i+ad));                        
                        end              
                  elseif coord(1)<=0 & coord(2)<=0                  
                            a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i)=b((abs(coord(1))+2):dims(1),(abs(coord(2))+2):dims(2),(i+ad));                                      
                  elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)<=s2              
                            a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i)=b(:,:,(i+ad));              
                  elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)<=s2                                          
                            a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i)=b(1:(1+s1-coord(1)),1:dims(2),(i+ad));                                          
                  elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)>s2                                         
                            a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i)=b(1:dims(1),1:(s2-coord(2)+1),(i+ad));                                           
                  elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)>s2                                         
                            a(coord(1):s1,coord(2):s2,i)=b(1:(s1-coord(1)+1),1:(s2-coord(2)+1),(i+ad));                                          
                  end
              end
          end
           
       elseif coord(3)>=1  
           ad=coord(3)-1;
           for i=coord(3):(coord(3)+dims(3)-1)
               if i>s3                   
               else                   
                   if coord(1)<=0 & coord(2)>0               
                         if (coord(2)+dims(2)-1)>s2                           
                            a(1:(coord(1)+dims(1)-1),coord(2):s2,i)=b((abs(coord(1))+2):dims(1),1:(s2-coord(2)+1),(i-ad));                           
                         else
                            a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i)=b((abs(coord(1))+2):dims(1),1:dims(2),(i-ad));                            
                         end             
                   elseif coord(1)>0 & coord(2)<=0               
                        if (coord(1)+dims(1)-1)>s1                                                                          
                             a(coord(1):s1,1:coord(2)+dims(2)-1,i)=b(1:(s1-coord(1)+1),(abs(coord(2))+2):dims(2),(i-ad));                                                    
                        else                                      
                            a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i)=b(1:dims(1),(abs(coord(2))+2):dims(2),(i-ad));                            
                        end               
                   elseif coord(1)<=0 & coord(2)<=0               
                        a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i)=b((abs(coord(1))+2):dims(1),(abs(coord(2))+2):dims(2),(i-ad));                                     
                   elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)<=s2               
                       a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i)=b(:,:,(i-ad));              
                   elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)<=s2                                      
                        a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i)=b(1:(1+s1-coord(1)),1:dims(2),(i-ad));                                    
                   elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)>s2                                    
                        a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i)=b(1:dims(1),1:(s2-coord(2)+1),(i-ad));                                
                   elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)>s2                                   
                        a(coord(1):s1,coord(2):s2,i)=b(1:(s1-coord(1)+1),1:(s2-coord(2)+1),(i-ad));                                   
                   end
               end
           end
       end       
   end
elseif nargin == 4
    a=varargin{1};
    b=varargin{2};
    coord=varargin{3};
    str=varargin{4};
    [dims(1) dims(2) dims(3)]=size(b);  
    [s1 s2 s3]=size(a);
    if s3 == 1
       if  coord(1)>s1 | coord(2)>s2
           error('Inproper Selection of Starting Pixel ');
       elseif coord(1)<=0 | coord(2)<=0
               if coord(1)<=0 & coord(2)>0
                   if coord(1)+dims(1)-1<1
                       error('Inproper Selection of Starting Pixel');
                   elseif (coord(2)+dims(2)-1)>s2 
                       if isequal(str,'max')                           
                            a(1:(coord(1)+dims(1)-1),coord(2):s2)=max(a(1:(coord(1)+dims(1)-1),coord(2):s2),b((abs(coord(1))+2):dims(1),1:(s2-coord(2)+1)));                     
                        elseif isequal(str,'min')
                            a(1:(coord(1)+dims(1)-1),coord(2):s2)=min(a(1:(coord(1)+dims(1)-1),coord(2):s2),b((abs(coord(1))+2):dims(1),1:(s2-coord(2)+1)));
                        else
                            error('Unknown parameter');
                        end
                   else
                       if isequal(str,'max')
                            a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1))=max(a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1)),b((abs(coord(1))+2):dims(1),1:dims(2)));                             
                       elseif isequal(str,'min')
                             a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1))=min(a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1)),b((abs(coord(1))+2):dims(1),1:dims(2)));        
                       else
                          error('Unknown parameter');
                       end
                   end
               elseif coord(1)>0 & coord(2)<=0
                   if coord(2)+dims(2)-1<1
                       error('Inproper Selection of Starting Pixel');
                   elseif (coord(1)+dims(1)-1)>s1   
                       if isequal(str,'max')
                            a(coord(1):s1,1:coord(2)+dims(2)-1)=max(a(coord(1):s1,1:coord(2)+dims(2)-1),b(1:(s1-coord(1)+1),(abs(coord(2))+2):dims(2)));                             
                       elseif isequal(str,'min')
                            a(coord(1):s1,1:coord(2)+dims(2)-1)=min(a(coord(1):s1,1:coord(2)+dims(2)-1),b(1:(s1-coord(1)+1),(abs(coord(2))+2):dims(2)));        
                       else
                          error('Unknown parameter');
                       end
                   else                    
                       if isequal(str,'max')
                            a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1))=max(a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1)),b(1:dims(1),(abs(coord(2))+2):dims(2)));                             
                       elseif isequal(str,'min')
                            a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1))=min(a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1)),b(1:dims(1),(abs(coord(2))+2):dims(2)));        
                       else
                          error('Unknown parameter');
                       end              
                   end
               elseif coord(1)<=0 & coord(2)<=0
                   if coord(1)+dims(1)-1<1 | coord(2)+dims(2)-1<1
                       error('Inproper Selection of Starting Pixel');
                   else                   
                       if isequal(str,'max')
                            a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1))=max(a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1)),b((abs(coord(1))+2):dims(1),(abs(coord(2))+2):dims(2)));                             
                       elseif isequal(str,'min')
                            a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1))=min(a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1)),b((abs(coord(1))+2):dims(1),(abs(coord(2))+2):dims(2)));        
                       else
                          error('Unknown parameter');
                       end
                   end
               end
                       
       else                     
                if (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)<=s2                     
                    if isequal(str,'max')
                         a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1))=max(a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1)),b);                             
                    elseif isequal(str,'min')
                         a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1))=min(a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1)),b);        
                    else
                         error('Unknown parameter');
                    end
                elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)<=s2                     
                    if isequal(str,'max')
                         a(coord(1):s1,coord(2):(coord(2)+dims(2)-1))=max(a(coord(1):s1,coord(2):(coord(2)+dims(2)-1)),b(1:(1+s1-coord(1)),1:dims(2)));                             
                    elseif isequal(str,'min')
                         a(coord(1):s1,coord(2):(coord(2)+dims(2)-1))=min(a(coord(1):s1,coord(2):(coord(2)+dims(2)-1)),b(1:(1+s1-coord(1)),1:dims(2)));        
                    else
                         error('Unknown parameter');
                    end
                elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)>s2                      
                    if isequal(str,'max')
                         a(coord(1):(coord(1)+dims(1)-1),coord(2):s2)=max(a(coord(1):(coord(1)+dims(1)-1),coord(2):s2),b(1:dims(1),1:(s2-coord(2)+1)));                             
                    elseif isequal(str,'min')
                         a(coord(1):(coord(1)+dims(1)-1),coord(2):s2)=min(a(coord(1):(coord(1)+dims(1)-1),coord(2):s2),b(1:dims(1),1:(s2-coord(2)+1)));        
                    else
                         error('Unknown parameter');
                    end
                elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)>s2                      
                    if isequal(str,'max')
                         a(coord(1):s1,coord(2):s2)=max(a(coord(1):s1,coord(2):s2),b(1:(s1-coord(1)+1),1:(s2-coord(2)+1)));                             
                    elseif isequal(str,'min')
                         a(coord(1):s1,coord(2):s2)=min(a(coord(1):s1,coord(2):s2),b(1:(s1-coord(1)+1),1:(s2-coord(2)+1)));        
                    else
                         error('Unknown parameter');
                    end
                end                                
        end                 
    else
        
        if  coord(1)>s1 | coord(2)>s2 | coord(3)>s3
           error('Inproper Selection of Starting Pixel ');
           return;
          
       elseif coord(1)<=0 | coord(2)<=0 | coord(3)<=0
           if (coord(1)+dims(1)-1)<1 | (coord(2)+dims(2)-1)<1 | (coord(3)+dims(3)-1)<1
               error('Inproper Selection Starting Pixel');
               return;
           end
       end
      
       if coord(3)<=0 & ((coord(3)+dims(3)-1)<=s3)
           ad=abs(coord(3))+1;
           for i=coord(3):(coord(3)+dims(3)-1)
               if i<=0                 
               else                 
                   if coord(1)<=0 & coord(2)>0     
                        if (coord(2)+dims(2)-1)>s2                            
                            if isequal(str,'max')
                                a(1:(coord(1)+dims(1)-1),coord(2):s2,i)=max(a(1:(coord(1)+dims(1)-1),coord(2):s2,i),b((abs(coord(1))+2):dims(1),1:(s2-coord(2)+1),(i+ad)));                             
                            elseif isequal(str,'min')
                                a(1:(coord(1)+dims(1)-1),coord(2):s2,i)=min(a(1:(coord(1)+dims(1)-1),coord(2):s2,i),b((abs(coord(1))+2):dims(1),1:(s2-coord(2)+1),(i+ad)));        
                            else
                                error('Unknown parameter');
                            end                        
                        else                           
                            if isequal(str,'max')
                                a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i)=max(a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i),b((abs(coord(1))+2):dims(1),1:dims(2),(i+ad)));                             
                            elseif isequal(str,'min')
                                a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i)=min(a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i),b((abs(coord(1))+2):dims(1),1:dims(2),(i+ad)));        
                            else
                                error('Unknown parameter');
                            end                        
                        end                      
                   elseif coord(1)>0 & coord(2)<=0              
                        if (coord(1)+dims(1)-1)>s1                               
                            if isequal(str,'max')
                                a(coord(1):s1,1:coord(2)+dims(2)-1,i)=max(a(coord(1):s1,1:coord(2)+dims(2)-1,i),b(1:(s1-coord(1)+1),(abs(coord(2))+2):dims(2),(i+ad)));                             
                            elseif isequal(str,'min')
                                a(coord(1):s1,1:coord(2)+dims(2)-1,i)=min(a(coord(1):s1,1:coord(2)+dims(2)-1,i),b(1:(s1-coord(1)+1),(abs(coord(2))+2):dims(2),(i+ad)));        
                            else
                                error('Unknown parameter');
                            end                          
                        else                              
                            if isequal(str,'max')
                                a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i)=max(a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i),b(1:dims(1),(abs(coord(2))+2):dims(2),(i+ad)));                             
                            elseif isequal(str,'min')
                                a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i)=min(a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i),b(1:dims(1),(abs(coord(2))+2):dims(2),(i+ad)));        
                            else
                                error('Unknown parameter');
                            end                      
                        end              
                  elseif coord(1)<=0 & coord(2)<=0                         
                            if isequal(str,'max')
                                a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i)=max(a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i),b((abs(coord(1))+2):dims(1),(abs(coord(2))+2):dims(2),(i+ad)));                             
                            elseif isequal(str,'min')
                                a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i)=min(a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i),b((abs(coord(1))+2):dims(1),(abs(coord(2))+2):dims(2),(i+ad)));        
                            else
                                error('Unknown parameter');
                            end
                  elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)<=s2                       
                            if isequal(str,'max')
                                a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i)=max(a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i),b(:,:,(i+ad)));                             
                            elseif isequal(str,'min')
                                a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i)=min(a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i),b(:,:,(i+ad)));        
                            else
                                error('Unknown parameter');
                            end           
                  elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)<=s2                       
                            if isequal(str,'max')
                                a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i)=max(a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i),b(1:(1+s1-coord(1)),1:dims(2),(i+ad)));                             
                            elseif isequal(str,'min')
                                a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i)=min(a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i),b(1:(1+s1-coord(1)),1:dims(2),(i+ad)));        
                            else
                                error('Unknown parameter');
                            end                                       
                  elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)>s2                       
                            if isequal(str,'max')
                                a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i)=max(a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i),b(1:dims(1),1:(s2-coord(2)+1),(i+ad)));                             
                            elseif isequal(str,'min')
                                a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i)=min(a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i),b(1:dims(1),1:(s2-coord(2)+1),(i+ad)));        
                            else
                                error('Unknown parameter');
                            end                                       
                  elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)>s2                        
                            if isequal(str,'max')
                                a(coord(1):s1,coord(2):s2,i)=max(a(coord(1):s1,coord(2):s2,i),b(1:(s1-coord(1)+1),1:(s2-coord(2)+1),(i+ad)));                             
                            elseif isequal(str,'min')
                                a(coord(1):s1,coord(2):s2,i)=min(a(coord(1):s1,coord(2):s2,i),b(1:(s1-coord(1)+1),1:(s2-coord(2)+1),(i+ad)));        
                            else
                                error('Unknown parameter');
                            end                                       
                  end
              end
          end
           
       elseif coord(3)>=1  
           ad=coord(3)-1;
           for i=coord(3):(coord(3)+dims(3)-1)
               if i>s3                   
               else                   
                   if coord(1)<=0 & coord(2)>0               
                         if (coord(2)+dims(2)-1)>s2                             
                            if isequal(str,'max')
                                a(1:(coord(1)+dims(1)-1),coord(2):s2,i)=max(a(1:(coord(1)+dims(1)-1),coord(2):s2,i),b((abs(coord(1))+2):dims(1),1:(s2-coord(2)+1),(i-ad)));                             
                            elseif isequal(str,'min')
                                a(1:(coord(1)+dims(1)-1),coord(2):s2,i)=min(a(1:(coord(1)+dims(1)-1),coord(2):s2,i),b((abs(coord(1))+2):dims(1),1:(s2-coord(2)+1),(i-ad)));        
                            else
                                error('Unknown parameter');
                            end                        
                         else                            
                            if isequal(str,'max')
                                a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i)=max(a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i),b((abs(coord(1))+2):dims(1),1:dims(2),(i-ad)));                             
                            elseif isequal(str,'min')
                                a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i)=min(a(1:(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i),b((abs(coord(1))+2):dims(1),1:dims(2),(i-ad)));        
                            else
                                error('Unknown parameter');
                            end                          
                         end             
                   elseif coord(1)>0 & coord(2)<=0               
                        if (coord(1)+dims(1)-1)>s1                               
                             if isequal(str,'max')
                                a(coord(1):s1,1:coord(2)+dims(2)-1,i)=max(a(coord(1):s1,1:coord(2)+dims(2)-1,i),b(1:(s1-coord(1)+1),(abs(coord(2))+2):dims(2),(i-ad)));                             
                            elseif isequal(str,'min')
                                a(coord(1):s1,1:coord(2)+dims(2)-1,i)=min(a(coord(1):s1,1:coord(2)+dims(2)-1,i),b(1:(s1-coord(1)+1),(abs(coord(2))+2):dims(2),(i-ad)));        
                            else
                                error('Unknown parameter');
                            end                                                  
                        else                              
                            if isequal(str,'max')
                                a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i)=max(a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i),b(1:dims(1),(abs(coord(2))+2):dims(2),(i-ad)));                             
                            elseif isequal(str,'min')
                                a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i)=min(a(coord(1):(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i),b(1:dims(1),(abs(coord(2))+2):dims(2),(i-ad)));        
                            else
                                error('Unknown parameter');
                            end                       
                        end               
                   elseif coord(1)<=0 & coord(2)<=0                       
                        if isequal(str,'max')
                            a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i)=max(a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i),b((abs(coord(1))+2):dims(1),(abs(coord(2))+2):dims(2),(i-ad)));                             
                        elseif isequal(str,'min')
                            a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i)=min(a(1:(coord(1)+dims(1)-1),1:(coord(2)+dims(2)-1),i),b((abs(coord(1))+2):dims(1),(abs(coord(2))+2):dims(2),(i-ad)));        
                        else
                            error('Unknown parameter');
                        end                                  
                   elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)<=s2                         
                       if isequal(str,'max')
                           a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i)=max(a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i),b(:,:,(i-ad)));                             
                       elseif isequal(str,'min')
                           a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i)=min(a(coord(1):(coord(1)+dims(1)-1),coord(2):(coord(2)+dims(2)-1),i),b(:,:,(i-ad)));        
                       else
                           error('Unknown parameter');
                       end             
                   elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)<=s2                           
                       if isequal(str,'max')
                           a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i)=max(a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i),b(1:(1+s1-coord(1)),1:dims(2),(i-ad)));                             
                       elseif isequal(str,'min')
                           a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i)=min(a(coord(1):s1,coord(2):coord(2)+dims(2)-1,i),b(1:(1+s1-coord(1)),1:dims(2),(i-ad)));        
                       else
                           error('Unknown parameter');
                       end                                  
                   elseif (coord(1)+dims(1)-1)<=s1 & (coord(2)+dims(2)-1)>s2                         
                       if isequal(str,'max')
                           a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i)=max(a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i),b(1:dims(1),1:(s2-coord(2)+1),(i-ad)));                             
                       elseif isequal(str,'min')
                           a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i)=min(a(coord(1):(coord(1)+dims(1)-1),coord(2):s2,i),b(1:dims(1),1:(s2-coord(2)+1),(i-ad)));        
                       else
                           error('Unknown parameter');
                       end                            
                   elseif (coord(1)+dims(1)-1)>s1 & (coord(2)+dims(2)-1)>s2                            
                       if isequal(str,'max')
                           a(coord(1):s1,coord(2):s2,i)=max(a(coord(1):s1,coord(2):s2,i),b(1:(s1-coord(1)+1),1:(s2-coord(2)+1),(i-ad)));                             
                       elseif isequal(str,'min')
                           a(coord(1):s1,coord(2):s2,i)=min(a(coord(1):s1,coord(2):s2,i),b(1:(s1-coord(1)+1),1:(s2-coord(2)+1),(i-ad)));        
                       else
                           error('Unknown parameter');
                       end                                
                   end
               end
           end
       end     
    end
      
else
    error('Too Many Input Arguments');
end