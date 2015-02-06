function [varargout] = tom_peak(varargin)
% TOM_PEAK determines the coordinates and value of the maximum of an array
% (2d, 3d).
% [c val] = tom_peak(A)
%
% [c val] = tom_peak(A,'spline')
% If the flag 'spline' is set the coordinates are determined with subpixel
% accuracy. When there are several pixels with the same maximum value, only the first will be reported.(A spline function is fitted between the next neighbours of the pixel with max value. 
% The maxima of the interpolation functions are determined for each dimension independently.) 
%
% [c val M] = tom_peak(A,R)
% After peak determination, all array elements within a circle (sphere for 3d matrix) of radius R centered
% at the maximum will be set to zero. R must be double.
% 
% PARAMETERS
%  INPUT
%   A       array - 1d, 2d or 3d
%   R       double number; Radius to set values to 0 
%   FLAG    can be set to 'spline' for subpixel accuracy
%
%  OUTPUT
%   c       coordinates of the maximum 
%   val     value of the maximum
%   M       array with all elements set to zero in the circle R
%
%
% EXAMPLE
%
%               1   2   3   4   17               
%               5   6   7   8   18
%           A = 10  9   40  6   20
%               12  16  17  20  30
%               5   23  35  12  6
%
%           [c val] = tom_peak(A)
%           c = [3 3]
%           val = 40
%
%           [c val] = tom_peak(A,'spline')
%           c = [3.0890    2.9770]
%           val = 40.2232
%           
%           [c val M] = tom_peak(A,2,'spline')
%           c = [3.0890    2.9770]
%           val = 40.2232  
%               1   2   0   4   17
%               5   0   0   0   18
%           M = 0   0   0   0   0
%               12  0   0   0   30
%               5   23  0   12  6
%
%REFERENCES
%
% SEE ALSO
%    TOM_LIMIT, TOM_PASTE
%   created by AL 11/09/02
%   updated by FF 10/11/04
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

 x_1 = 1:0.001:3;
 x   = 1:3;

switch (nargin)
  case  0
      error('Not enough Arguments');
      
   %------------------------------------------------------------------------------ 
   
   
  case  1                                    % no radius given, no interpolation
      
                                             % output arguments: c    position
                                             %                   v    value of maximum
                                                
      a=varargin{1};
      
      [dim1 dim2 dim3]=size(a);
      if dim3 == 1  
          v=max(max(a));
          [s1 s2]=find(a==v);
          c=[s1 s2];
          varargout{1}=c;
          varargout{2}=v;
      else
          v=max(max(max(a)));
          [s1 s2 s3]=find(a==v);
          %if floor(s2/dim2)==(s2/dim2)
          %    s3=s2/dim2;
          %    ss2=dim2;
          %else
              %s3 = floor((s2)/dim2) + 1;
              s3 = floor((s2-1)/dim2) + 1;
              ss2=rem(s2-1,dim2)+1;
          %end
          c=[s1 ss2 s3];
          varargout{1}=c;
          varargout{2}=v;
      end
      
      %--------------------------------------------------------------------
      
      
  case  2     
      a=varargin{1};
      [dim1 dim2 dim3]=size(a) ;
      
      if strcmp(varargin{2},'spline')               % no radius given, spline interpolation used
                                                    % for subpixel accuracy
                                                    % output arguments: c    position
                                                    %                   v    value of maximum val
                                                   
          [c,v] = tom_intpeak(a); 
          varargout{1}=c;
          varargout{2}=v;
        
   
              
              
         %-----------------------------------------------------------------
      
         
       elseif class(varargin{2})== 'double'               % no interpolation 
                                                          % radius given
                                                          % output values:  c   position
                                                          %                 v   value of maximum
                                                          %                 a   new matrix with zero-circle
          
          radius=varargin{2};       
              if dim3 == 1                                % 2d Matrix
                  v=max(max(a));
                  [s1 s2]=find(a==v);
                  s1=s1(1);s2=s2(1);
                  c=[s1 s2];
                  
                  if (s1-radius)<=0
                      chk=s1-radius;
                      xradl=radius; % xradl : x radius lower border
                      while chk<=0
                          chk=chk+1;
                          xradl=xradl-1;   
                          
                      end
                  else
                      xradl=radius;
                  end
                  if (s1+radius)>dim1
                      chk=s1+radius;
                      xradu=radius;  % xradu : x radius upper border
                      while chk>dim1
                          chk=chk-1;
                          xradu=xradu-1;   
                          
                      end
                  else
                      xradu=radius;
                  end
                  
                  if (s2-radius)<=0
                      chk=s2-radius;
                      yradl=radius; % yradl : y radius lower border
                      while chk<=0
                          chk=chk+1;
                          yradl=yradl-1;
                      end
                  else
                      yradl=radius;
                  end
                  if (s2+radius)>dim2
                      chk=s2+radius;
                      yradu=radius; % yradu : y radius upper border
                      while chk>dim2
                          chk=chk-1;
                          yradu=yradu-1;
                      end
                  else
                      yradu=radius;
                  end
                  
                  %----Setting the circular or spherical neighborhood to zero
                  
                  tmp1=a((s1-xradl):(s1+xradu),(s2-yradl):(s2+yradu));
                  circle=tom_circle(radius);%circle=tom_circle2(radius);changed 19/07 WDN
                  circle=1-circle;
                  center=radius+1;
                  tmp2=circle((center-xradl):(center+xradu),(center-yradl):(center+yradu));
                  tmp=tmp1.*tmp2;
                  a((s1-xradl):(s1+xradu),(s2-yradl):(s2+yradu))=tmp;
                     
                  varargout{1}=c;
                  varargout{2}=v;
                  varargout{3}=a;
                      
              else                                                   % 3d Matrix                  
                  v=max(max(max(a)));
                  [s1 s2 s3]=find(a==v);  
                  s1=s1(1);s2=s2(1);s3=s3(1);
                  s3 = floor((s2-1)/dim2) + 1;
                  s2=rem(s2-1,dim2)+1;
%                   if floor(s2/dim2)==(s2/dim2)
%                       s3=s2/dim2;
%                       ss2=dim2;
%                   else
%                       s3 = floor(s2/dim2) + 1;
%                       ss2=rem(s2,dim2);
%                   end
                  c=[s1 s2 s3];
                  
                  if (s1-radius)<=0
                      chk=s1-radius;
                      xradl=radius; % xradl : x radius lower border
                      while chk<=0
                          chk=chk+1;
                          xradl=xradl-1;   
                          
                      end
                  else
                      xradl=radius;
                  end
                  if (s1+radius)>dim1
                      chk=s1+radius;
                      xradu=radius;  % xradu : x radius upper border
                      while chk>dim1
                          chk=chk-1;
                          xradu=xradu-1;                           
                      end
                  else
                      xradu=radius;
                  end
                  
                  if (s2-radius)<=0
                      chk=s2-radius;
                      yradl=radius; % yradl : y radius lower border
                      while chk<=0
                          chk=chk+1;
                          yradl=yradl-1;
                      end
                  else
                      yradl=radius;
                  end
                  if (s2+radius)>dim2
                      chk=s2+radius;
                      yradu=radius; % yradu : y radius upper border
                      while chk>dim2
                          chk=chk-1;
                          yradu=yradu-1;
                      end
                  else
                      yradu=radius;
                  end
                  
                  
                  
                  sphere=tom_sphere([2*radius+1 2*radius+1 2*radius+1],radius);
                  sphere=1-sphere;
                  center=radius+1;
                  for j=(-radius):radius
                      if (s3+j)<=0 | (s3+j)>dim3
                          
                      else                
                          tmp1=a((s1-xradl):(s1+xradu),(s2-yradl):(s2+yradu),s3+j);
                          tmp2=sphere((center-xradl):(center+xradu),(center-yradl):(center+yradu),((radius+1)+j));
                          tmp=tmp1.*tmp2;
                          a((s1-xradl):(s1+xradu),(s2-yradl):(s2+yradu),s3+j)=tmp;
                      end
                  end
                  varargout{1}=c;
                  varargout{2}=v;
                  varargout{3}=a;
              end
          else
              error('r must be double');
          end
          
         
          %----------------------------------------------------------------
          
      case 3                                              % interpolation used to determine max value and position
                                                          % radius given
                                                          % output values:  c   position
                                                          %                 v   value of maximum
                                                          %                 a   new matrix with zero-circle
          
          a=varargin{1};
          radius=varargin{2}; 
          
          [c,val]=tom_intpeak(a); 
          [dim1 dim2 dim3]=size(a);
          
        
          if dim3 == 1                                %2d Matrix
              v=max(max(a));
              [s1 s2]=find(a==v);
              s1=s1(1);s2=s2(1);
             
                  %----- Checking if the Radius exceeds Matrix's Dimensions
                  
                  if (s1-radius)<=0
                      chk=s1-radius;
                      xradl=radius; % xradl : x radius lower border
                      while chk<=0
                          chk=chk+1;
                          xradl=xradl-1;   
                          
                      end
                  else
                      xradl=radius;
                  end
                  if (s1+radius)>dim1
                      chk=s1+radius;
                      xradu=radius;  % xradu : x radius upper border
                      while chk>dim1
                          chk=chk-1;
                          xradu=xradu-1;   
                          
                      end
                  else
                      xradu=radius;
                  end
                  
                  if (s2-radius)<=0
                      chk=s2-radius;
                      yradl=radius; % yradl : y radius lower border
                      while chk<=0
                          chk=chk+1;
                          yradl=yradl-1;
                      end
                  else
                      yradl=radius;
                  end
                  if (s2+radius)>dim2
                      chk=s2+radius;
                      yradu=radius; % yradu : y radius upper border
                      while chk>dim2
                          chk=chk-1;
                          yradu=yradu-1;
                      end
                  else
                      yradu=radius;
                  end
                  
                  %----Setting the circular or spherical neighborhood to zero
                  
                  tmp1=a((s1-xradl):(s1+xradu),(s2-yradl):(s2+yradu));
                  circle=tom_circle(radius);%circle=tom_circle2(radius);changed 19/07 WDN
                  circle=1-circle;
                  center=radius+1;
                  tmp2=circle((center-xradl):(center+xradu),(center-yradl):(center+yradu));
                  tmp=tmp1.*tmp2;
                  a((s1-xradl):(s1+xradu),(s2-yradl):(s2+yradu))=tmp;
                  varargout{1}=c;
                  varargout{2}=v;
                  varargout{3}=a;
              else                                                   % 3d Matrix                  
                  v = max(max(max(a)));
                  [s1 s2 s3] = find(a==v); 
                  s1 = s1(1);s2=s2(1);s3=s3(1);
                  s3 = floor((s2-1)/dim2) + 1;
                  s2=rem((s2-1),dim2)+1;
%                   if floor(s2/dim2)==(s2/dim2)
%                       s3=s2/dim2;
%                       ss2=dim2;
%                   else
%                       s3 = floor(s2/dim2) + 1;
%                       ss2=rem(s2,dim2);
%                   end
                                   
                  if (s1-radius)<=0
                      chk=s1-radius;
                      xradl=radius; % xradl : x radius lower border
                      while chk<=0
                          chk=chk+1;
                          xradl=xradl-1;   
                      end
                  else
                      xradl=radius;
                  end
                  if (s1+radius)>dim1
                      chk=s1+radius;
                      xradu=radius;  % xradu : x radius upper border
                      while chk>dim1
                          chk=chk-1;
                          xradu=xradu-1;  
                      end
                  else
                      xradu=radius;
                  end
                  if (s2-radius)<=0
                      chk=s2-radius;
                      yradl=radius; % yradl : y radius lower border
                      while chk<=0
                          chk=chk+1;
                          yradl=yradl-1;
                      end
                  else
                      yradl=radius;
                  end
                  if (s2+radius)>dim2
                      chk=s2+radius;
                      yradu=radius; % yradu : y radius upper border
                      while chk>dim2
                          chk=chk-1;
                          yradu=yradu-1;
                      end
                  else
                      yradu=radius;
                  end
                  sphere = tom_sphere([2*radius+1 2*radius+1 2*radius+1],radius);
                  sphere=1-sphere;
                  center=radius+1;
                  for j=(-radius):radius
                      if (s3+j)<=0 | (s3+j)>dim3
                          
                      else                
                          tmp1=a((s1-xradl):(s1+xradu),(s2-yradl):(s2+yradu),s3+j);
                          tmp2=sphere((center-xradl):(center+xradu),(center-yradl):(center+yradu),((radius+1)+j));
                          tmp=tmp1.*tmp2;
                          a((s1-xradl):(s1+xradu),(s2-yradl):(s2+yradu),s3+j)=tmp;
                      end
                  end
                  varargout{1}=c;
                  varargout{2}=v;
                  varargout{3}=a;
              end
end
 

function [c,val]=tom_intpeak(a)

% TOM_INTPEAK determines Value and position of maximum of a matrix to
%   subpixel accuracy(a spline function is fitted between the next neighbours of the pixel with maximum value. 
%   The maximums of the interpolation functions are determined for each dimension independently.)  
%  
%   07.04.04 GS


 [dim1 dim2 dim3]=size(a);
 
 x_1 = 1:0.001:3;
 x   = 1:3;

 
 if dim3 == 1                              %2d Matrix
     [v, s1] = max(a(:));
     [s1 s2] = ind2sub(size(a), s1);    % changed: Thomas Haller, 03.12.2008
     %v=max(max(a));
     %[s1 s2]=find(a==v);
     %s1=s1(1); s2=s2(1);     
     if s1 == 1 | s1==dim1
         c(1)=s1;
         val(1)=v;
     else
         peak(1:3)=a(s1-1:s1+1,s2);
         y_1=interp1( x ,peak, x_1,'spline');
         [val(1),pos_peak]=max(y_1);
         c(1)=s1-1+(pos_peak-1)*0.001;
     end           
     if s2 == 1 | s2==dim2
         c(2)=s2;
         val(2)=v;
     else
         peak(1:3)=a(s1,s2-1:s2+1);
         y_1=interp1( x ,peak, x_1,'spline');
         [val(2),pos_peak]=max(y_1);
         c(2)=s2-1+(pos_peak-1)*0.001;
     end
     
     val=max(val);
     
 else                                                          % 3d Matrix
     v=max(max(max(a)));
     [s1 s2 s3] = find(a==v);    
     s1=s1(1); s2=s2(1); s3=s3(1);
     s3=floor((s2-1)/dim2) + 1;
     s2=rem((s2-1),dim2)+1;
%      
%      if floor(s2/dim2)==(s2/dim2)
%          s3=s2/dim2;
%          ss2=dim2;
%      else
%          s3 = floor(s2/dim2) + 1;
%          ss2=rem(s2,dim2);
%      end
     
     if s1 == 1 | s1 == dim1
         c(1)=s1;
         val(1)=v;
     else
         peak(1:3)=a(s1-1:s1+1,s2,s3);
         y_1=interp1( x ,peak, x_1,'spline');
         [val(1),pos_peak]=max(y_1);
         c(1)=s1-1+(pos_peak-1)*0.001;
     end
     if s2 == 1 | s2 == dim2
         c(2) = s2;
         val(2)=v;
     else 
         peak(1:3)=a(s1,s2-1:s2+1,s3);
         y_1=interp1( x ,peak, x_1,'spline');
         [val(2),pos_peak]=max(y_1);
         c(2)=s2-1+(pos_peak-1)*0.001;
     end
     
     if s3 == 1 | s3==dim3
         c(3)=s3;
         val(3)=v;
     else 
         peak(1:3)=a(s1,s2,s3-1:s3+1);
         y_1=interp1( x ,peak, x_1,'spline');
         [val(3),pos_peak]=max(y_1);
         c(3)=s3-1+(pos_peak-1)*0.001;    
         
     end
     
     val=max(val);
 end
           
