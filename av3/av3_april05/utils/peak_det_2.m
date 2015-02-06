function [varargout]=peak_det_2(varargin)

%PEAK_DET_2 Determination of the maximum value
%   [COORD VAL]=TOM_PEAK(A) Value and position of maximum will  be  determined  by          
%   this call and the result will be reported.
%  
%
%
%   Example :
%  -----------
%
%               1   2   3   4   17               
%               5   6   7   8   18
%           a = 10  9   40  6   20
%               12  16  17  20  30
%               5   23  35  12  6
%
%           [c v m]=tom_peak(a)
%
%           c = [3.0890    2.9770]
%           v = 40.2232
%
%
%
%
% see also : TOM_LIMIT, TOM_PASTE
%       
%  05/07/03  


 x_1 = 1:0.001:3;
 x   = 1:3;

switch (nargin)
  case  0
      error('Not enough Arguments');
      
   %------------------------------------------------------------------------------
  
   
   
  case  1                                       %no radius given    
      
      a=varargin{1};
      [dim1 dim2 dim3]=size(a);
            
      if dim3 == 1                              %2d Matrix
          b=max(max(a));
          [s1 s2]=find(a==b);
          if s1 == 1 | s1==dim1
              c(1)=s1;
              val(1)=b;
          else
              peak(1:3)=a(s1-1:s1+1,s2);
              y_1=interp1( x ,peak, x_1,'spline');
              [val(1),pos_peak]=max(y_1);
              c(1)=s1-1+(pos_peak-1)*0.001;
          end           
          if s2 == 1 | s2==dim2
              c(2)=s2;
              val(2)=b;
          else
              peak(1:3)=a(s1,s2-1:s2+1);
              y_1=interp1( x ,peak, x_1,'spline');
              [val(2),pos_peak]=max(y_1);
              c(2)=s2-1+(pos_peak-1)*0.001;
              val=max(val);
          end

          varargout{1}=c;
          varargout{2}=val;
          
          
    else                                                          % 3d Matrix
        b=max(max(max(a)));
        [s1 s2 s3]=find(a==b);         
        if floor(s2/dim2)==(s2/dim2)
            s3=s2/dim2;
            ss2=dim2;
        else
            s3 = floor(s2/dim2) + 1;
            ss2=rem(s2,dim2);
        end
        c=[s1 ss2 s3];
          
        if s1 == 1 | s1 == dim1
            c(1)=s1;
            val(1)=b;
        else
            peak(1:3)=a(s1-1:s1+1,ss2,s3);
            y_1=interp1( x ,peak, x_1,'spline');
            [val(1),pos_peak]=max(y_1);
            c(1)=s1-1+(pos_peak-1)*0.001;
        end
        if ss2 == 1 | ss2 == dim2
            c(2)=ss2;
            val(2)=b;
        else 
            peak(1:3)=a(s1,ss2-1:ss2+1,s3);
            y_1=interp1( x ,peak, x_1,'spline');
            [val(2),pos_peak]=max(y_1);
            c(2)=ss2-1+(pos_peak-1)*0.001;
        end
        
        if s3 == 1 | s3==dim3
            c(3)=s3;
            val(3)=b;
        else 
            peak(1:3)=a(s1,ss2,s3-1:s3+1);
            y_1=interp1( x ,peak, x_1,'spline');
            [val(3),pos_peak]=max(y_1);
            c(3)=s3-1+(pos_peak-1)*0.001;            
        end
        
        val=max(val);
        varargout{1}=c;
        varargout{2}=val;
       
        
    end
  
 

%--------------------------------------------------------------------------


  case  2                   % Two Input Arguments, the Matrix and the desired Radius
      a=varargin{1};
      radius=varargin{2};
      [dim1 dim2 dim3]=size(a);
         
      if dim3 == 1                                %2d Matrix
          b=max(max(a));
          [s1 s2]=find(a==b);
          if s1 == 1 | s1==dim1
              c(1)=s1;
              val(1)=b;
          else
              peak(1:3)=a(s1-1:s1+1,s2);
              y_1=interp1( x ,peak, x_1,'spline');
              [val(1),pos_peak]=max(y_1);
              c(1)=s1-1+(pos_peak-1)*0.001;
          end     
          
          if s2 == 1 | s2==dim2
              c(2)=s2;
              val(2)=b;
          else
              peak(1:3)=a(s1,s2-1:s2+1);
              y_1=interp1( x ,peak, x_1,'spline');
              [val(2),pos_peak]=max(y_1);
              c(2)=s2-1+(pos_peak-1)*0.001;
              val=max(val);
          end
          b=max(val);
          
          
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
        circle=tom_circle2(radius);
        circle=1-circle;
        center=radius+1;
        tmp2=circle((center-xradl):(center+xradu),(center-yradl):(center+yradu));
        tmp=tmp1.*tmp2;
        a((s1-xradl):(s1+xradu),(s2-yradl):(s2+yradu))=tmp;
        
        varargout{1}=c;
        varargout{2}=b;
        varargout{3}=a,
        
        
    else                  % 3d Matrix
        
        b=max(max(max(a)));
        [s1 s2 s3]=find(a==b);         
        if floor(s2/dim2)==(s2/dim2)
            s3=s2/dim2;
            ss2=dim2;
        else
            s3 = floor(s2/dim2) + 1;
            ss2=rem(s2,dim2);
        end
        c=[s1 ss2 s3];
          
        if s1 == 1 | s1 == dim1
            c(1)=s1;
            val(1)=b;
        else
            peak(1:3)=a(s1-1:s1+1,ss2,s3);
            y_1=interp1( x ,peak, x_1,'spline');
            [val(1),pos_peak]=max(y_1);
            c(1)=s1-1+(pos_peak-1)*0.001;
        end
        if ss2 == 1 | ss2 == dim2
            c(2)=ss2;
            val(2)=b;
        else 
            peak(1:3)=a(s1,ss2-1:ss2+1,s3);
            y_1=interp1( x ,peak, x_1,'spline');
            [val(2),pos_peak]=max(y_1);
            c(2)=ss2-1+(pos_peak-1)*0.001;
        end
        
        if s3 == 1 | s3==dim3
            c(3)=s3;
            val(3)=b;
        else 
            peak(1:3)=a(s1,ss2,s3-1:s3+1);
            y_1=interp1( x ,peak, x_1,'spline');
            [val(3),pos_peak]=max(y_1);
            c(3)=s3-1+(pos_peak-1)*0.001;            
        end
        
        
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
        
        s2=ss2;
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
        
        
        
        sphere=tom_sphere(radius);
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
        varargout{2}=b;
        varargout{3}=a;
               
        
    end
        
        
end


