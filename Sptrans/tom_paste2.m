function a=tom_paste(varargin)
%TOM_PASTE2 paste one array in another one
%
%   a=TOM_PASTE2(A,B,coord)
%   a=TOM_PASTE2(A,B,coord,option)
%
%   The array B will be pasted into the array A so that its upper left corner
%   will become pixel (coord(1) coord(2)) in the array A.coord(1), coord(2) 
%   may be outside A. The scaling of the rows/columns is: -n -(n-1) ...... -1 0 1 2 3
%   which means that the coord(1), coord(2) can also take the 0 value. This function
%   is compatible with 3-D arrays also. In that case the paste command is applied
%   in every slice of the 3-D image
%
%PARAMETERS
%
%  INPUT
%   A                   destination image/data/matrix
%   B                   original image/data/matrix
%   coord               [x y] or [x y z] coordinate where the original data
%                       is going to be paste
%   option              which method is going to be used to paste
%
%  OUTPUT
%   a                   image/data/matrix after pasting
%
%OPTION
%   The user can use the extra option:
%   'min'       new values will be pasted only if they are smaller than the original data.
%   'max'       new values will be pasted only if they are bigger than the original data.
%   'mean'      The mean value of the original and destination is outputed
%   'trans'     The elements with value=0 of the original data are treated as
%               transparent.(Will not be pasted onto the destination)
%   'meantrans' The elements with value=0 of the original data are treated as
%               transparent.All other elements are pasted with the option
%               'mean'
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
%   updated by FR 01/09/04
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
else
    a=varargin{1};
    b=varargin{2};
    co=varargin{3};
    [d1 d2 d3]=size(b);  
    [s1 s2 s3]=size(a);
    
    if nargin>=4 ,opt=varargin{4}; else opt='keine'; end;
    if nargin==5 ,par=varargin{5};else par=0;  end;
    
    vx=max(1,co(1)):min(s1,co(1)+d1-1);vy=max(1,co(2)):min(s2,co(2)+d2-1);
    wx=max(-co(1)+2,1):min(s1-co(1)+1,d1);wy=max(-co(2)+2,1):min(s2-co(2)+1,d2); 
    
    if not(isempty(vx) | isempty(vy) )
       if d3==1,a(vx,vy)=mix(a(vx,vy),b(wx,wy),opt,par);
        else
            vz=max(1,co(3)):min(s3,co(3)+d3-1);wz=max(-co(3)+2,1):min(s3-co(3)+1,d3);
            if not(isempty(vz)),a(vx,vy,vz)=mix(a(vx,vy,vz),b(wx,wy,wz),opt,par);
            end;
        end;   
    end;
end;

function c=mix(a,b,opt,par);
switch lower(opt)
    case 'max',c=max(a,b);
    case 'min',c=min(a,b);
    case 'mean',c=(a+b)/2;
    case 'power_mean',
        c=(a+b)/2;
        mask=a>0; b_overlap=b.*mask;
        mask2=abs(mask-1); d=c.*mask+b.*mask2;
        c=d;
    case 'trans',mask=(b~=par);c=mask.*b+(1-mask).*a;
    case 'meantrans',
        mask=(b~=par);
        mask2=(mask .* ((-666)==a));
        c=mask.*(b/2)+mask2.*(b/2)+(1-mask).*a +((1-mask2).*mask.*(a/2));
    otherwise,c=b;
end;

        
    