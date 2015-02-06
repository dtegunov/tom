function a=tom_limit(varargin)
%TOM_LIMIT Setting limits to the values of an image
%
%   a=tom_limit(varargin)
%
%   A=TOM_LIMIT(A,LOW,UP,PROPERTIES) In case of the option PROPERTIES
%   is ommited, the values of the image A bigger than UP are set to UP,
%   and the values smaller than LOW are set to LOW. PROPERTIES is a cell 
%   array containing one string. The only acceptable string are 'z' or 'az'.
%   In case of 'z', the values smaller than LOW are set to zero and the
%   values bigger than UP are set to UP. In case of 'az', the values outside
%   of the limits LOW and UP are set to zero 
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
%            0  -1   1   5
%           -1   2   1   3
%       a =  1  -3   7   4
%           -3   8   1   1
%           -5   2   0   6
%
%   c = tom_limit(a,-1,3,'z')
%            0  -1   1   3
%           -1   2   1   3
%       c =  1   0   3   3  
%            0   3   1   1
%            0   2   0   3
%
%   c = tom_limit(a,-1,3,'az')
%            0  -1   1   0
%           -1   2   1   3
%       c =  1   0   0   0  
%            0   0   1   1
%            0   2   0   0
%
%REFERENCES
%
%SEE ALSO
%   TOM_PEAK, TOM_PASTE, TOM_CIRCLE, TOM_SPHERE, TOM_MIRROR
%
%   created by AL 08/10/02
%   updated by FF 09/24/02
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

switch (nargin)
    
case 0
    error('Not enough Arguments');
    
case 1
    error('Not enough Arguments');
    
case 2 
    error('Not enough Arguments');
    
case 3
    a=varargin{1};
    low=varargin{2};
    up=varargin{3};
    [s1 s2 s3]=size(a);
    ilow = find(a<low);
    a(ilow) = low;
    iup = find(a>up);
    a(iup) = up;
    
case 4
    a=varargin{1};
    low=varargin{2};
    up=varargin{3};
    [s1 s2 s3]=size(a);
    string=varargin{4};
    if string == 'z'
        ilow = find(a<low);
        a(ilow) = 0;
        iup = find(a>up);
        a(iup) = up;
    elseif string == 'az'
        ilow = find(a<low);
        a(ilow) = 0;
        iup = find(a>up);
        a(iup) = 0;
    else
        error('Unknown Argument');
    end
    
case 5
    error('Too many Arguments');
end
