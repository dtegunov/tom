function vol_out=tom_point_mirror(vol,point_coord)
%TOM_POINT_MIRROR Creation of the mirror of the input Image
%
%   c=tom_mirror(varargin)
%
%   C=TOM_POINT_MIRROR(A,'axe')  
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
%       c=tom_mirror(a,'x')
%
%            1   2   3   4   5          21  22  23  24  25
%            6   7   8   9  10          16  17  18  19  20
%       a = 11  12  13  14  15     c =  11  12  13  14  15
%           16  17  18  19  20           6   7   8   9  10
%           21  22  23  24  25           1   2   3   4   5
%
%
%
%   c=tom_mirror(a,'y')
%
%            1   2   3   4   5          5    4   3   2   1
%            6   7   8   9  10          10   9   8   7   6
%       a = 11  12  13  14  15     c =  15  14  13  12  11
%           16  17  18  19  20          20  19  18  17  16
%           21  22  23  24  25          25  24  23  22  21
%
%REFERENCES
%
%SEE ALSO
%   TOM_LIMIT, TOM_PASTE, TOM_MOVE, TOM_RED
%
%   created by AL 08/09/02
%   optimized for speed by AK 10/08/07
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


if (nargin<2)
    point_coord=floor(size(vol)./2)+1;
end;



    
    




