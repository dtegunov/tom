function c=tom_error(varargin)
%TOM_ERROR creates ...
%
%   c=tom_error(varargin)
%
%PARAMETERS
%
%  INPUT
%   C=TOM_ERROR(A,'G',M,V) This function adds Gaussian noise to the input image A. To each  
%   element of the array A a random number is added. The random numbers follow a Gaussian 
%   distribution with mean value M and variance V. Deafault values for M and V are 0 and
%   0,01.
%   C=TOM_ERROR(A,'P') This function adds Poisson Noise to the input image A. To each
%   element of the array A is written an integer random number. The random number follows
%   a Poisson distribution. After the generation of the noise, the noise is
%   added to the image.
%  
%  OUTPUT
%   c                   ...
%
%EXAMPLE
%     c=tom_error(a,'G',0,1)
%
%            225  225  225  226  226          212  186  189  249  221
%            225  225  225  226  226          182  255  218  236  249
%       a =  226  226  225  226  226   c =    255  242  235  182  232
%            226  226  225  226  226          215  238  230  255  192
%            225  226  225  225  226          155  229  203  209  255
%
%REFERENCES
%
%SEE ALSO
%   TOM_PEAK, TOM_LIMIT, TOM_FILTER, TOM_PASTE
%
%   created by AL 08/09/02
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

switch(nargin)
    
case 0
    error('Not Enough Arguments');
    
case 1
    error('Not Enough Arguments');
    
case 2    
    a=varargin{1};
    str=varargin{2};
    [s1 s2 s3]=size(a);
    if str == 'G'
        error('Not enough Arguments for Gaussian noise');
    elseif str == 'P'
        for i=1:s3
            c=imnoise(a(:,:,i),'poisson');
            a(:,:,i)=a(:,:,i)+c;
        end
    else
        error('Unknown Parameter');
    end
case 3
    error('Not enough Arguments');
case 4
    a=varargin{1};
    str=varargin{2};
    [s1 s2 s3]=size(a);
    m=varargin{3};
    v=varargin{4};
    if str == 'G'
        %for i=1:s3            
            c=a+[sqrt(v)*randn(s1,s2,s3)+m];
        %end removed funny loop - YC and FF
    elseif str == 'P'
        ('Too many Arguments for Poisson noise');
    else
        error('Unknown Parameter');
    end
    
end

