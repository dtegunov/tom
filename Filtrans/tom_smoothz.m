function image = tom_smoothz(image, border, sigma, flag)
%TOM_SMOOTHZ performs smoothing of borders along z only
%
%   image = tom_smoothz(image, border, sigma, flag)
%
%PARAMETERS
%
%  INPUT
%   image               1,2 or 3D data
%   border              width of border
%   sigma               steepness of smoothing function (default: BORDER)
%   flag                 - 'zero' set values to Zero at edges
%  
%  OUTPUT
%   image               ...
%
%   TOM_SMOOTH(IMAGE, BORDER, SIGMA, FLAG) performs smoothing of borders for 1D, 2D, and 3D data.
%   edges are set to mean value of all border pixels, values get damped by
%   a gaussian function F(n) = (Exp(-2*n**2/SIGMA**2)) -
%   Exp(-2)/(1-Exp(-2)).
%
%EXAMPLE
%   rand('state',sum(100*clock));
%   image = rand(64, 64, 64);
%   out = tom_smoothz(image, 10, 8, 'zero');
%   tom_dspcub(out);
%
%   In the example a random image is smoothened, values are set zero at
%   border.
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FF 02/27/03
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


%   check dimensions
[ix_up iy_up iz_up] = size(image);
%   calculate mean (zero if flag is set to zero)
error(nargchk(2,4,nargin));
if (nargin < 3)
    sigma = border;
end;
if (nargin > 3)
    if strmatch(flag,'zero')
        imdev = 0;
    end;
else
    imdev = 1;
end;
%   calculate mean values at edges
if imdev > 0
    if iz_up == 1 & iy_up == 1
    imdev = sum([image(1),image(ix_up)])/2;
    elseif iz_up == 1 & iy_up > 1
        imdev = (sum(image(1:ix_up,1)) + sum(image(1:ix_up,iy_up))) / 2*size(image,1);
    else
        imdev = (sum(sum(image(1:ix_up,1:iy_up,1))) + sum(sum(image(1:ix_up,1:iy_up,iz_up))))  / (2*size(image,1)*size(image,2));
    end;
end;


%   design filter
%   gaussian
xx = 1:border;
fact = (exp(-2*(xx.^2)./(sigma^2))-exp(-2))./(1-exp(-2));
%   if sigma < border : set values to zero for xx > sigma
if sigma < border
    fact(floor(sigma):border)=0;
end;
%   cosine^2 
%   xx = 1:border;
%   fact(:) = 0.5*((cos((border-xx)/border*pi))+1);

for ix = 1:border,
    if iz_up == 1 & iy_up == 1
        image(ix) = imdev + (image(ix)-imdev)*fact(border-ix+1);
        image(ix_up) = imdev + (image(ix_up)-imdev)*fact(border-ix+1);
    elseif iz_up == 1 & iy_up > 1
        image(ix:ix_up,ix)= imdev + (image(ix:ix_up,ix)-imdev)*fact(border-ix+1);
        image(ix:ix_up,iy_up)= imdev + (image(ix:ix_up,iy_up)-imdev)*fact(border-ix+1);
    else
        image(1:ix_up,1:iy_up,ix)= imdev + ( image(1:ix_up,1:iy_up,ix) -imdev) *fact(border-ix+1);
        image(1:ix_up,1:iy_up,iz_up)= imdev + (image(1:ix_up,1:iy_up,iz_up) -imdev)*fact(border-ix+1);
        iz_up = iz_up-1;
    end;
end;