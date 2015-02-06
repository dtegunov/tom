function [imageSTD iSTD fftMask] = tom_os3_std(img,imageMean,mask,maskSize,fftImage,fftMask)
%tom_os3_std
%   
% 
%   tom_os3_std(img,imageMean,mask,maskSize,fftImg,fftMask)
%   
%
%PARAMETERS
%
%  INPUT
%   img         - the search image / volume
%   imageMean   - the mean volume 
%   mask        - the mask used for normalization
%   fftImage    - the image in fourierspace
%   fftMask     - optional, the mask in fourierspace
%
%  OUTPUT
%   imageSTD    - the resulting mean volume
%   iSTD        - the center value of meanVol
%   fftMask     - optional. The mask already transformed to fourierspace
%
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by TH2 07/07/07
%   updated by ..
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
    center = floor(size(img)/2)+1;
    if(length(center) == 2)
        center(3) = 1;
    end;

    if(~isequal(size(img),size(mask)))
        mask = tom_os3_pasteCenter(zeros(size(img)),mask);
    end;    
    
    if(nargin<4)
        maskSize    = sum(sum(sum(mask ~= 0)));
    end;
    
    if(nargin < 5)
        fftImage = fftn(img);
    end;
    
    if(nargin < 6)
        fftMask = fftn(mask);
    end;
    
    imageSTD = real(ifftshift(ifftn(fftn(img.*img).*fftMask)));
    imageVAR = (imageSTD / (maskSize) - imageMean.*imageMean);
    imageVAR = max(imageVAR,0);
    imageSTD = sqrt(imageVAR);
    iSTD     = imageSTD(center(1),center(2),center(3));
    
%%  checked and working