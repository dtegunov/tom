function [meanVol iMean fftImg fftMask] = tom_os3_mean(img,mask,maskSize,fftImg,fftMask)
%tom_os3_mean
%   
% 
% 
%   tom_os3_mean(img,mask,maskSize,fftImg,fftMask)
%
%PARAMETERS
%
%  INPUT
%   img         - the search image / volume
%   mask        - the mask used for normalization
%   fftImg      - optional. The image already transformed to fourierspace
%   fftMask     - optional. The mask already transformed to fourierspace
%  
%  OUTPUT
%   meanVol     - the resulting mean volume
%   iMean       - the center value of meanVol
%   fftImg      - the image already transformed to fourierspace
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
        mask = tom_os3_pasteCenter(zeros(size(img),'single'),mask);
    end;    
    if(nargin<3)
        maskSize    = sum(sum(sum(mask ~= 0)));
    end;
    
    if(nargin < 4)
        fftImg = fftn(img);
    end;
    if(nargin < 5)
        fftMask = fftn(mask);
    end;
    meanVol = real(ifftshift(ifftn(fftImg.*fftMask)))/maskSize;
    iMean = meanVol(center(1),center(2),center(3));
    
    
%%  checked and working