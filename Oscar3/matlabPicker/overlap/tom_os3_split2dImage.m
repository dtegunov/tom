function subImages = tom_os3_split2dImage(image,template,subVolumeSize,mode)
%tom_os3_split2dImage
%   
%   tom_os3_split2dImage(image,template,subVolumeSize)
%
%PARAMETERS
%
%  INPUT
%   image       - the search image 
%   template    - the template searched
%   subVolumeSize- number of subImages per image axis (x and y)
%   
%  
%  OUTPUT
%   subImages   - stack of subImages
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

if(subVolumeSize == 1)
    subImages = {image};
    return;
end;

imageSize = size(image);
templateSize = size(template);

subImages = {};
offset = templateSize(1)/2;
for x = 0:subVolumeSize-1
    for y = 0:subVolumeSize-1
        X1 = 1 + x*imageSize(1)/subVolumeSize;
        if(X1 >1)
            X1 = X1 - templateSize(1)/2-offset;
        end;
        X2 = x*imageSize(1)/subVolumeSize + imageSize(1)/subVolumeSize;
        if(X2 < imageSize(1))
            X2 = X2 + templateSize(1)/2+offset;
        end;
        
        
        Y1 = 1 + y*imageSize(2)/subVolumeSize;
        if(Y1 >1)
            Y1 = Y1 - templateSize(2)/2-offset;
        end;
        Y2 = y*imageSize(2)/subVolumeSize + imageSize(2)/subVolumeSize;
        if(Y2 < imageSize(2))
            Y2 =Y2 + templateSize(2)/2+offset;
        end;        
        if(strcmp(mode,'images'))
            subImages{length(subImages)+1} = image(X1:X2,Y1:Y2);
        else
            subImages{length(subImages)+1} = [X1 X2 Y1 Y2];
        end;
    end;
end;
