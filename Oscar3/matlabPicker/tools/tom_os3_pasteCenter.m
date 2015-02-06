function res = tom_os3_pasteCenter(img,img2)
%pastes img2 into center  of img  

%% tested && works
dim = tom_os3_fftSelect(img);
if(dim == 2)
    if(size(img,1) == 1)
        img = img';
        img2 = img2';
        dim = 1;
    elseif(size(img,2) == 1)
        dim = 1;
    end;
end;


imsize = size(img);
res = img;


if(dim == 1)
    centerX = floor(imsize(1)/2) - floor(size(img2,1)/2) +1 ;
    
    res(centerX:centerX+size(img2,1)-1) = img2;
elseif(dim == 2)

    
    centerX = floor(imsize(1)/2) - floor(size(img2,1)/2)  +1;
    centerY = floor(imsize(2)/2) - floor(size(img2,2)/2)  +1;
    
    res(centerX:centerX+size(img2,1)-1,centerY:centerY+size(img2,2)-1) = img2;
    
else
    
    
    centerX = floor(imsize(1)/2)  - floor(size(img2,1)/2) +1;
    centerY = floor(imsize(2)/2)  - floor(size(img2,2)/2) +1;
    centerZ = floor(imsize(3)/2)  - floor(size(img2,3)/2) +1;
    res(centerX:centerX+size(img2,1)-1,centerY:centerY+size(img2,2)-1,centerZ:centerZ+size(img2,3)-1) = img2;
    
end;