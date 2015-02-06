%%
function [res pos]= tom_os3_subVolume(pos,img1,img2,mode)
%reads subvolume of img1 of same size as img2 at position pos depending on
%mode - center or edge 
%center means pos is in the center of img2, edge points to the very first
%element of img2
dim = tom_os3_fftSelect(img1);
res = 0;
if(length(pos) == 2)
    pos = [pos 1];
end;

%%
if(pos(1) > size(img1,1) || pos(2) > size(img1,2) || pos(3) > size(img1,3) || pos(1) <= 0 || pos(2) <= 0 || pos(3) <= 0)
    return;
elseif(size(img1) == size(img2))
    res = img1;
    pos = [1 size(img1,1) 1 size(img1,2) 1 size(img1,3)];
    return;
end;    


%%
imgSize = size(img2);
if(strcmp(mode,'edge'))
    lowerX = pos(1) ;
    upperX = pos(1) + imgSize(1);
    
    
    if(lowerX < 1 || upperX > size(img1,1))
        return;
    end;
    
    
    lowerY = pos(2) ;
    upperY = pos(2) + imgSize(2);
    
    if(lowerY < 1 || upperY > size(img1,2))
        return;
    end;
    
    
    lowerZ = 1;
    upperZ = 1;
    if(dim == 3)
        lowerZ = pos(3) ;
        upperZ = pos(3) + imgSize(3);
    end;
    if(lowerZ < 1 || upperZ > size(img1,3))
        return;
    end;    

elseif(strcmp(mode,'center'))

    lowerX = pos(1) - imgSize(1)/2;
    upperX = pos(1) + imgSize(1)/2-1;
    
    if(lowerX < 1) 
        lowerX = 1;
    end;
    if(upperX > size(img1,1))
        upperX = size(img1,1);
    end;
    
    lowerY = pos(2) - imgSize(2)/2;
    upperY = pos(2) + imgSize(2)/2-1;
    
    if(lowerY < 1)
        lowerY = 1;
    end;
    if(upperY > size(img1,2))
        upperY = size(img1,2);
    end;
    
    lowerZ = 1;
    upperZ = 1;
    if(dim == 3)
        lowerZ = pos(3) - imgSize(3)/2;
        upperZ = pos(3) + imgSize(3)/2-1;
    end;
    
    if(lowerZ < 1) 
        lowerZ = 1;
    end;
    if(upperZ > size(img1,3))
       upperZ = size(img1,3);
    end;
    
end;

res = img1(lowerX:upperX,lowerY:upperY,lowerZ:upperZ);
pos = [lowerX upperX lowerY upperY lowerZ upperZ];
