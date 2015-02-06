%%
%merge subVolumes back into big volume
function volume = tom_os3_merge(volumeSize,template,subVolumes,subVolSize)


if(isempty(subVolSize) || isequal(subVolSize,[1 1 1]))
    volume = subVolumes{1};
    return;
end;
%%
%volume size
volume = zeros(volumeSize);
dim = tom_os3_fftSelect(volume);
%%
%calculate the center of the template
templateCenterX = floor(size(template,1)/2);
templateCenterY = floor(size(template,2)/2);
if(dim == 3)
    templateCenterZ = floor(size(template,3)/2);
else
    templateCenterZ = 0;
end;

%%
%create return value storage
res = {};
counter = 1;
%%
%generate subvolumes of the image
x=1;
while(x <=size(volume,1))
%chech if x and the needed add on of the templateSize /2 
%is of the image
%if so, set lower and upper values to 1 or volumeSize
%%cut out then
    if(x - templateCenterX < 1)
        lowerX = 1;
    else
        lowerX = x - templateCenterX;
    end;
    
    if(x + subVolSize(1) + templateCenterX > volumeSize(1))
        upperX = volumeSize(1);
    else
        upperX = x + subVolSize(1) + templateCenterX;
    end;
%%  
    y=1;      
    while(y <=size(volume,2))
        if(y - templateCenterY < 1)
            lowerY = 1;
        else
            lowerY = y - templateCenterY;
        end;

        if(y + subVolSize(2) + templateCenterY > volumeSize(2))
            upperY = volumeSize(2);
        else
            upperY = y + subVolSize(2) + templateCenterY;
        end;
%%      
        z=1;  
        while(z <=size(volume,3))
            if(z - templateCenterZ < 1)
                lowerZ = 1;
            else
                lowerZ = z - templateCenterZ;
            end;

            if(z + subVolSize(3) + templateCenterZ > volumeSize(3))
                upperZ = volumeSize(3);
            else
                upperZ = z + subVolSize(3) + templateCenterZ;
            end;
            
            z = upperZ+1;
            
%%            
%          merge back
            
           volume(lowerX:upperX,lowerY:upperY,lowerZ:upperZ) = subVolumes{counter};
           counter = counter +1;
        end;
        
%%        
        y = upperY+1;
    end;
%%    
    x = upperX+1;
    
    
end;


