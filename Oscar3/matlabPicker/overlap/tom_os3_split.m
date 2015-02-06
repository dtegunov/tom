%% split a volume
function [res coordinates]= tom_os3_split(volume,template,subVolSize , mode)
% volume - volume to be split
% tempalate - tempalte volume 
% subVolumeSize -  size of subvolume of img
% mode - return list of subVolumes if mode == volumes
%      - return list of coordinates of subVolumes in the large volume if mode == coordinates  
if( isempty(subVolSize) || isequal(subVolSize,[1 1 1]) || isequal(size(volume),size(template)))
    
    if(strcmp(mode,'volume'))
        res = {volume};
    elseif(isequal(size(volume),size(template)))
        res{1}.coordinatesWithExtension   = [1,size(volume,1),1,size(volume,2),1,size(volume,3)];
        res{1}.shiftVector                = [0,0,0];
    end;
    
    return;
end;

dim = tom_os3_fftSelect(volume);
%% volume size
volumeSize = size(volume);
%% calculate the center of the template
templateCenterX = floor(size(template,1)/2)+1;
templateCenterY = floor(size(template,2)/2)+1;
if(dim == 3)
    templateCenterZ = floor(size(template,3)/2)+1;
else
    templateCenterZ = 0;
end;

%% create return value storage
res = {};
counter = 1;
%% generate subvolumes of the image
x=1;
while(x <=size(volume,1))
%check if x and the needed add on of the templateSize /2 +1
%is out off the image
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
            if(z - templateCenterZ < 1 || templateCenterZ == 0)
                %additional constraint is for 
                lowerZ = 1;
            else
                lowerZ = z - templateCenterZ;
            end;

            if(z + subVolSize(3) + templateCenterZ > volumeSize(3) || templateCenterZ == 0)
                upperZ = volumeSize(3);
            else
                upperZ = z + subVolSize(3) + templateCenterZ;
            end;
            
%%	cut out and paste
            
            if(strcmp(mode,'volumes'))
                %return list of subVolumes
                res{counter} = volume(lowerX:upperX,lowerY:upperY,lowerZ:upperZ);
            else
                %return list of coordinates of subvolumes in the large
                %volume
                
                
                lowerX = lowerX - (lowerX>1 && mod(upperX - lowerX-1,2) ~= 0);
                upperX = upperX + (upperX<volumeSize(1) &&  mod(upperX - lowerX-1,2) ~= 0);
                
                lowerY = lowerY - (lowerY>1 && mod(upperY - lowerY-1,2) ~= 0);
                upperY = upperY + (upperY<volumeSize(2) &&  mod(upperY - lowerY-1,2) ~= 0);
                
                lowerZ = lowerZ - (lowerZ>1 && mod(upperZ - lowerZ-1,2) ~= 0);
                upperZ = upperZ + (upperZ<volumeSize(3) &&  mod(upperZ- lowerZ-1,2) ~= 0);
                
                
                res{counter}.coordinatesWithExtension   = [lowerX,upperX,lowerY,upperY,lowerZ,upperZ];
                res{counter}.shiftVector                = [x - lowerX, y - lowerY , z - lowerZ];
                
                res{counter}.subVolumeSize= subVolSize;
                
                if(numel(lowerX:upperX) < subVolSize)
                    res{counter}.subVolumeSize(1) = numel(lowerX:upperX) -(x - lowerX);
                end;
                
                if(numel(lowerY:upperY) < subVolSize)
                    res{counter}.subVolumeSize(2) = numel(lowerY:upperY);
                end;
                
                if(numel(lowerZ:upperZ) < subVolSize)
                    res{counter}.subVolumeSize(3) = numel(lowerZ:upperZ);
                end;
                
            end;
%%  loop            
            counter = counter +1;
            z = z + subVolSize(3);
        end;
        
%%        
        y = y + subVolSize(2);
    end;
%%    
    x = x + subVolSize(1);
    
    
end;