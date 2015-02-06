function [noiseVolume noisyVolume]= tom_os3_generateNoise(volume,noiseType,SNR,mask,volumeMean,volumeSTD)
% volume = the volume
% noisetype = additive | multiplicative(not jet)
% SNR iterator
% volumeMean - mean of the real signal (either all values ~= 0 or all values under a mask)
% volumeSTD - std of the real signal (same)

volumeSize = size(volume);
% volumeMean = mean(volume(find(volume ~= 0)));

normalisation = true;
if(strcmp(noiseType,'additive'))
    
    
    if(SNR == 0)
        error('Using SNR 0 is not possible');
    end;
    
    if(nargin <= 4)
        volumeMean= mean(volume(:));
        volumeSTD = std(volume(:));
        mask=ones(size(volume));
        normalisation = false;
    end;

    if(sum(size(SNR))>2)
        vs = volumeSize/2;
        counter =1;
        coor = zeros(4);
%         for x=1:2
%             for y=1:2
% %                 coor(counter,1) = x+(x-1)*vs(1);
% %                 coor(counter,2) = y+(y-1)*vs(2);
% %                 coor(counter,3) = x+(x-1)*vs(1)+vs(1)-1;
% %                 coor(counter,4) = y+(y-1)*vs(2)+vs(2)-1;
%                 
%                 counter = counter + 1;
%             end;
%         end;    
%         
        
        coor(:,1) = [1;1;512;512];
        coor(:,2) = [1;513;512;1024];
        coor(:,3) = [513;1;1024;512];
        coor(:,4) = [513;513;1024;1024];

        
        noisyVolume = zeros(volumeSize);
        
        counter = 1;
        noisyVolume = zeros(volumeSize);
        for i=SNR(1):SNR(2):SNR(3)
            varNoise = volumeSTD^2 / i;
            
            noiseVolume= random('norm',volumeMean,sqrt(varNoise),vs);
            m = find(mask == 0);
            mask(m) = 1;
            if(normalisation)
                noisyVolume(coor(1,counter):coor(3,counter),coor(2,counter):coor(4,counter))= (volume(coor(1,counter):coor(3,counter),coor(2,counter):coor(4,counter))+ noiseVolume.*mask(coor(1,counter):coor(3,counter),coor(2,counter):coor(4,counter)))./mask(coor(1,counter):coor(3,counter),coor(2,counter):coor(4,counter));
            else
                noisyVolume(coor(1,counter):coor(3,counter),coor(2,counter):coor(4,counter))= volume(coor(1,counter):coor(3,counter),coor(2,counter):coor(4,counter)) + noiseVolume;
            end;
            counter = counter + 1;
        end;
    else
        varNoise = volumeSTD^2 / SNR;
    
        noiseVolume= random('norm',volumeMean,sqrt(varNoise),volumeSize);
        m = find(mask == 0);
        mask(m) = 1;
        
        if(normalisation)
            noisyVolume= (volume+ noiseVolume.*mask)./mask;
        else
                noisyVolume= volume + noiseVolume;
        end;
    end;
end;


