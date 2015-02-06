function [subVolume vMean vSTD] = tom_os3_subVolumeStatistics(volume,coordinates,templateSize, mask)
%tom_os3_subVolumeStatistics reads a subVolume defined by coordinates out
%of a give volume. The volume borders are extended by templateSize/2+1 to
%avoid frequency (durch periodizitaet erzeugte artefakte am rand - name?)
%artefacts in later applied frequency domain applications.
%
%The statistics of the subVolume under the given mask are calculated by tom_os3_mean and std in
%frequency domain. 
%
%If the distance of the subvolume to the volume ends is less than
%templateSize/2, 
%
%  tom_os3_subVolumeStatistics(volume,coordinates,templateSize, mask)
%
%PARAMETERS
%
%  INPUT
%   volume      - the source volume
%   coordinates - coordinates of the subVolume in the volume
%   templateSize- the size of the template [x y z]
%   mask        - optional, if not given a sphere mask of template/2 radius
%                 is generated automaticaly
%  OUTPUT
%   subVolume
%   vMean
%   vSTD
% 
% 
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   tom_os3_subVolume, tom_os3_mean, tom_os3_std
%
%   created by TH2 07/11/07
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

%%
%   create a template appropriate mask
    if(~exist('mask'))
         mask = tom_os3_sphereMask(templateSize);
        %use mask below for testing 
%         mask = ones(templateSize(1),templateSize(2),templateSize(3),'single');
    end;

%%  calculate center of the subvolume relative to the subvolume    
    centerX = floor(length(coordinates(1):coordinates(2))/2);
    centerY = floor(length(coordinates(3):coordinates(4))/2);
    centerZ = floor(length(coordinates(5):coordinates(6))/2);
    
%   create a cube of same size as the later used subvolume, can be improved
    subVol  = volume(coordinates(1):coordinates(2), ...
                     coordinates(3):coordinates(4), ...
                     coordinates(5):coordinates(6));
    
%   extend the subvolume to calculate volume statistics under the (given) mask accurately             
    [subVolume pos]= tom_os3_subVolume([coordinates(1) coordinates(3) coordinates(5)]+[centerX centerY centerZ], ... 
                                        volume,zeros(size(subVol) + templateSize,'single'), ...
                                        'center');
                                    
%%  calculate volume statistics according to mask    
    vMean  = tom_os3_mean(subVolume,mask);
    vSTD   = tom_os3_std(subVolume,vMean,mask);          
                               
%%  paste volumes into zeros if neccessary to make sure the really used subVolume is in the center of the resulting volume                                
    if(coordinates(1) < floor(templateSize(1)/2)+1 || size(volume,1) - coordinates(2) < floor(templateSize(1)/2)+1 || ...
       coordinates(3) < floor(templateSize(2)/2)+1 || size(volume,2) - coordinates(4) < floor(templateSize(2)/2)+1 || ...
       coordinates(5) < floor(templateSize(3)/2)+1 || size(volume,3) - coordinates(6) < floor(templateSize(3)/2)+1 )                           
    
        lowerX = (coordinates(1) < floor(templateSize(1)/2)+1) * (floor(templateSize(1)/2)+1-coordinates(1));
        lowerY = (coordinates(3) < floor(templateSize(2)/2)+1) * (floor(templateSize(2)/2)+1-coordinates(3));
        lowerZ = (coordinates(5) < floor(templateSize(3)/2)+1) * (floor(templateSize(3)/2)+1-coordinates(5));
        
        
        upperX = (size(volume,1) - coordinates(2) < floor(templateSize(1)/2)+1) * (floor(templateSize(1)/2)+1 - (size(volume,1) - coordinates(2)));
        upperY = (size(volume,2) - coordinates(4) < floor(templateSize(2)/2)+1) * (floor(templateSize(2)/2)+1 - (size(volume,2) - coordinates(4)));
        upperZ = (size(volume,3) - coordinates(6) < floor(templateSize(3)/2)+1) * (floor(templateSize(3)/2)+1 - (size(volume,3) - coordinates(6)));
        
        pos = [lowerX+(lowerX == 0) ... 
               lowerY+(lowerY == 0) ...
               lowerZ+(lowerZ == 0)]; 
        
        subVolume = tom_paste(zeros(lowerX + size(subVolume,1) + upperX, ...
                                    lowerY + size(subVolume,2) + upperY, ...
                                    lowerZ + size(subVolume,3) + upperZ, 'single'),subVolume,pos+1);
        
        vMean     = tom_paste(zeros(lowerX + size(vMean,1) + upperX, ...
                                    lowerY + size(vMean,2) + upperY, ...
                                    lowerZ + size(vMean,3) + upperZ, 'single'),vMean,pos+1);
                                
        vSTD      = tom_paste(zeros(lowerX + size(vSTD,1) + upperX, ...
                                    lowerY + size(vSTD,2) + upperY, ...
                                    lowerZ + size(vSTD,3) + upperZ, 'single'),vSTD,pos+1);                        
    end;