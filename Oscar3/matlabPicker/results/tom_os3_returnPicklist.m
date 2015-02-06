function [pickList align2d] = tom_os3_returnPicklist(peaks,angles,job,n)
%tom_os3_returnPicklist
%
%   tom_os3_returnPicklist
%
%PARAMETERS
%
%  INPUT
%       peaks
%       angles
%       job
%       n 
%  OUTPUT
%       picklist
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
%set the radius of the peak area - this area will be set to zero after one
%of the peaks has been selected

templateSize = job.templateSize;
binning = job.options.modifications.binning;
volumeName = job.volumeFile;


templateSize = templateSize / 2^binning;

if(length(templateSize) >2 && templateSize(3) <1)
    templateSize(3) = 1;
    
end;
radius = ceil(templateSize(1)/2);



pickList = {};
if(~isequal(size(peaks),templateSize))
    peaks = tom_os3_eraseBorders(peaks,zeros(templateSize),binning);
end;

counter = 1;
while(counter <= n && std(peaks(:))>0)
    
    [coordinates value peaks] = tom_peakc(peaks,radius);

    
    
    if(length(coordinates) == 2)
        %if is 2d template
        coordinates(3) = 1;
    end;

    if(iscell(angles))
        a = angles.value;
        angl = a(coordinates(1),coordinates(2),coordinates(3));
    else
        angl = angles(coordinates(1),coordinates(2),coordinates(3));
    end;
    coordinates = coordinates*2^binning;
    
    if(size(peaks,3) == 1)
        %if is 2d template
        coordinates(3) = 1;
    end;
    
    pick.coordinates  = coordinates;
    
    pick.angle        = angl;
    pick.value        = value;
    pick.templateSize = templateSize;
    pick.filename     = volumeName;    
    %pick.job = job;
    pickList{counter} = pick;

    counter = counter+1;
end;


align2d = tom_os3_pickList2Align2d(pickList,binning);