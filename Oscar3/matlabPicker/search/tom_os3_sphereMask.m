%%
%just a wrapper for tom_spheremask
%returns the sphere mask and number of voxels > 0
function [mask maskSize] = tom_os3_sphereMask(template,blurring);

    templateSize = size(template);
% 
    if( (length(templateSize) == 3 && ~isequal(templateSize(1), templateSize(2), templateSize(3)) ) || ...
        (length(templateSize) == 2 && ~isequal(templateSize(1), templateSize(2) )))
        mask = ones(templateSize);
        maskSize = numel(mask);
        return;
    end;
    
    if(~exist('blurring'))
        if(size(template,1) < 50)
            blurring = floor(size(template,1)/10);
        else
            blurring = 5;
        end;
    end;
    blurring = 0;
    radius = floor(min(size(template))/2)-blurring-1;
    center = floor(size(template)/2) +1;
    if(length(center) == 2)
        center = [center 1];
    end;

    mask = single(sphMask(ones(size(template)), ...
                            radius, ...
                            blurring, ...
                            center));


    maskSize = length(find(mask > 0));
    %maskSize = numel(mask);
    
    
%%    
function mask = sphMask(vol,radius,sigma,center)

    if (nargin < 4)
        center=[floor(size(vol,1)/2)+1, floor(size(vol,2)/2)+1, floor(size(vol,3)/2)+1];
    end;
    [x,y,z]=ndgrid(0:size(vol,1)-1,0:size(vol,2)-1,0:size(vol,3)-1);
    if (nargin < 2)
        radius = floor((min(min(size(vol,1),size(vol,2)),size(vol,3))-1)/2) ;
    end;
    x = sqrt((x+1-center(1)).^2+(y+1-center(2)).^2+(z+1-center(3)).^2);
    ind = find(x>=radius);
    clear y z;
    mask = ones(size(vol,1), size(vol,2), size(vol,3));

    mask(ind) = 0;
    if (nargin > 2) 
        if (sigma > 0)
            mask(ind) = exp(-((x(ind) -radius)/(0.5*sigma)));
            ind = find(x >= radius + sigma);
            mask(ind) = 0;
        end;
    end;
 
   
   