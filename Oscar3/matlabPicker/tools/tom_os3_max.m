function [res coor] = tom_os3_max(img)

if(size(img,3) == 1 && size(img,2) == 1)
    %1d
    [res coor] = max(img);
    coor = [coor 1 1];
elseif(size(img,3) == 1)
    %2d
    [res x y] = tom_os3_max2d(img);
    coor = [x y 1];
else
    %3d
    [res x y z] = tom_os3_max3d(img);
    coor = [ x y z];
end;
