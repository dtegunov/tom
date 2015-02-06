function [val x y ] = tom_os3_max2d(plane)



    [val pos] = max(plane);
    
    [val y] = max(val);
    
    x = pos(y);