function [ factors ] = tom_shiftfactors( dims, delta )

if size(delta,1) > 1
    delta = delta';
end

dimx = dims(1);
dimy = 1;
dimz = 1;

if (numel(dims) > 1)
    dimy = dims(2);
end;
if (numel(dims) > 2)
    dimz = dims(3);
end;

[x,y,z]=ndgrid( -floor(dimx/2):-floor(dimx/2)+(dimx-1),...
                -floor(dimy/2):-floor(dimy/2)+dimy-1, ...
                -floor(dimz/2):-floor(dimz/2)+dimz-1);
indx = find([dimx,dimy,dimz] == 1);
delta(indx)=0;
delta = delta./[dimx dimy dimz];
x = delta(1)*x + delta(2)*y + delta(3)*z; clear y; clear z;

factors = exp(-2*pi*i*x);

end