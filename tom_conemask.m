function [ map ] = tom_conemask( map, direction, coneangle )

side = size(map, 1);
[xx, yy, zz] = ndgrid(-side/2:side/2-1,-side/2:side/2-1,-side/2:side/2-1);
vx = [xx(:), yy(:), zz(:)];

dotprod = vx*direction';
posondir = repmat(direction, [numel(xx), 1]).*repmat(dotprod, [1 3]);
conewidth = tand(coneangle).*abs(dotprod);
postocone = posondir - vx;
distfromcone = sqrt(sum(postocone.^2, 2))./max(10, conewidth)./1.2;
weight = max(0, min(1, distfromcone));
weight = (cos(weight * pi) + 1) / 2;

map = reshape(weight, side, side, side);

end

