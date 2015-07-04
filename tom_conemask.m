function [ map ] = tom_conemask( map, direction, coneangle )

side = size(map, 1);
[xx, yy, zz] = ndgrid(0:side-1,0:side-1,0:side-1);
xx = xx - side / 2;
yy = yy - side / 2;
zz = zz - side / 2;
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

