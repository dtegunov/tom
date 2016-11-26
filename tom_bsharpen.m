function [ sharpened ] = tom_bsharpen( map, angpix, bfac, mtfslope )

m = -size(map,1) / 2;
p = size(map,1) / 2 - 1;
[x, y, z] = ndgrid(m:p, m:p, m:p);
r = sqrt(x.^2+y.^2+z.^2)./size(map,1);

mtf = 1 + r.*mtfslope;

r = (r./angpix).^2;
r = exp(0.25.*bfac.*r).*mtf;

map = fftshift(fftn(map));
map = map./r;
sharpened = real(ifftn(ifftshift(map)));

end

