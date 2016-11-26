function [ averaged1, averaged2 ] = tom_prevent_divergence( h1, h2, nshells )

side = size(h1, 1);
mi = -side / 2;
ma = mi + side - 1;

[x, y, z] = ndgrid(mi:ma,mi:ma,mi:ma);
r2 = x.^2+y.^2+z.^2;
mask = r2 < nshells*nshells;

h1ft = fftshift(fftn(h1));
h2ft = fftshift(fftn(h2));

averaged1 = h1ft;
averaged1(mask) = (h1ft(mask) + h2ft(mask)).*0.5;
averaged2 = h2ft;
averaged2(mask) = (h1ft(mask) + h2ft(mask)).*0.5;

averaged1 = real(ifftn(ifftshift(averaged1)));
averaged2 = real(ifftn(ifftshift(averaged2)));

end