function [ ctf, acurve, pcurve ] = tom_ctf2d( length, pixelsize, pixeldelta, pixelangle, voltage, cs, defocus, defocusdelta, defocusangle, amplitude, phaseshift, bfactor )

lambda = 12.2643247 / sqrt(voltage * (1.0 + voltage * 0.978466e-6)) * 1e-10;
lambda2 = lambda * 2;

[x, y] = ndgrid(-length:length-1, -length:length-1);
alpha = atan2(y, x);
pixeldelta = pixeldelta / 2;
pixelsize = pixelsize + pixeldelta.*cos((alpha - pixelangle).*2);
ny = 1./(pixelsize * 2 * length);
points = sqrt(x.^2+y.^2);
points = points.*ny;
k2 = points.^2;
term1 = lambda.^3.*cs.*k2.^2;

defocusdelta = defocusdelta / 2;
defocus = defocus + defocusdelta.*cos((alpha - defocusangle).*2);
w = pi./2.*(term1 + lambda2.*defocus.*k2) - phaseshift;

acurve = cos(w).*amplitude;
pcurve = -sqrt(1 - amplitude.^2).*sin(w);
bfactor = exp(-bfactor.*k2.*0.25);
ctf = (pcurve + acurve).*bfactor;


end

