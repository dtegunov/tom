function [ ctf, acurve, pcurve ] = tom_ctf2d( length, pixelsize, voltage, cs, defocus, astigmatism, angle, amplitude, phaseshift, bfactor )

ny = 0.5 / pixelsize;
voltagerel = voltage * (1 + voltage / 1022000);
lambda = sqrt(150.4 / voltagerel) * 1e-10;
lambda2 = lambda * 2;

[x, y] = ndgrid(-length:length-1, -length:length-1);
alpha = atan2(y, x);
points = sqrt(x.^2+y.^2);
points = points./(length).*ny;
k2 = points.^2;
term1 = lambda.^3.*cs.*k2.^2;

astigmatism = astigmatism / 2;
defocus = defocus + astigmatism.*cos((alpha - angle).*2);
w = pi./2.*(term1 + lambda2.*defocus.*k2) - phaseshift;

acurve = cos(w).*amplitude;
pcurve = -sqrt(1 - amplitude.^2).*sin(w);
bfactor = exp(-bfactor.*k2.*0.25);
ctf = (pcurve + acurve).*bfactor;


end

