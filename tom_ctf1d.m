function [ ctf, acurve, pcurve ] = tom_ctf1d( length, pixelsize, voltage, cs, defocus, amplitude, phaseshift, bfactor )

ny = 1 / pixelsize;
lambda = 12.2643247 / sqrt(voltage * (1.0 + voltage * 0.978466e-6)) * 1e-10;
lambda2 = lambda * 2;

points = 0:length-1;
points = points./(2 * length).*ny;
k2 = points.^2;
term1 = lambda.^3.*cs.*k2.^2;

w = pi./2.*(term1 + lambda2.*defocus.*k2) - phaseshift;

acurve = cos(w).*amplitude;
pcurve = -sqrt(1 - amplitude.^2).*sin(w);
bfactor = exp(-bfactor.*k2.*0.25);
ctf = (pcurve + acurve).*bfactor;


end

