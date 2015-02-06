function [ ctf ] = tom_ctf1d( length, pixelsize, voltage, cc, cs, defocus, amplitude, decayksquared, decaycohill, decayspread )

ny = 0.5 / pixelsize;
voltagerel = voltage * (1 + voltage / 1022000);
lambda = sqrt(150.4 / voltagerel) * 1e-10;
lambda2 = lambda * 2;
ccvoltage = cc / voltage;
decayksquared = decayksquared * ny;

points = 0:length;
points = points./(length).*ny;
k = points;
k2 = points.^2;
term1 = lambda.^3.*cs.*k2.^2;
if decayksquared == 0
    decayksquared = 1;
else
    decayksquared = exp(-(k2./decayksquared.^2));
end;

e_i_k = 1;
if decaycohill ~= 0
    q0 = decaycohill./(cs.*lambda.^3).^0.25;
    e_i_k = exp(-(pi * pi).*q0.^2.*(cs * lambda.^3 * k.^3 - defocus.*lambda.*k).^2);
end;

e_e_k = 1;
if decayspread ~= 0
    delta_z = ccvoltage.*decayspread;
    e_e_k = exp(-(pi.*delta_z.*lambda2.*k2));
end;

w = pi./2.*(term1 - lambda2.*defocus.*k2);

ctf = e_e_k.*e_i_k.*decayksquared.*(sqrt(1 - amplitude.^2).*sin(w) + cos(w).*amplitude);

end

