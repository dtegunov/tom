function [ deconv ] = tom_deconv_tomo( vol, angpix, defocus, snrfalloff, deconvstrength, highpassnyquist, phaseflipped )

% Parameters:
%
% vol:
% tomogram volume (or 2D image)
%
% angpix:
% angstrom per pixel
%
% defocus:
% defocus in micrometers, positive = underfocus
%
% snrfalloff:
% how fast does SNR fall off, i. e. higher values will downweight high frequencies; 
% values like 1.0 or 1.2 seem reasonable
%
% deconvstrength:
% how much will the signal be deconvoluted overall, i. e. a global scale for SNR; 
% exponential scale: 1.0 is SNR = 1000 at zero frequency, 0.67 is SNR = 100, and so on
%
% highpassnyquist:
% fraction of Nyquist frequency to be cut off on the lower end (since it will be boosted the most)
%
% phaseflipped:
% whether the data are already phase-flipped
%
%
%
% Usage example:
% deconv = tom_deconv_tomo(mytomo, 3.42, 6, 1.1, 1, 0.02, false);
%

highpass = 0:1/2047:1;
highpass = min(1, highpass./highpassnyquist).*pi;
highpass = 1-cos(highpass);

snr = exp((0:-1/2047:-1).* snrfalloff.* 100./ angpix).* (10^(3 * deconvstrength)).* highpass;
ctf = tom_ctf1d(2048, angpix*1e-10, 300e3, 2.7e-3, -defocus*1e-6, 0.07, 0,0);
if phaseflipped
    ctf = abs(ctf);
end;
wiener = ctf./(ctf.*ctf+1./snr);

s1 = -floor(size(vol,1)/2);
f1 = s1 + size(vol,1) - 1;
s2 = -floor(size(vol,2)/2);
f2 = s2 + size(vol,2) - 1;
s3 = -floor(size(vol,3)/2);
f3 = s3 + size(vol,3) - 1;

[x, y, z] = ndgrid(s1:f1,s2:f2,s3:f3);
x = x./abs(s1);
y = y./abs(s2);
z = z./max(1, abs(s3));
r = sqrt(x.^2+y.^2+z.^2);
r = min(1, r);
r = ifftshift(r);

x = 0:1/2047:1;

ramp = interp1(x,wiener,r);

deconv = real(ifftn(fftn(single(vol)).*ramp));

end