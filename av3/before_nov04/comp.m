function out = tom_compare(vol1, vol2, num)
% TOM_COMPARE comparison of 2 fourier transforms
%
% results = tom_compare(vol1, vol2, num)
%
% The Fourier space is subdivided into num rings  of  spherical-shell
% zones.  For  each  zone  values  are  calculated  which  serve  as
% measures for the agreement of Fourier coefficients, e.g., the mean
% square deviation.
%
%  PARAMETERS
%   VOL1        VOLUME 1 - in real space!
%   VOL2        VOLUME 2 - in real space!
%   NUM         Number of shells in reciprocal space
%
%   RESULTS     Array of dimension 10xNUM
%
% Oct. 18, 2003
% M. Riedlberger
% 

if (ndims(vol1) ~= ndims(vol2));
   error('Number of dimensions of input volumes do not equal!');
   return;
end;
if (size(vol1) ~= size(vol2));
   error('Dimensions of input volumes do not equal!');
end;
if (size(vol1,1) ~= size(vol1,2) || size(vol1,3) ~= size(vol1,2));
   error('x, y and z dimensions must be equal!');
end;


s = size(vol1,1);
f1 = fftshift(fftn(vol1))/(sqrt(s^3));
f2 = fftshift(fftn(vol2))/(sqrt(s^3));
f12 = f1-f2;
f1 = single(f1.*conj(f1));
f2 = single(f2.*conj(f2));
f12 = single(f12.*conj(f12));
out = single(zeros(num,11));
compare(f1, f2, f12, num, 1, num, out);
