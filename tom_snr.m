function [ im ] = tom_snr( im, snr, mask )

if nargin < 3
    mask=ones(size(im));
end;
idx = find(mask>0);

vi = var(im(idx));
noisetarget = vi/snr;

noise = randn(size(im));
noise = noise - mean(noise(:));
vn = var(noise(:));
noise = noise.*sqrt((noisetarget / vn));

im = im + noise;

end

