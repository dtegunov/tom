function filtered = tom_BandpassNeat( image, low, high )

padded = padarray(image, floor(size(image)./2), 'symmetric');
padded = tom_bandpass(padded, low*2, high*2);
filtered = tom_cut_out(padded,'center',size(image));

end

