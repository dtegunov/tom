function [ stddev, avg ] = tom_localstd( map, radius )

mask = tom_spheremask(ones(size(map)),radius);
maskinner = tom_spheremask(ones(size(map)),2);
mask = mask - maskinner;
masksum = sum(mask(:));
mask = fftn(fftshift(mask));
mapft = fftn(map);
map2ft = fftn(map.^2);

conv1 = real(ifftn(mapft.*conj(mask)));
conv2 = real(ifftn(map2ft.*conj(mask)));

stddev = sqrt(max(0, masksum.*conv2 - conv1.^2))./masksum;
avg = conv1./masksum;

end

