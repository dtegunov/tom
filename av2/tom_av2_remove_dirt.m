function image = tom_av2_remove_dirt(image,mask,dev)

if nargin < 3
    dev = 1;
end

image_mask = image.*abs(1-mask);

[mean max min std] = tom_dev(image_mask);

idx = find(image_mask<mean-dev.*std);

image(idx) = mean;