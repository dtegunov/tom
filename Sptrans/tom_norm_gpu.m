function [ map ] = tom_norm_gpu( map )

m = mean(mean(mean(map)));
s = std(map(:), 1);

map = (map - m)./s;

end

