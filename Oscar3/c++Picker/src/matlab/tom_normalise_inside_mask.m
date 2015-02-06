function v = tom_normalise_inside_mask(v, mask, truely_inside)


if (~exist('truely_inside','var'))
    truely_inside = false;
end;

v = double(v);

idx_mask = mask~=0;
nmask = sum(idx_mask(:));

t2 = v;
t2_mean = sum(t2(idx_mask(:))) / nmask;
t2 = t2 - t2_mean;

t2 = t2 .* mask;


if (truely_inside)
    ndivisor = nmask;
else
    ndivisor = numel(t2);
end;


t2_mean = sum(t2(idx_mask(:))) / nmask;
t2_variance = (sum(t2(idx_mask(:)).^2) - nmask*t2_mean^2) / ndivisor;
t2 = (t2 - t2_mean) / sqrt(t2_variance);
t2(~idx_mask) = 0;

v = t2;















