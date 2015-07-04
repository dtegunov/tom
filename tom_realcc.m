function [ score ] = tom_realcc( im1, im2, mask )

if nargin < 3
    mask = ones(size(im1));
end;

im1 = tom_norm(im1, 'mean0+1std', mask);
im2 = tom_norm(im2, 'mean0+1std', mask);

score = im1.*im2;
score = sum(score(mask == 1));
score = score / sum(sum(mask==1));

end

