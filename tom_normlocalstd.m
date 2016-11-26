function [ normalized ] = tom_normlocalstd( sharpened, unsharpened, localradius )

stdfilt = tom_localstd(sharpened, localradius);
stdunfilt = tom_localstd(unsharpened, localradius);

diff = ones(size(stdfilt));
diff(stdfilt>0) = stdunfilt(stdfilt>0)./stdfilt(stdfilt>0);
% diff = diff.^2;
normalized = sharpened.*diff;

end

