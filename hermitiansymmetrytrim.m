function trimmed = hermitiansymmetrytrim(z)

n = size(z);
trimmed = z(1:floor(n(1)/2)+1,:,:);