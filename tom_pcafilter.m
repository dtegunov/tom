function [ filtered, fractionexplained ] = tom_pcafilter( data, numcomponents )

[eigenvecs, eigenvals, latent, tsquared, explained, mu] = pca(data, 'NumComponents', numcomponents, 'Centered', true);

filtered = repmat(mu, [size(data, 1), 1]);
for v=1:numcomponents
    eigenvec = eigenvecs(:,v)';
    eigenvec = repmat(eigenvec, [size(data, 1), 1]);
    eigenvec = eigenvec.*repmat(eigenvals(:,v), [1, size(data, 2)]);
    filtered = filtered + eigenvec;
end;

fractionexplained = sum(explained(1:numcomponents)) / 100;

end

