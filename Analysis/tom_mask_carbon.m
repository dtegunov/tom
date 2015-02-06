function mask = tom_mask_carbon(I)
%CREATE_MICROGRAPH_MASK Create mask to remove dark strips and other noise
% NOTE:
%   Input and output will have same dimensions. If input is not 128x128,
%   then it is downsampled to this resolution, and the resulting mask
%   upsampled again.

N = size(I, 1);
scale = N / 128;
if (N ~= 128) 
    I = imresize(I, 1/scale);
end

edges = create_edge_probability_map(I) > 0.2;

se = strel('disk', 4);
se2 = strel('disk', 10);
area_threshold = 100;

A = imdilate(edges, se);
% figure; imagesc(A); axis image; colormap gray;

% remove small components
cc = bwconncomp(A);
    for j = 1:length(cc.PixelIdxList)
        if length(cc.PixelIdxList{j}) < area_threshold
            A(cc.PixelIdxList{j}) = false;
        end
    end

% figure; imagesc(A); axis image; colormap gray;

B = imdilate(A, se2);
% figure; imagesc(B); axis image; colormap gray;

C = ~A;
cc = bwconncomp(C);
for index = 1:cc.NumObjects
    selection = B(cc.PixelIdxList{index});
    if all(selection) && size(selection, 1) > 0
        A(cc.PixelIdxList{index}) = true;
    end
end

% mask = ~A;
D = ~A;

% Try to remove components of the mask that correspond to carbon strips
cc = bwconncomp(D);
stats = regionprops(cc, 'Area', 'ConvexArea');
areas = [stats.Area];
convex_areas = [stats.ConvexArea];
[~, largest_indices] = sort(areas, 2, 'descend');
n = min(5, length(largest_indices));

area_ratios = zeros(n, 1);
for i = 1:n
    index = largest_indices(i);
    area_ratios(i) = areas(index) / convex_areas(index);
    if area_ratios(i) < 0.8     % probably part of carbon strip
        D(cc.PixelIdxList{index}) = false;
    end
end
% disp(area_ratios);


mask = D;

% Resize output to same size as input

if (N ~= 128) 
    mask = imresize(mask, scale) > 0.5;
end

end


function I3 = create_edge_probability_map(I)
%CREATE_EDGE_PROBABILITY_MAP Probability of each pixel of being an edge
%
%   NOTE:
%       Input image is assumed to be 128x128
%
%   Just for convenience, copy of part of script 'n08d26_edge_detection'.

% Create stack of rotated filters

N = size(I, 1);

sigma = 2;

n = 15*sigma;
ratio = 2.5;
thetas = linspace(0, 180, 13);

f1 = zeros(n);
f2 = zeros(n);
for i = 1:length(thetas)
    theta = thetas(i);
    [f1(:, :, i) f2(:, :, i)] = oriented_filter_pair(n, sigma, ratio, theta);
end

% Apply all filters to image

% extra:
% directions = false(N);

images = zeros(N);
for i = 1:length(thetas)
    y1 = imfilter(I, f1(:, :, i), 'symmetric', 'same');
    y2 = imfilter(I, f2(:, :, i), 'symmetric', 'same');
    y = y1.*y1 + y2.*y2;
    images(:, :, i) = y;
    
    % extra
%     directions(:, :, i) = y2 > 0;
end

[I, indices] = max(images, [], 3);

% % extra
% direction = false(N);
% for i = 1:N
%     for j = 1:N
%         direction(i, j) = directions(i, j, indices(i, j));
%     end
% end

% Apply oriented non-maximal suppression

I2 = zeros(size(I));
is_vertical = mod(thetas(indices) + 45, ones(size(I))*180) < 90;
N = size(I, 1);
for i = 2:N-1
    for j = 2:N-1
        if is_vertical(i, j)
            if I(i, j) >= I(i, j-1) && I(i, j) >= I(i, j+1)
                I2(i, j) = I(i, j);
            end
        else
            if I(i, j) >= I(i-1, j) && I(i, j) >= I(i+1, j)
                I2(i, j) = I(i, j);
            end
        end
    end
end

% Convert to probability-like number
% Assume image is 128x128
intensity_sigma = 1.5e7;

I3 = 1 - exp(-I2/intensity_sigma);

end


function [f1, f2] = oriented_filter_pair(n, sigma, ratio, theta)
%ORIENTED_FILTER_PAIR Create pair of oriented filters to detect edges
%   Creates pair of filters for detecting edges. In the x direction, 'f1'
%   is the second Gaussian derivative and 'f2' its Hilbert transform.
%   In the y direction is just a normal gaussian with a larger sigma
%   (=sigma*ratio). This is if theta=0, otherwise its just rotated by
%   theta.

scale = 10;

x = gaussian_function(2, sigma*scale, n*scale);
h = imag(hilbert(x));
y = gaussian_function(0, ratio*sigma*scale, n*scale);

f1 = repmat(x, [n*scale 1]) .* repmat(y', [1 n*scale]);
f2 = repmat(h, [n*scale 1]) .* repmat(y', [1 n*scale]);

f1 = imrotate(f1, theta, 'crop');
f2 = imrotate(f2, theta, 'crop');

f1 = imresize(f1, 1/scale);
f2 = imresize(f2, 1/scale);

end


function x = gaussian_function(derivative, sigma, size)
%GAUSSIAN Return sampled gaussian, or derivative of it
%   INPUT:
%       derivative      0, 1 or 2
%       sigma           3 for example
%       size            5*sigma + 1for example
%   
%   For derivative = 0, returns sampled exp(-x^2/(sigma^2)), so no normalisation.

t = (1:size) - (size+1)/2;    % symmetric about 0, step size = 1

x = [];

% s2 = sigma*sigma;
t = t/sigma;

if derivative == 0
    x = exp(-(t.*t));
end
if derivative == 1
    x = exp(-(t.*t))*2.*t;
end
if derivative == 2
    x = exp(-(t.*t)).*(-2+4*t.*t);
end

end


