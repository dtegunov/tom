function mask = tom_mask_carbon2(I,edge_thr,area_threshold,dil_cer)
%CREATE_MICROGRAPH_MASK Create mask to remove dark strips and other noise
% NOTE:
%   Input and output will have same dimensions. If input is not 128x128,
%   then it is downsampled to this resolution, and the resulting mask
%   upsampled again.

if (nargin < 2)
    edge_thr=[0.2 0.05 0.6 7];
end;

if (nargin < 3)
    area_threshold = 800;
end;

if (nargin < 3)
    dil_cer =2;
end;

if (nargin < 4)
    img_sz =[512 512];
end;

sz_org=size(I);
if (size(I,1)~=img_sz(1))
   I=imresize(I,img_sz);
end;
    
    
se = strel('disk', dil_cer);

%% filter org img
%edges = create_edge_probability_map(I) > edge_thr;
edges= edge(double(tom_filter(I,'aniso_linearPsi',1,2)),'canny',edge_thr(1),edge_thr(2)); % caron film

edges_sharp= imdilate(edge(double(I),'canny',edge_thr(3),edge_thr(4)),se); %carbon edge

if (std2(remove_small_comp(edges_sharp,area_threshold))==0)
    mask=ones(sz_org);
    return;
end;

A = imdilate(edges, se);
A = remove_small_comp(A,area_threshold);



%% check for obj without sharp edges 
tmp=tom_os3_eraseBorders(A,ones(2,2));
all_sz_org=get_comp_size(tmp);

tmp=tom_os3_eraseBorders(A.*(edges_sharp==0),ones(2,2));
[all_sz_sharp cc]=get_comp_size(tmp);

for i=1:length(all_sz_sharp)
    if isempty(find(all_sz_sharp(i)==all_sz_org)) ~=1
        A(cc.PixelIdxList{i})=0;
    end;
end;

A=remove_small_comp(A,area_threshold);

%% sharpen boarders

mask=zeros(size(A));

[comp cc]=get_comp_size(A);

for ii=1:length(comp)
    
    tmp=zeros(size(A));
    tmp(cc.PixelIdxList{ii})=1;
    
    tmpB = tom_os3_eraseBorders(tmp.*(edges_sharp==0),ones(2,2));
    
    [v cct]=get_comp_size(tmpB);
    tmpBsz_sort=sort(v,'descend');
    
    area_threshold=tmpBsz_sort(2);
    
    for j = 1:length(cct.PixelIdxList)
        if length(cct.PixelIdxList{j}) <= area_threshold
            tmpB(cct.PixelIdxList{j}) = false;
        end
    end
    mask = mask + imfill(imdilate(tmpB,se));
    
end;

mask=(mask+(edges_sharp))>0;
mask=imfill(mask,'holes');



%%  Resize output to same size as input
mask=imresize((mask==0),sz_org)>0.4;






function [v cc]=get_comp_size(I)
cc = bwconncomp(I);
for i=1:length(cc.PixelIdxList)
    v(i)=length(cc.PixelIdxList{i});
end;


function [I cc]=remove_small_comp(I,thr)

cc = bwconncomp(I);
for j = 1:length(cc.PixelIdxList)
    if length(cc.PixelIdxList{j}) < thr
        I(cc.PixelIdxList{j}) = false;
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




