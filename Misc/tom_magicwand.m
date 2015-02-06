function bin_mask = tom_magicwand(im, ylist, xlist, zlist,tolerance)
%TOM_MAGICWAND simulates the Photoshop's magic wand tool
%
%   bin_mask = tom_magicwand(im, ylist, xlist, zlist,tolerance)
%
%PARAMETERS
%
%  INPUT
%   im                  input image RGB
%   ylist               vector of row cordinates    (reference pixels)
%   xlist               vector of column cordinates (reference pixels)
%   zlist               ...
%   tolerance           distance to reference pixels
%  
%  OUTPUT
%   bin_mask            binary mask of selected regions
%
%EXAMPLE
%   The following code selects the girl's face:
%   im = imread('trees.tif');
%   bin_mask = tom_magicwand(im, [160], [220], 6);
%   subplot(1, 2, 1); imshow(im);
%   subplot(1, 2, 2); imshow(bin_mask);
%
%   for i=1:100;bin_mask = tom_magicwand(im, [160], [220],i);subplot(1, 2,1);
%       imshow(im);
%       subplot(1, 2, 2); imshow(bin_mask);drawnow;pause(.3);
%   end;
%
%   for i=.1:.05:1;
%       bin_mask = tom_magicwand(tom_spheremask(ones(64,64,64),16,15), [33], [33],  [33],i);
%       tom_dspcub(bin_mask);drawnow;
%   end;
%
% NOTES
%   * Tested with MATLAB R13 & Image Processing Toolbox
%   * C++ version by Daniel Lau and Yoram Tal are also available in MATLAB
%   Central
%
%REFERENCES
%
%SEE ALSO
%   ...
%
% (C) 2004 Son Lam Phung
% Email: s.phung@ecu.edu.au
% Edith Cowan University

H = size(im, 1); % image height
W = size(im, 2); % image width

% Check arguments
if any(ylist > H)
    s = sprintf('Row cordinates greater than the image height (%g).', H);
    error(s);
end

if any(xlist > W)
    s = sprintf('Column cordinates greater than the image height (%g).', W); 
    error(s);
end

c_r = double(im);

N = length(ylist); % Number of reference pixels

if size(im,3)==1

    % Find all pixels whose colors fall within the specified tolerance
    color_mask = false(H, W);
    for idx = 1:length(ylist)
        ref_r = double(im(ylist(idx), xlist(idx), 1));
        color_mask = color_mask | ((c_r - ref_r) .^ 2) <= tolerance ^ 2;
    end

    % Connected component labelling
    [objects, count] = bwlabel(color_mask, 8);

    % Initialize output mask
    bin_mask = false(H, W);

    % Linear indices of reference pixels
    pos_idxs = (xlist - 1) * H + ylist;

    for idx = 1:count
        object = (objects == idx); % an object

        % Add to output mask if the object contains a reference pixel
        if any(object(pos_idxs))
            bin_mask = bin_mask | object;
        end
    end
else
    L= size(im, 3); % image width
    % Find all pixels whose colors fall within the specified tolerance
    color_mask = false(H, W,L);
    for idx = 1:length(ylist)
        ref_r = double(im(ylist(idx), xlist(idx), zlist(idx)));
        color_mask = color_mask | ((c_r - ref_r) .^ 2) <= tolerance ^ 2;
    end

    % Connected component labelling
    [objects, count] = bwlabeln(color_mask, 18);

    % Initialize output mask
    bin_mask = false(H, W, L);

    % Linear indices of reference pixels
    pos_idxs = (zlist - 1) * H * W +(xlist - 1) * H + ylist;

    for idx = 1:count
        object = (objects == idx); % an object

        % Add to output mask if the object contains a reference pixel
        if any(object(pos_idxs))
            bin_mask = bin_mask | object;
        end
    end

end;
