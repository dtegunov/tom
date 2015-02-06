function [mpt, mval] = tom_mark_matchncc(im0, im1, center, w, max_disparity, im_shift, min_val, mutual)
%TOM_MARK_MATCHNCC finds best match using normalize cross correlation
% 
%  [mpt, mval] = TOM_MARK_MATCHNCC(im0, im1, center, w, 
%                                  max_disparity, im_shift, 
%                                  min_val, mutual)
%   Matches the region of size w around 
%   "center" of the image "im0" with the image "im1"
%   If max_disparity > 0, only in the region around
%   "center+im_shift" (regionsize = "2*max_disparity+1")
%   is searched.
%   Calles normxcorr2!
%   All parameters after center, can be ommittet or set to 'default'
%   to use the default. Look the source for details...
%
%PARAMETERS
%  INPUT
%     im0: Image from where the template is taken.
%     im1: Image where the template correlates.
%     center: Coordinates of the center of the template in 
%         im0. The center has the form [x, y] where im0(y,x) is the center.
%     w:  Sets the size of the used template around "center".
%         The finally used windowsize is always odd and integer and greater 
%         then 2;
%         Setting w to <0, 'default', [] or NaN, causes to use the standard
%         value. 
%         Setting it to a value in the range [0, 1], the actually
%         used window size is w = w*min(size(im0))-1).
%         Otherwise window w is used.
%         To obtain a valid value w is modified finally to w = bitor(fix(w),1).
%     max_disparity: if greater 0, only in a part of im1 is searched for
%         the template. To fast up the computation, if its known, that the
%         template in the new image is around a certain point.
%     im_shift: together with max_disparity it gives where in im1 the NCC
%         is computed (im_shift = [x y]). Precisely the window in im1 around
%         center+im_shift +- max_disparity is taken.
%         Setting im_shift to [] or leaving one component nan, is the same
%         as [0 0].
%         Can be computed using TOM_MARK_GETAVGSHIFT
%     min_val: If the correlation coeff. mval is not greater or equal than
%         min_val, the found correlated point is dropped.
%     mutual: Tries to refind the original point from the matched point, by
%         computing the NCC again, with the same parameters. If the points
%         not corrispond, no result is returned.
%  OUTPUT:
%     mpt: absolute coordinates ([x z]) of the best match in image im1.
%     mval: norm. coorelationcoeffizient of the match [-1..1].
%   If no point is found, mpt and mval are set to [].
%
%   See also TOM_MARK_FINDMATCHES, TOM_MARK_GETAVGSHIFT
%
%REFERENCES
%
%   16/03/07
%
%   created by Thomas Haller
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check and prepare input parameters.
mpt = [];
mval = [];

sim0 = size(im0);
sim1 = size(im1);

center = round(center);

if (~exist('w', 'var') || strcmp(w, 'default') || isempty(w) || ~(w>=0))
    w = 0.1;
end;
if (w <= 1)
    w = w*(min(size(im0))-1);
end;
w = bitor(fix(w), 1);
w_half = fix(w/2);
    
    
if (~exist('max_disparity', 'var') || strcmp(max_disparity, 'default') || numel(max_disparity)~=1 || ~(max_disparity>0))
    max_disparity = 0;
end;
max_disparity = ceil(max_disparity);

if (~exist('im_shift', 'var') || strcmp(im_shift, 'default') || numel(im_shift)~=2 || any(isnan(im_shift)))
    im_shift = [0 0];
else
    im_shift = round(im_shift);
end;


if (~exist('min_val', 'var') || strcmp(min_val, 'default') || isnan(min_val) || isempty(min_val)) 
    min_val = -2;
end;

if (~exist('mutual', 'var') || strcmp('mutual', 'default'))
    mutual = true;
end;

% do it...
[mpt, mval] = matchncc(im0, im1, sim0, sim1, center, w, w_half, max_disparity, im_shift, min_val, mutual);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the computation (with the clean, expected parameters,
% as prepeared from tom_mark_matchncc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mpt, mval] = matchncc(im0, im1, sim0, sim1, center, w, w_half, max_disparity, im_shift, min_val, mutual)

mpt = [];
mval = [];

if (any(center <= w_half) || any(sim0-center < w_half) || any(sim1<w))
    return;
end;



if (max_disparity > 0)
    im1center = center + im_shift;
    rect = [max(1, im1center-max_disparity-w_half); min(sim1([2,1]), im1center+max_disparity+w_half)];
    if (rect(2,1)-rect(1,1)+1 < w || rect(2,2)-rect(1,2)+1 < w)
        if (rect(2,2)-rect(1,2) < w || rect(2,1)-rect(1,1)+1 < w)
            return;
        end;
        error('TODO: CHANGE THIS CASE TO "return"?????');
    end;
    im1_window = im1(rect(1,2):rect(2,2),rect(1,1):rect(2,1));
else
    im1_window = im1;
end;


im0_window = im0(center(2)-w_half:center(2)+w_half, center(1)-w_half:center(1)+w_half);



[tmpt, tmval] = getBestNCCMatch( im0_window, im1_window, w_half);

if (isempty(tmpt))
    return;
end;
if (min_val > tmval) 
    return;
end;

if (max_disparity > 0) 
    tmpt = tmpt + rect(1,:) - 1;%% TODO
end;

if (mutual)
   [tmpt2, tmval2] = matchncc(im1, im0, sim1, sim0, tmpt, w, w_half, max_disparity, -im_shift, min_val, false);
   if (isempty(tmpt2) || (max(abs(tmpt2 - center)) > 1))
       return;
   end;
end;
mpt = tmpt;
mval = tmval;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matches the template of size (2*w_half+1) with 
% the image im (using normxcorr2).
% Returned is the imagecoordinate in im, where the 
% best match occured. value is the corralation-value
% there (between -1 and +1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos, value] = getBestNCCMatch(templ, im, w_half);

c = normxcorr2(templ, im);
c_old = c;
c = c(w_half+1 : size(c,1)-w_half, w_half+1 : size(c,2)-w_half);
[value, index] = max(c(:));    
[pos(2), pos(1)] = ind2sub(size(c), index);

c(pos(2), pos(1)) = -2;
[value2, index2] = max(c(:));    

if (value2 == value)
    warning('HALLO');
    pos = [];
    value = [];
end;
