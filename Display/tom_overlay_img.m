function h=tom_overlay_img(img,mask,weight,colour_rgb,figure_h)
%TOM_OVERLAY_IMG overlays an image with a mask 
%
%   h=tom_overlay_img(img,mask,weight,colour_rgb,figure_h)
%
%PARAMETERS
%
%  INPUT
%   img               image
%   mask              mask 2 overlay
%   colour_rgb        colour for overlay mask
%   weight            weight for transperency  
%   figure_h          figure handle
%
%
%EXAMPLE
%
% my_img=rand(64,64);
% mask=tom_spheremask(ones(64,64),12);
%
% tom_overlay_img(my_img,mask,0.4,[0 1 0]);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN/FB 01/24/06
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


if (nargin < 3)
    weight=0.3;
end;


if (nargin < 4)
    colour_rgb=[1 0 0];
end;


if (nargin < 5)
    figure_h='';
end;



if (isempty(figure_h)==0)
    figure(figure_h);
end;

mask=mask.*weight;
mask(find(mask==0))=1;

col_img(:,:,1)=ones(size(mask)).*colour_rgb(1);
col_img(:,:,2)=ones(size(mask)).*colour_rgb(2);
col_img(:,:,3)=ones(size(mask)).*colour_rgb(3);

imagesc(col_img); 

hold on; h=imagesc(img'); axis image; colormap gray;  hold off;

set(h, 'AlphaData', mask);



