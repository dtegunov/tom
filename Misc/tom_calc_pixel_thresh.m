function [thresh pix_nr im_back]=tom_calc_pixel_thresh(im,num_of_pixel,accuracy)
%TOM_CALC_PIXEL_THRESH creates ...
%
%   [thresh pix_nr im_back]=tom_calc_pixel_thresh(im,num_of_pixel,accuracy)
%
%PARAMETERS
%
%  INPUT
%   im                  ...
%   num_of_pixel        ...
%   accuracy            ...
%  
%  OUTPUT
%   thresh              ...
%   pix_nr              ...
%   im_back             ...
%
%EXAMPLE
%   .. = tom_calc_pixel_thresh(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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

pix_nr=size(im,1).*size(im,2);

lauf=max(max(im));


%while pix_nr < num_of_pixel
for i=1:10000 

    lauf=lauf-accuracy;
    pix_nr=sum(sum(im>=lauf));
    if (pix_nr > num_of_pixel)
        break;
    end;

end;
thresh=lauf;

im_back=(im>=thresh).*im;