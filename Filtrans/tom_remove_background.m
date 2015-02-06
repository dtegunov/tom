function [im_out, bg] = tom_remove_background(image,nofrequencies)
%TOM_REMOVE_BACKGROUND creates ...
%
%   [im_out, bg] = tom_remove_background(image,nofrequencies)
%
%PARAMETERS
%
%  INPUT
%   image               ...
%   nofrequencies       ...
%  
%  OUTPUT
%   im_out              ...
%   bg                  ...
%
%EXAMPLE
%   ... = tom_remove_background(...);
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




bg = tom_bandpass(image,0,nofrequencies,10);

im_out = image - bg;

