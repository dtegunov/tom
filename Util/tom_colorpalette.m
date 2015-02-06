function rgb = tom_colorpalette(no_colors)
%TOM_COLORPALETTE creates ...
%
%   rgb = tom_colorpalette(no_colors)
%
%PARAMETERS
%
%  INPUT
%   no_colors           ...
%  
%  OUTPUT
%   rgb                 ...
%
%EXAMPLE
%   ... = tom_colorpalette(...);
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


hsv = ones(no_colors,3);
hsv(:,1) = ([1:no_colors]-1)'/no_colors;
hsv(:,2) = 0.5;
hsv(:,3) = 0.8;

rgb = hsv2rgb(hsv);
