function out=tom_filter2resolution(in,pixsize_in_Angstrom,cutoff_in_Angstrom)
%TOM_FILTER2RESOLUTION filters a volume with known object-pixelsize to a
%cutoff resolution limit.
%
%   out=tom_filter2resolution(in,pixsize_in_Angstrom,cutoff_in_Angstrom)
%
%PARAMETERS
%
%  INPUT
%   in                  input 3D volume
%   pixsize_in_Angstrom objectpixelsize in Angstrom
%   cutoff_in_Angstrom  resolution cutoff in Angstrom
%  
%  OUTPUT
%   out                 filtered volume
%
%EXAMPLE
%REFERENCES
%
%SEE ALSO
%
%   created by SN 24/09/10 
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

nyquist=pixsize_in_Angstrom.*2;
mask_radius=(nyquist./cutoff_in_Angstrom).*size(in,1)./2;
%mask=tom_spheremask(ones(size(in)),mask_radius,0,[size(in)./2+1]);
mask=tom_spheremask(ones(size(in)),mask_radius);
out=tom_apply_weight_function(in,mask);
