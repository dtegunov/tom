function [out,n]=tom_hist3d(in)

% TOM_HIST3D calculates a histogram of 3D data.
%
%  [out,n]=tom_hist3d(in)
%
%  Calculates the histogram for 3D data, using the
%  matlab hist function.
%
% PARAMETERS
%
%  INPUT
%   in       3D array
%
%  OUTPUT
%   out      array with channel occupancy values
%   n        array with correspnoding values of each column
%
% EXAMPLE:
%
% [out,n]=tom_hist3d(in);
% plot(n,out);
%
%REFERENCES
%
% See also
%  HIST, HISTC
%
%   created by 02/01/03
%   updated by FF 03/28/03
%   updated by AK 05/23/06
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

error(nargchk(0, 1, nargin, 'struct'))


[out n]=hist(reshape(single(in),1,size(in,1).*size(in,2).*size(in,3).*size(in,4)),100);
% auf out: Kanalnr. vs. Haeufigkeit, auf n: Kanalnr. -> Value
