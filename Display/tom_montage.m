function h = tom_montage(a, nMontageCols,range)
%TOM_MONTAGE creates a visualization of 3D image in a gallery
%
%   tom_montage(image_in, columns, rows)
%
%PARAMETERS
%
%  INPUT
%   image_in            3D image
%   columns             number of images per column
%   rows                number of images per rows
%  
%  OUTPUT
%   h           		...
%
%EXAMPLE
%   ... = tom_montage(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   %TOM_VOLXY, TOM_VOLXYZ, TOM_EMBROWSE, TOM_DSPCUB
%
%   created by AK 11/08/05
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

error(nargchk(0, 3, nargin, 'struct'))

%[a,cm, nMontageCols] = parse_inputs(varargin{:});


[nRows, nCols, nBands, nFrames] = size(a);

% Estimate nMontageColumns and nMontageRows given the desired ratio of
% Columns to Rows to be one (square montage).
%aspectRatio = 1; 
%nMontageCols = sqrt(aspectRatio * nRows * nFrames / nCols);

% Make sure montage rows and columns are integers. The order in the adjustment
% matters because the montage image is created horizontally across columns.
%nMontageCols = ceil(nMontageCols); 
nMontageRows = ceil(nFrames / nMontageCols);

% Create the montage image.
b = a(1,1); % to inherit type 
b(1,1) = 0; % from a
b = repmat(b, [nMontageRows*nRows, nMontageCols*nCols, nBands, 1]);

rows = 1 : nRows; 
cols = 1 : nCols;

for i = 0:nMontageRows-1
  for j = 0:nMontageCols-1,
    k = j + i * nMontageCols + 1;
    if k <= nFrames
      b(rows + i * nRows, cols + j * nCols, :) = a(:,:,:,k);
    else
      break;
    end
  end
end

if exist('range') ~= 1
  hh = imshow(b);
else
  hh = imshow(b,range);
end

if nargout > 0
    h = hh;
end
