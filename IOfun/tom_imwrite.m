function tom_imwrite(image, filename, fmt)
%TOM_IMWRITE write image to graphics file
%
%   tom_imwrite(image, filename, fmt)
%
%PARAMETERS
%
%  INPUT
%   image               2D array
%   filename            filename of file
%   fmt                 format (e.g., 'TIFF')
%  
%  OUTPUT
%
%    IMWRITE(A,FILENAME,FMT) writes the image A to FILENAME.
%    FILENAME is a string that specifies the name of the output
%    file, and FMT is a string the specifies the format of the
%    file.  A can be either a grayscale image (M-by-N) or a
%    truecolor image (M-by-N-by-3)
%
%EXAMPLE
%   tom_imwrite(tom_norm(xxx, '3std'), 'xxx.tif', 'TIFF');
%   
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FF 01/02/03
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

image=tom_norm(image', 255);
imwrite(uint8(image), filename, fmt);