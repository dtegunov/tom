%TOM_MEX_EMWRITE_PASTE Writes a volume to EM-file by pasting it into
% an existing EM-file.
% [size, magic, comment, emdata, userdata] =
%   tom_mex_emwrite_paste(filename, volume, [subregion_start, reverse_sampling, allow_conversion]);
%
%INPUT:
%  filename: the name of the em-file to read. It must already exist.
%  volume: The volume which will be saved to file.
%  subregion_start: The zero based index of the first voxel where to start writing.
%    It is a 3-vector. Setting to [] is the same as [0,0,0].
%  reverse_sampling: Optional parameter to not paste the volume contigously
%    into the emfile. Setting to [] is the same as [1,1,1]. Otherwise it must
%    be an integer 3-vector saying It is the inverse of the sampling in TOM_EMREADC3
%  allow_conversion: Optional parameter which can be true, false or []
%    Omitting this parameter or setting to [] defaults to true. If false,
%    the data-type of the volume must correspond to the type saved in the
%    em-file or an error occures. If true, conversions are allowed.
%
%OUTPUT:
%  size: The original size of the volume as saved in the em-file after
%   appending the volume.
%  magic: The first 4 bytes from the em-header.
%  comment: Data from the emheader (see TOM_EMREAD)
%  emdata: Data from the emheader (see TOM_EMREAD)
%  userdata: Data from the emheader (see TOM_EMREAD)
%
% This is a mex-function for the C-function defined in io.c
%
% The function checks all input parameters and allocates all memory before
% starting to write the file. The only error which can happen then, is an
% IO-error. The header of the EM-file is not touched.
%
%SEE ALSO
%  TOM_EMREAD, TOM_EMREADC3.
%
%   created by Thomas Haller Jan. 22 2008
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

error('This m-file can not be executed. Compile the underlying mex-routine.');
