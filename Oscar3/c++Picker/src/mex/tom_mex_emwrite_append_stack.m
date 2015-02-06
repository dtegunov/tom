%TOM_MEX_EMWRITE_APPEND_STACK Writes a volume to em-file by appending
% it along the 3th dimension.
% [size, magic, comment, emdata, userdata] =
%   tom_mex_emwrite_append_stack(filename, volume, [allow_conversion]);
%
%INPUT:
%  filename: the name of the em-file to read. It must already exist.
%  volume: The volume to append to the file.
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
% This is a mex-function for the C-function tom_io_em_write_append_stack 
% defined in tom_io.c
%
% The function checks all input parameters and allocates all memory before
% starting to write the file. The only error which can happen then, is an 
% IO-error. In that case the file may be damaged :)
% However, you can recover the em-file by truncating it at its original
% size, and check the Z-dimension in the em-header.
% 
%SEE ALSO
%  TOM_EMREAD
%
%   created by Thomas Haller Nov. 28 2007
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
