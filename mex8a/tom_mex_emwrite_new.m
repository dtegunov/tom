%TOM_MEX_EMWRITE_NEW Creates an new em-file with data set to zero
%  tom_mex_emwrite_new(filename, size, magic, [comment, emdata, userdata]);
%
%INPUT:
%  filename: the name of the em-file to read. It must already exist.
%  size: The size of the volume. Must be a 3-vector of positive integer.
%  magic: The first 4 bytes from the em-header. Must be int8 type.
%    It can also be a string describing the data-type used.
%    Possible values are 'double', 'single', 'int32', 'int16', 'int8', 'complex32', 'complex64'
%  comment: Optional data from the emheader (80 int8) (see TOM_EMREAD)
%  emdata: Optional data from the emheader (40 int32) (see TOM_EMREAD)
%  userdata: Optional data from the emheader (256 int8) (see TOM_EMREAD)
%
%EXAMPLE
%  tom_mex_emwrite_new('empty_file.em', [30,40,50], int8([6,0,0,9]), zeros(1,80,'int8'), zeros(1,40,'int32'), zeros(1,256,'int8');
%
%SEE ALSO
%  TOM_EMREAD
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

error('This m-file can not be executed. Compile the underlying mex-routine (maybe the mex-file is shadowed in the path).');


