%TOM_MEX_EMWRITE Saves an em-file to disk.
% tom_mex_emwrite(filename, data, subregion_in, magic, comment, emdata, userdata);
%
%INPUT:
%  filename: the name of the em-file.
%  data: The data to be saved.
%  subregion_in: Can be either [] for saving the whole volume or must be a 6-vector
%     for saving only a subvolume (same format as in TOM_EMREADC3).
%  magic: The first 4 bytes of the header. If the datatype of the data differs from
%     what is given in magic, the data will be converted.
%  comment: See the definition of the EM-Header in TOM_EMREAD.
%  emdata: See the definition of the EM-Header in TOM_EMREAD.
%  userdata: See the definition of the EM-Header in TOM_EMREAD.
%
% This is a mex-function for the C-function defined in libtomc.
%
%SEE ALSO
%  TOM_EMREAD3, TOM_MEX_EMREAD
%
%   created by Thomas Haller Jan. 21 2008
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

