%TOM_MEX_EMREAD Reads em-file from disk.
% [vol, size, magic, comment, emdata, userdata] =
%   tom_mex_emread(filename, [subregion, sampling, binning]);
%
%INPUT:
%  filename: the name of the em-file to read.
%  subregion: optional parameter for reading only a subregion of the
%     volume/image. If ommitted or set to an empty array, the entire volume
%     is read. Otherwise it must be a 6 vector containing the "first"
%     corner and the size of the subregion. The "first" corner
%     is given at position subregion(1:3) and is ZERO-based (i.e. the first
%     element of the volume has coorinate [0,0,0] despite the MATLAB
%     convention [1,1,1]). The size of the region is given in
%     subregion(4:6). The subregion must be entirely inside the volume. Thus if
%     you read an 2D-image subregion(3) must be 0 and subregion(5) must be
%     1.
%  sampling: optional parameter for sub-sampling. It can be an empty array
%     for no sampling. Otherwise it must be a 3-vector with the sampling
%     factor along each direction independently. 0 or 1 means no sampling.
%     A value of n means that every n-th value along the direction is
%     taken. If the subregion size is not a multiple of the sampling
%     factor, the last trailing elements are sampled too.
%  binning: optional parameter for binning. A binning factor of n means
%     that each n elements along the direction is averaged. If the
%     subregion size is not a multiple of the binning factor, the trailing
%     voxels are not taken.
%
%  All parameters can be combined. If sampling and binning is combined,
%  the volume is sampled first and then binned.
%  Notice that the sampling and binning factors are the number of voxels
%  to be combined along a certain dimension. The size of the resulting volume
%  is the size of the subregion divided by the sampling/binning factor.
%OUTPUT:
%  vol: The read volume/image. The data-type is the same as the original one,
%    except in case of binning, where doubles are returned.
%  size: The original size of the volume as saved in the em-file.
%  magic: The first 4 bytes from the em-header.
%  comment: Data from the emheader (see TOM_EMREAD)
%  emdata: Data from the emheader (see TOM_EMREAD)
%  userdata: Data from the emheader (see TOM_EMREAD)
%
% This is a mex-function for the C-function read_vol defined in tom_io.c
%
%SEE ALSO
%  TOM_EMREAD
%
%   created by Thomas Haller Nov. 09 2007
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

