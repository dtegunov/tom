function vol = tom_emreadc3(varargin)
%TOM_EMREADC3 reads data in EM-file format.
%
%vol = tom_emreadc3(filename, [subregion, sampling, binning])
%
%  Reads an EM-Image from file.
%  see TOM_EMREAD for a documentation of the EM-format.
%  The input parameters are the same as for TOM_MEX_EMREAD,
%  and this function acts as a wrapper.
%
%  Returns a structure containing the volume/image.
%  It is nearly the same as returned by TOM_EMREADC
%
%  The parameters subregion, sampling and binning are optional and default
%  to []. Otherwise subregion must be a 6-vector and sampling and binning a
%  3-vector. 
%  Subregion has the values [first_corner_x, first_corner_y,
%  first_corner_z, size_x, size_y, sizez]. The first corner is the zero
%  based index of the volume. All parameters are relative to the original
%  volume size.
%  Sampling is for each direction the integral number how many voxels are 
%  skiped. Binning is the number along each direction how many voxels are 
%  combined.
%  All parameters can be combined. In case of sampling and binning at the
%  same time, the volume is sampled first and than binned over the sampled
%  values.
%
%  Example:
%  To read the first slice (an 2D image) of size 100x100x1, with first
%  sampling if 3 and a binning of 2 (i.e. combine 2 voxels).
%  v = tom_emreadc3(filename, [0,0,0, 100,100,1], [3,3,1], [2,2,1]);
%
%
%REFERENCES
%
%SEE ALSO
%   TOM_EMREAD, TOM_MEX_EMREAD
%
%   created by Thomas Haller Nov. 13 2007
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

select_file = nargin < 1;
if (select_file)
    [filename, pathname] = uigetfile({'*.*';'*.em';'*.vol';'*.norm';'*.mask'}, 'Pick an EM-file');
    if (isequal(filename,0) || isequal(pathname,0))
        return;
    end;
    varargin{1} = [pathname filename];
end;

if (~ischar(varargin{1}))
    error('The first parameter must be the filename.');
end;

varargin{1} = getAbsFilename(varargin{1}, ~select_file);

[vol, em_size, magic, comment, emdata, userdata] = tom_mex_emread(varargin{:});
%vol = ones(1,1,1); em_size=size(vol); magic=zeros(1,4); comment=ones(1,80); emdata=zeros(1,40); userdata=zeros(1,256);
%n = 50000; tic; for i=(1:n) s = tom_emreadc3('/fs/pool/pool-bmsan/haller/DA/data/vols/testvol.em', [],[1000,1000,1000],[0,0,1]); end; t = toc; t/n



vol = tom_emheader_data_2_struct(vol, em_size, magic, comment, emdata, userdata, varargin{1});








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = getAbsFilename(s, call_which)


[success, message] = fileattrib(s);
if (~success)
    if (call_which)
        fid = fopen(s);
        if (fid ~= -1) 
            s = fopen(fid);
            fclose(fid);
        else
            call_which = false;
        end;
    end;
    if (~call_which)
        error(message);
    end;
else
    if (message.directory)
        error(['The given name is not a file, but a directory. (' s ')']);
    end;
    s = message.Name;
end;










