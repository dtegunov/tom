function tom_emwritec3(varargin)
%TOM_EMWRITEC3 write data in EM-file format.
%
%tom_emwritec3([filename], data, [command], ...)
%
%  Saves the data as a new EM-Image. If the file already existed,
%  it will be overwritten.
%  see TOM_EMREAD for a documentation of the EM-format.
%  The function is a wrapper for TOM_MEX_EMWRITE and prepares
%  the EM-header for the mex-function.
%
%PARAMETERS
%filename: The name of the em-file. If ommited, the use can browse a file.
%data: The data to save to file.
%  This can be either a numeric array or a em-structure.
%  If a numeric array is passed, an empty header will be used (i.e. all set to 0).
%subregion: Specify to only save a subregion of the volume.
%  This is to cut out a subregion of the input volume (data), not to paste
%  it into an already existing EM-file (for that see TOM_EMWRITE_PASTE).
%  It can be either [] or a 6 vector containing the zero based coordinate of the "first"
%  voxel and the size of the subregion.
%type: A string specifying the data-type to save to em. Per default the
%  type of data (or data.Value) will be used. This is to force a conversion.
%  Possible values are 'int8', 'int16', 'int32', 'single', 'double'
%  Note that only conversions are possible, which can occure without loss of
%  precision and the conversion single->double.
%  Complex data can only be saved as 'single'.
%
%EXAMPLE
%   tom_emwritec3(vol);
%        a fileselect-box appears and the EM-file can be selected to
%        save vol to file. vol can be a EM-structure or a numerical array.
%
%   tom_emwritec3(filename, vol);
%       Writes the volume vol to EM-file filename.
%
%   tom_emwritec3(filename, vol, 'standard', 'single');
%       Writes the data 'vol' to the file 'filename' with type single.
%       The parameter standard is optional.
%
%   tom_emwritec3(filename, [1024 2048 512], 'new', 'double', header);
%       creates a new file of [1024 2048 512] dimensions in x, y and z in
%       datatype double. It initialises the EM-header with data from
%       header, where header is an EM-header or an EM-structure.
%       header and the datatype are optional.
%
%   tom_emwritec3(filename, vol, 'subregion', [10 20 30], reverse_sampling, allow_conversion);
%       Writes the data with an offset of [10 20 30] to the already existing
%       file 'filename'. The offset is zerobased, i.e. the first corner has
%       coordinate [0,0,0].
%       There are two optional parameters. reverse_sampling can be either
%       [] or a 3-vector specifying that the volume will not be pasted
%       contigously into the file but with gaps. It is exactly the opposite
%       of the sampling factor in TOM_EMREADC3. Using a sampling different
%       from 1 along X slows the process down significantly.
%       allow_conversion is a boolean flag, so that simple type conversions
%       of the volume are allowed. Even in that case only some conversions
%       are allowed. Only those where the conversion can happen without
%       loss of data and the conversions int32->single, double->single,
%       double->single_complex.
%       The EM-Header of the file remains unchanged.
%
%   tom_emwritec3(filename, vol, 'append', allow_conversion);
%       Appends the volume to the already existing EM-file filename, by
%       extending its Z-dimension. allow_conversion is a boolean flag
%       allowing some simple data-type conversions. The EM-Header of the
%       file remains unchanged (except a change of the volume size).
%
%REFERENCES
%
%SEE ALSO
%   TOM_EMREAD, TOM_MEX_EMWRITE, TOM_EMWRITE_PASTE
%
%   created by Thomas Haller Jan. 21 2007
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


if (isempty(varargin))
    error('No volume given to save as em.');
end;
if (~ischar(varargin{1}))
    filename = '';
else 
    filename = varargin{1};
    varargin = varargin(2:end);
end;
if (isempty(varargin))
    error('No volume given to save as em.');
end;


data = varargin{1};
varargin = varargin(2:end);

type = '';
if (isempty(varargin))
    command = 'standard';
else
    command = varargin{1};
    if (~ischar(command))
        error('The parameter command must be a string');
    end;
    if (~any(strcmp({'int8', 'int16', 'int32,' 'single', 'float', 'single_complex', 'float_complex', 'double'}, command)))
        varargin = varargin(2:end);
    else 
        command = 'standard';
    end;
end;

switch (lower(command))
    case 'new'
        if (~isnumeric(data) || ~any(numel(data)==[2,3]))
            error('For creating a new volume, the parameter vol must be the size of the em-file (i.e. a 3-vector with sizes [x,y,z])');
        end;
        if (numel(data) == 2)
            data(3) = 1;
        end;
        type = '';
        if (~isempty(varargin))
            if (ischar(varargin{1}))
                type = varargin{1};
                varargin = varargin(2:end);
            end;                
        end;
        old_emstruct = [];
        if (~isempty(varargin))
            if (~isstruct(varargin{1}) || length(varargin)~=1)
                error('Creating a new EM-file can have two optional parameters. First a string specifying the datatype, then a structure containing the header');
            end;
            old_emstruct = varargin{1};
            varargin = varargin(2:end);
        end;

        filename = get_filename(filename);
        if (isempty(filename))
            return;
        end;

        [magic, comment, emdata, userdata] = tom_emheader_struct_2_data(old_emstruct, type);
        tom_mex_emwrite_new(filename, data, magic, comment, emdata, userdata);
    case 'append'
        allow_conversion = false;
        if (~isempty(varargin))
            if (~(isnumeric(varargin{1})||islogical(varargin{1})) || all(numel(varargin{1})~=[0 1]))
                error('allow_conversion must be either [] or a boolean scalar .');
            end;
            if (~isempty(varargin{1}))
                allow_conversion = varargin{1};
            end;
        end;
        
        filename = get_filename(filename);
        if (isempty(filename))
            return;
        end;
        if (isnumeric(data))
            datan = data;
        elseif (isstruct(data) && numel(data)==1 && isfield(data,'Value') && isnumeric(data.Value))
            datan = data.Value;
        else
            error('Data must be a numeric array or an EM-structure');
        end;        
        %[size, magic, comment, emdata, userdata] = ...
        tom_mex_emwrite_append_stack(filename, datan, allow_conversion);
        
    case 'subregion'
        if (isempty(varargin) || ~isnumeric(varargin{1}) || numel(varargin{1})~=3)
            error('Subregion needs the first voxel where to begin pasting the volume.');
        end;
        subregion_start = varargin{1};
        reverse_sampling = [];
        allow_conversion = false;
        varargin = varargin(2:end);
        if (~isempty(varargin))
            reverse_sampling = varargin{1};
            if (~isnumeric(reverse_sampling) || (numel(reverse_sampling)~=3 && ~isempty(reverse_sampling)))
                error('Reverse sampling must be either [] or a 3-vector.');
            end;
            varargin = varargin(2:end);
        end;
        if (~isempty(varargin))
            if (~(isnumeric(varargin{1})||islogical(varargin{1})) || all(numel(varargin{1})~=[0 1]))
                error('allow_conversion must be either [] or a boolean scalar .');
            end;
            if (~isempty(varargin{1}))
                allow_conversion = varargin{1};
            end;
        end;
        
        filename = get_filename(filename);
        if (isempty(filename))
            return;
        end;
        if (isnumeric(data))
            datan = data;
        elseif (isstruct(data) && numel(data)==1 && isfield(data,'Value') && isnumeric(data.Value))
            datan = data.Value;
        else
            error('Data must be a numeric array or an EM-structure');
        end;
        
        %[size, magic, comment, emdata, userdata] = ...
        tom_mex_emwrite_paste(filename, datan, subregion_start, reverse_sampling, allow_conversion);            
        
    case 'standard'
        type = '';
        if (~isempty(varargin))
            if (ischar(varargin{1}))
                type = varargin{1};
                varargin = varargin(2:end);
            end;                
        end;
        if (~isempty(varargin))
            error('Working in standard mode has only one additional parameter: The data-type (as a string).');
        end;
        filename = get_filename(filename);
        if (isempty(filename))
            return;
        end;
        
        
        if (isnumeric(data))
            datan = data;
        elseif (isstruct(data) && numel(data)==1 && isfield(data,'Value') && isnumeric(data.Value))
            datan = data.Value;
        else
            error('Data must be a numeric array or an EM-structure');
        end;
        [magic, comment, emdata, userdata] = tom_emheader_struct_2_data(data, type);
        tom_mex_emwrite(filename, datan, [], magic, comment, emdata, userdata);
    otherwise
        error(['Unknown command ' command ]);
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filename = get_filename(filename)

if (isempty(filename))
    [filename, pathname] = uiputfile({'*.em';'*.vol';'*.*'}, 'Save as EM-file');
    if (~isequal(filename,0) && ~isequal(pathname,0))
        filename = [pathname filename];
    else
        filename = '';
    end;    
end;

