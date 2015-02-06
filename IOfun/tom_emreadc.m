function Data = tom_emreadc(varargin)
%TOM_EMREADC reads data in EM-file format. 
%
%   Data = tom_emreadc(varargin)
%
%    Reads an EM-Image File (V-Format) 
%	 a raw format with a 512 Byte Header.
%    If no header was provided in EMWRITE
%    the header information will be
%    abandoned in EMREAD.
%    That way data can be saved 
%    and loaded without providing a header
%    but with compatible file-format to EM.
%    The keyword 'layer' followed by a number
%    reads only one layer from a 3D volume.
%    Keyword 'subregion' followed by a 3D vector 
%    and a subregion 3D vector reads only a subregion
%    from a 3D volume.
%
%    Structure of EM-Data Files:
%    -Byte 1: Machine Coding:       Machine:    Value:
%                                   OS-9         0
%                                   VAX          1
%                                   Convex       2
%                                   SGI          3
%                                   Mac          5
%                                   PC           6
%    -Byte 2: General purpose. On OS-9 system: 0 old version 1 is new version
%    -Byte 3: Not used in standard EM-format, if this byte is 1 the header is abandoned. 
%    -Byte 4: Data Type Coding:         Image Type:     No. of Bytes:   Value:
%                                       byte            1               1
%                                       short           2               2
%                                       long int        4               4
%                                       float           4               5
%                                       float complex   8               8
%                                       double          8               9
%                                       double complex  16              10
%    -Three long integers (3x4 bytes) are image size in x, y, z Dimension
%    -80 Characters as comment
%    -40 long integers (4 x 40 bytes) are user defined parameters
%    -256 Byte with userdata, i.e. the username 
%    -Raw data following with the x variable as the fastest dimension, then y and z
%
%    -The parameters are coded as follwing:
%       No.  |  Name  |  Value  |  Factor  |  Comment
%       1       U        Volt      1000       accelerating voltage
%       2       COE      ???m        1000       Cs of objective lense
%       3       APE      mrad      1000       aperture
%       4       VE       x         1          end magnification
%       5       VN       -         1000       postmagnification of CCD
%       6       ET       s         1000       exposure time in seconds
%       7       XG       -         1          pixelsize in object-plane
%       8       DG       nm        1000       EM-Code:
%                                               EM420=1;CM12=2;CM200=3;
%                                               CM120/Biofilter=4;CM300=5;
%                                               Polara=6;Titan=7;extern=0;
%       9       APD      nm        1000       photometer aperture
%       10      L        nm        1000       phys_pixel_size * nr_of_pixels
%       11      DF       Angstr.   1          defocus, underfocus is neg.
%       12      FA       Angstr.   1          astigmatism
%       13      PHI      deg/1000  1000       angle of astigmatism 
%       14      DR       Angstr.   1          drift in Angstr.
%       15      DELT     deg/1000  1000       direction of drift
%       16      DDF      Angstr.   1          focusincr. for focus-series
%       17      X0       -         1          obsolete
%       18      Y0       -         1          obsolete
%       19      KW       deg/1000  1000       tiltangle 
%       20      KR       deg/1000  1000       axis perpend. to tiltaxis
%       21      -        Angstr.   1           
%       22      SC       ASCII     1
%       23      -        -         -
%       24      -        pixel     1          internal: subframe X0
%       25      -        pixel     1          internal: subframe Y0
%       26      -        Angstr.   1000       internal: resolution
%       27      -        -         -          internal: density
%       28      -        -         -          internal: contrast
%       29      -        -         -          internal: unknown
%       30      SP       -         1000       mass centre X
%       31      SP       -         1000       mass centre Y
%       32      SP       -         1000       mass centre Z
%       33      H        -         1000       height
%       34      -        -         1000       internal: unknown
%       35      D1       -         1000       width 'Dreistrahlbereich'
%       36      D2       -         1000       width 'Achrom. Ring'
%       37      -        -         1          internal: lambda
%       38      -        -         1          internal: delta theta
%       39      -        -         1          internal: unknown
%       40      -        -         1          internal: unknown
%
%  INPUT
%   em_name             filename
%    
%   one or more of the following parameters:
%   'subregion',[startx starty startz],[sizex sizey sizez]
%   'resample',[x y z]
%   'binning',nr_binning (can be a scalar or a [x y z] vector for
%                          different binning in each direction)
%
%    
%   OUTPUT:
%    Data               Structure of Image Data
%    Data.Value         Raw data of image, or stack
%    Data.Header        Header information
%
%EXAMPLE
%    i=tom_emreadc;
%       a fileselect-box appears and the EM-file can be picked
%
%    i=tom_emreadc('Proj.em');
%
%    i=tom_emreadc('HPIEMV','subregion',[102 162 1],[23 29 0]);
%       reads a subregion starting from position (102, 162,1) of an image
%       and gives back an image of size (23,29). Alternative to read
%       whole volume and reduce by redi=i.Value(102:124,162:190);
%
%    i=tom_emreadc('TRIPODV','subregion',[14 17 19],[19 16 9]);
%       reads a subregion starting from position (14,17,19) in the volume
%       and gives back a volume of size (20,17,10). Alternative to read
%       whole volume and reduce by redi=i.Value(14:33,17:33,19:28);
%
%    i=tom_emreadc('TRIPODV','resample',[1 2 3]);
%       reads the file and resamples the file using every value
%       in x direction, every second value in y direction and
%       every third value in z direction
%    i=tom_emreadc('TRIPODV','binning',2);
%       reads the file and bins the file 2 times
%
%REFERENCES
%
%SEE ASLO
%   TOM_EMREAD, TOM_EMWRITE, TOM_EMHEADER, TOM_EMREADEMHEADER
%
%   created by SN 09/23/02
%   updated by AK 12/22/05 added resampling
%   updated by AK 02/05/06 added completely new c backend with support for subregion
%       and resampling at the same time
%   updated by AK 03/08/06 added on-disk binning and byte swapping and support for
%       files > 2G on 64 bit Linux systems
%   updated by AK 08/31/06 fixed endian swap when binning is used, added binning
%       values for each dimension
%   updated by SN 24/11/06 new microscope code added
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

emtype=cellstr(['extern         '; 'EM420          '; 'CM12           '; 'CM200          '; 'CM120/Biofilter'; 'CM300          '; 'Polara         '; 'Titan          '; 'Tecnai F20     ']);
%                                               EM420=1;CM12=2;CM200=3;
%                                               CM120/Biofilter=4;CM300=5;
%                                               Polara=6;Titan=7;extern=0;


error(nargchk(0,8,nargin))

[em_name, resample, binning, subregion_start, subregion_size] = ParseInputs(varargin{:});

if nargin <1 
    [filename, pathname] = uigetfile({'*.*';'*.em';'*.vol';'*.norm';'*.mask'}, 'Pick an EM-file');
    if isequal(filename,0) || isequal(pathname,0); disp('No data loaded.'); return; end;
    em_name=[pathname filename];
end;

if findstr(em_name,'/')
    filename=em_name(max(findstr(em_name,'/'))+1:size(em_name,2));
    pathname=em_name(1:max(findstr(em_name,'/')));
else
    filename=em_name;
    pathname=[pwd '/'];
    
end

if size(binning,1) == 1 && size(binning,2) == 1
    binning = [binning binning binning];
end

if binning == 0
    binning = [];
else
    binning = 2.^binning;
end

if ~isempty(subregion_start)
    if subregion_start(1) <= 0 || subregion_start(2) <= 0 || subregion_start(3) <= 0
        error('subregion start must be > 0');
    end
end

if isempty(resample) && isempty(binning) && isempty(subregion_start)   
    [Data.Value Data.Header.Magic Data.Header.Size Data.Header.Comment, Data.Header.Parameter Data.Header.Fillup]=tom_emreadinc_resample2(em_name);
elseif isempty(resample) && isempty(binning) && ~isempty(subregion_start)
    [Data.Value Data.Header.Magic Data.Header.Size Data.Header.Comment, Data.Header.Parameter Data.Header.Fillup]=tom_emreadinc_resample2(em_name,[1 1 1],subregion_start,subregion_size+1);
elseif isempty(subregion_start) && ~isempty(resample) && isempty(binning)
    [Data.Value Data.Header.Magic Data.Header.Size Data.Header.Comment, Data.Header.Parameter Data.Header.Fillup]=tom_emreadinc_resample2(em_name,resample);
elseif isempty(subregion_start) && ~isempty(resample) && ~isempty(binning)
    [Data.Value Data.Header.Magic Data.Header.Size Data.Header.Comment, Data.Header.Parameter Data.Header.Fillup]=tom_emreadinc_resample3(em_name,resample,binning);
elseif isempty(subregion_start) && isempty(resample) && ~isempty(binning)
    [Data.Value Data.Header.Magic Data.Header.Size Data.Header.Comment, Data.Header.Parameter Data.Header.Fillup]=tom_emreadinc_resample3(em_name,[1 1 1],binning);
elseif ~isempty(subregion_start) && isempty(resample) && ~isempty(binning)
    [Data.Value Data.Header.Magic Data.Header.Size Data.Header.Comment, Data.Header.Parameter Data.Header.Fillup]=tom_emreadinc_resample3(em_name,[1 1 1],binning,subregion_start,subregion_size+1);
else 
    [Data.Value Data.Header.Magic Data.Header.Size Data.Header.Comment, Data.Header.Parameter Data.Header.Fillup]=tom_emreadinc_resample3(em_name,resample,binning,subregion_start,subregion_size+1);
end


Data.Header.EM=struct('Magic',Data.Header.Magic,'Size',Data.Header.Size,'Comment',Data.Header.Comment,'Parameter',Data.Header.Parameter,'Fillup',Data.Header.Fillup);
Data.Header.Voltage=Data.Header.Parameter(1);
Data.Header.Cs=Data.Header.Parameter(2)./1000;
Data.Header.Aperture=Data.Header.Parameter(3);
Data.Header.Magnification=Data.Header.Parameter(4);
Data.Header.Postmagnification=Data.Header.Parameter(5)./1000;
Data.Header.Exposuretime=Data.Header.Parameter(6)./1000;
Data.Header.Objectpixelsize=Data.Header.Parameter(7)./1000;
if Data.Header.Parameter(8)<0 || Data.Header.Parameter(8)>8; Data.Header.Parameter(8)=0; end;
Data.Header.Microscope=emtype{Data.Header.Parameter(8)+1};
Data.Header.Pixelsize=Data.Header.Parameter(9)./1000;
Data.Header.CCDArea=Data.Header.Parameter(10)./1000;
Data.Header.Defocus=Data.Header.Parameter(11);
Data.Header.Astigmatism=Data.Header.Parameter(12);
Data.Header.AstigmatismAngle=Data.Header.Parameter(13)./1000;
Data.Header.FocusIncrement=Data.Header.Parameter(14); % changed by SN 060906
Data.Header.CountsPerElectron=Data.Header.Parameter(15)./1000;
Data.Header.Intensity=Data.Header.Parameter(16)./1000;
Data.Header.EnergySlitwidth=Data.Header.Parameter(17);
Data.Header.EnergyOffset=Data.Header.Parameter(18);
Data.Header.Tiltangle=Data.Header.Parameter(19)./1000;
Data.Header.Tiltaxis=Data.Header.Parameter(20)./1000;
Data.Header.Marker_X=Data.Header.Parameter(24);
Data.Header.Marker_Y=Data.Header.Parameter(25);
Data.Header.Username=num2str(Data.Header.Fillup(1:20));
Data.Header.Date=num2str(Data.Header.Fillup(11:18));
Data.Header.Filename=filename;
Data.Header.Pathname=pathname;
Data.Header.Size=size(Data.Value);if ndims(Data.Value)==2; Data.Header.Size(3)=1;end;

if ~isempty(binning)
    Data.Header.Objectpixelsize = Data.Header.Objectpixelsize .* binning(1);
end


%%%
%%% Subfunction ParseInputs
%%%
function [filename, resample, binning, subregion_start, subregion_size] = ParseInputs(varargin)

filename = '';
resample = [];
binning = 0;
subregion_start = [];
subregion_size = [];

switch nargin

case 0
    

case 1
    filename = varargin{1};
    

case 3
    filename = varargin{1};
    if strcmp(varargin{2}, 'resample')
        if size(varargin{3},2) == 3 && size(varargin{3},1) == 1
            resample = varargin{3};
        else
            error('Resampling value must be a 3 element vector');
        end
    elseif strcmp(varargin{2}, 'binning')
        if isnumeric(varargin{3}) && size(varargin{3},1) == 1 && size(varargin{3},2) == 1
            binning = varargin{3};
        elseif isnumeric(varargin{3}) && length(varargin{3}) == 3
            binning = varargin{3};
        else
            error('Binning value must be a scalar');
        end
    else
        error('Illegal arguments.');
    end
    

case 4
    filename = varargin{1};
    if strcmp(varargin{2},'subregion')
        if size(varargin{3},2) == 3 && size(varargin{3},1) == 1 && size(varargin{4},2) == 3 && size(varargin{4},1) == 1
            subregion_start = varargin{3};
            subregion_size = varargin{4};
        else
            error('Subregion definition must be two 3 element vectors');
        end
    end

case 5
    filename = varargin{1};
    if strcmp(varargin{2}, 'resample')
        if size(varargin{3},2) == 3 && size(varargin{3},1) == 1
            resample = varargin{3};
        else
            error('Resampling value must be a 3 element vector');
        end
    elseif strcmp(varargin{2}, 'binning')
        if isnumeric(varargin{3}) && size(varargin{3},1) == 1 && size(varargin{3},2) == 1
            binning = varargin{3};
        elseif isnumeric(varargin{3}) && length(varargin{3}) == 3
            binning = varargin{3};
        else
            error('Binning value must be a scalar');
        end
    else
        error('Illegal arguments.');
    end

    if strcmp(varargin{4}, 'resample')
        if size(varargin{5},2) == 3 && size(varargin{5},1) == 1
            resample = varargin{5};
        else
            error('Resampling value must be a 3 element vector');
        end
    elseif strcmp(varargin{4}, 'binning')
        if isnumeric(varargin{5}) && size(varargin{5},1) == 1 && size(varargin{5},2) == 1
            binning = varargin{5};
        elseif isnumeric(varargin{5}) && length(varargin{5}) == 3
            binning = varargin{5};
        else
            error('Binning value must be a scalar');
        end
    else
        error('Illegal arguments.');
    end
    

case 6
    filename = varargin{1};
    if strcmp(varargin{2}, 'resample')
        if size(varargin{3},2) == 3 && size(varargin{3},1) == 1
            resample = varargin{3};
        else
            error('Resampling value must be a 3 element vector');
        end
        if strcmp(varargin{4}, 'subregion')
            if size(varargin{5},2) == 3 && size(varargin{5},1) == 1 && size(varargin{6},2) == 3 && size(varargin{6},1) == 1
                subregion_start = varargin{5};
                subregion_size = varargin{6};
            else
                error('Subregion definition must be two 3 element vectors');
            end
        else
            error('Illegal arguments.');
        end
    elseif strcmp(varargin{2}, 'subregion')
        if size(varargin{3},2) == 3 && size(varargin{3},1) == 1 && size(varargin{4},2) == 3 && size(varargin{4},1) == 1
            subregion_start = varargin{3};
            subregion_size = varargin{4};
        else
            error('Subregion definition must be two 3 element vectors');
        end
        if strcmp(varargin{5},'resample')
            if size(varargin{6},2) == 3 && size(varargin{6},1) == 1
                resample = varargin{6};
            else
                error('Resampling value must be a 3 element vector');
            end
        elseif strcmp(varargin{5}, 'binning')
            if isnumeric(varargin{6}) && size(varargin{6},1) == 1 && size(varargin{6},2) == 1
                binning = varargin{6};
            elseif isnumeric(varargin{6}) && length(varargin{6}) == 3
                binning = varargin{6};
            else
                error('Binning value must be a scalar');
            end
        else
            error('Illegal arguments.');
        end
    elseif strcmp(varargin{2}, 'binning')
        if isnumeric(varargin{3}) && size(varargin{3},1) == 1 && size(varargin{3},2) == 1
            binning = varargin{3};
        elseif isnumeric(varargin{3}) && length(varargin{3}) == 3
            binning = varargin{3};
        else
            error('Binning value must be a scalar');
        end
        if strcmp(varargin{4}, 'subregion')
            if size(varargin{5},2) == 3 && size(varargin{5},1) == 1 && size(varargin{6},2) == 3 && size(varargin{6},1) == 1
                subregion_start = varargin{5};
                subregion_size = varargin{6};
            else
                error('Subregion definition must be two 3 element vectors');
            end
        else
            error('Illegal arguments.');
        end

    else
        error('Illegal arguments.');
    end
    

case 8
    filename = varargin{1};
    if strcmp(varargin{2}, 'resample')
        if size(varargin{3},2) == 3 && size(varargin{3},1) == 1
            resample = varargin{3};
        else
            error('Resampling value must be a 3 element vector');
        end
        if strcmp(varargin{4}, 'subregion')
            if size(varargin{5},2) == 3 && size(varargin{5},1) == 1 && size(varargin{6},2) == 3 && size(varargin{6},1) == 1
                subregion_start = varargin{5};
                subregion_size = varargin{6};
            else
                error('Subregion definition must be two 3 element vectors');
            end
            if strcmp(varargin{7}, 'binning')
                if isnumeric(varargin{8}) && size(varargin{8},1) == 1 && size(varargin{8},2) == 1
                    binning = varargin{8};
                elseif isnumeric(varargin{8}) && length(varargin{8}) == 3
                    binning = varargin{8};
                else
                    error('Binning value must be a scalar');
                end
            else
                error('Illegal arguments.');
            end

        elseif strcmp(varargin{4},'binning')
            if isnumeric(varargin{5}) && size(varargin{5},1) == 1 && size(varargin{5},2) == 1
                binning = varargin{5};
            elseif isnumeric(varargin{5}) && length(varargin{5}) == 3
               binning = varargin{5};
            else
                error('Binning value must be a scalar');
            end

            if strcmp(varargin{6}, 'subregion')
                if size(varargin{7},2) == 3 && size(varargin{7},1) == 1 && size(varargin{8},2) == 3 && size(varargin{8},1) == 1
                    subregion_start = varargin{7};
                    subregion_size = varargin{8};
                else
                    error('Subregion definition must be two 3 element vectors');
                end
            else
                error('Illegal arguments.');
            end

        else
            error('Illegal arguments.');
        end



    elseif strcmp(varargin{2}, 'binning')
        if isnumeric(varargin{3}) && size(varargin{3},1) == 1 && size(varargin{3},2) == 1
            binning = varargin{3};
        elseif isnumeric(varargin{3}) && length(varargin{3}) == 3
            binning = varargin{3};
        else
            error('Binning value must be a scalar');
        end
        if strcmp(varargin{4}, 'subregion')
            if size(varargin{5},2) == 3 && size(varargin{5},1) == 1 && size(varargin{6},2) == 3 && size(varargin{6},1) == 1
                subregion_start = varargin{5};
                subregion_size = varargin{6};
            else
                error('Subregion definition must be two 3 element vectors');
            end
            if strcmp(varargin{7},'resample')
                if size(varargin{8},2) == 3 && size(varargin{8},1) == 1
                    resample = varargin{8};
                else
                    error('Resampling value must be a 3 element vector');
                end
            else
                error('Illegal arguments.');
            end
        elseif strcmp(varargin{4}, 'resample')
            if size(varargin{5},2) == 3 && size(varargin{5},1) == 1
                resample = varargin{5};
            else
                error('Resampling value must be a 3 element vector');
            end
            if strcmp(varargin{6},'subregion')
                if size(varargin{7},2) == 3 && size(varargin{7},1) == 1 && size(varargin{8},2) == 3 && size(varargin{8},1) == 1
                    subregion_start = varargin{7};
                    subregion_size = varargin{8};
                else
                    error('Subregion definition must be two 3 element vectors');
                end
            end
        else
            error('Illegal arguments.');
        end



    elseif strcmp(varargin{2},'subregion')
        if size(varargin{3},2) == 3 && size(varargin{3},1) == 1 && size(varargin{4},2) == 3 && size(varargin{4},1) == 1
            subregion_start = varargin{3};
            subregion_size = varargin{4};
        else
            error('Subregion definition must be two 3 element vectors');
        end
        if strcmp(varargin{5},'resample')
            if size(varargin{6},2) == 3 && size(varargin{6},1) == 1
                resample = varargin{6};
            else
                error('Resampling value must be a 3 element vector');
            end

            if strcmp(varargin{7},'binning')
                if isnumeric(varargin{8}) && size(varargin{8},1) == 1 && size(varargin{8},2) == 1
                    binning = varargin{8};
                elseif isnumeric(varargin{8}) && length(varargin{8}) == 3
                    binning = varargin{8};
                else
                    error('Binning value must be a scalar');
                end
            else
                error('Illegal arguments.');
            end

        elseif strcmp(varargin{5},'binning')
            if isnumeric(varargin{6}) && size(varargin{6},1) == 1 && size(varargin{6},2) == 1
                binning = varargin{6};
            elseif isnumeric(varargin{6}) && length(varargin{6}) == 3
                binning = varargin{6};
            else
                error('Binning value must be a scalar');
            end
            if strcmp(varargin{7},'resample')
                if size(varargin{8},2) == 3 && size(varargin{8},1) == 1
                    resample = varargin{8};
                else
                    error('Resampling value must be a 3 element vector');
                end
            else
                error('Illegal arguments.');
            end
        else
            error('Illegal arguments.');
        end

    else
        error('Illegal arguments.');
    end
    

otherwise
    error('Illegal arguments.');
end

