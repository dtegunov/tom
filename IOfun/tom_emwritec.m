function tom_emwritec(em_name,Data,form,nr,nrarea)
% TOM_EMWRITEC writes data in EM-file format.
%
%    Writes an EM-Image File (V-Format) 
%	 a raw format with a 512 Byte Header.
%    Keyword 'subregion' followed by a 3D vector 
%    and a subregion 3D vector writes only a subregion
%    to an existing file. The keyword 'new' creates an 
%    empty file of arbitrary size (also huge ones).
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
%                                       complex         8               8
%                                       double          8               9
%    -Three long integers (3x4 bytes) are image size in x, y, z Dimension
%    -80 Characters as comment
%    -40 long integers (4 x 40 bytes) are user defined parameters
%    -256 Byte with userdata, i.e. the username 
%    -Raw data following with the x variable as the fastest dimension, then y and z
%
%    -The parameters are coded as follwing:
%       No.  |  Name  |  Value  |  Factor  |  Comment
%       1       U        Volt      1000       accelerating voltage
%       2       COE      �m        1000       Cs of objective lense
%       3       APE      mrad      1000       aperture
%       4       VE       x         1          end magnification
%       5       VN       -         1000       postmagnification of CCD
%       6       ET       s         1000       exposure time in seconds
%       7       XG       -         1          pixelsize in object-plane
%       8       DG       nm        1000       EM-Code:
%                                               EM420=1;CM12=2;CM200=3;
%                                               CM120/Biofilter=4;CM300=5;
%                                               CM300/Tecnai=6;extern=0;
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
%       25      -        pixel     1          internal: subframe Y0header.Header.Magic[4] == 1
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
%    tom_emwritec(em_name)
%
%  INPUT
%   em_name             Filename
%   Data                matrix or file struct containing data to write
%                       optional format identifier: 'char' 'short' 'float' 'double'
%
%EXAMPLE
%   tom_emwritec(out);
%        a fileselect-box appears and the EM-file can be selected
%
%   tom_emwritec(out,Data);
%       writes the data 'Data' to the file 'out', format float;
%
%   tom_emwritec(out,Data,'standard','short');
%       writes the data 'Data' to the file 'out', format short;
%
%   tom_emwritec('Newvolume',[1024 2048 512],'new');
%       creates a new file of [1024 2048 512] dimensions in x, y and z.
%
%   tom_emwritec('Newvolume',[1024 2048 512],'new','double');
%       creates a new file of [1024 2048 512] dimensions in x, y and z, file format is double.
%
%   tom_emwritec('Volume',out,'subregion',[10 20 30],[63 63 15]);
%       writes the data with an offset of [10 20 30] to the file 'Volume'.
%       The last input argument is the size of the volume.
%
%REFERENCES
%
%SEE ALSO
%     TOM_EMREAD, TOM_EMREADC, TOM_EMWRITE, TOM_EMHEADER, TOM_EMREADEMHEADER
%
%   created by SN 12/09/02
%   updated by AK 08/08/06 added endianess swapping, support for files > 2G,
%    compatibily for tom_emwrite, support for data types other than float
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


error(nargchk(0,5,nargin))
if nargin <2 form='standard';  nr='float';
[filename, pathname] = uiputfile({'*.em';'*.*'}, 'Save as EM-file');
if isequal(filename,0) | isequal(pathname,0) disp('No data saved.'); return; end;
em_name=[pathname filename];
end;

if nargin == 2
    form='standard';
    if ~isstruct(Data)
        nr='float';
    else
        if Data.Header.EM.Magic(4) == 1
            nr = 'char';
        elseif Data.Header.EM.Magic(4) == 2
            nr = 'short';
        elseif Data.Header.EM.Magic(4) == 5
            nr = 'float';
        elseif Data.Header.EM.Magic(4) == 9
            nr = 'double';
        end
    end
end

if isequal(form,'standard') && nargin == 3
    nr = 'float';
end

if ~isequal(form,'subregion')
    %if nargin > 3
    try
        switch (nr)
            case 'char'
                nr = 1;
            case 'short'
                nr = 2;
            case 'float'
                nr = 4;
            case 'single'
                nr = 4;
            case 'double'
                nr = 8;
            otherwise
                error('format unknown');
        end
    catch
        nr=4;
    end
    %else
    %    nr = 4;
    %end
end

if isequal(form,'subregion')     
    
    if tom_isemfile(em_name) ~= 1
       error('File does not exist or is not readable');
    end
    
    header = tom_reademheader(em_name);
    %char
    if header.Header.Magic(4) == 1
        Data = int8(Data);
    %short    
    elseif header.Header.Magic(4) == 2
        Data = int16(Data);
    %float
    elseif header.Header.Magic(4) == 5
        Data = single(Data);
    %double
    elseif header.Header.Magic(4) == 9
        Data = double(Data);
    else
        error('Data format in file is not supported!');
    end
end



if isequal(form,'standard')
    
    if isstruct(Data)~=1
        Data=tom_emheader2(Data);
    end
    
    %char
    if nr == 1
        Data.Value = int8(Data.Value);
        Data.Header.EM.Magic(4) = 1;
    %short    
    elseif nr == 2
        Data.Value = int16(Data.Value);
        Data.Header.EM.Magic(4) = 2;
    %float
    elseif nr == 4
        Data.Value = single(Data.Value);
        Data.Header.EM.Magic(4) = 5;
    %double
    elseif nr == 8
        Data.Value = double(Data.Value);
        Data.Header.EM.Magic(4) = 9;
    else
        error('Data format in file is not supported!');
    end
    
    if length(Data.Header.EM.Size) == 2
        Data.Header.EM.Size(3) = 1;
    end
    
    Data.Header.EM.Magic(1) = 6;
    
end


if isequal(form,'standard')     

    Data.Header.EM.Parameter(1)=Data.Header.Voltage;
    Data.Header.EM.Parameter(2)=Data.Header.Cs.*1000;
    Data.Header.EM.Parameter(3)=Data.Header.Aperture;
    Data.Header.EM.Parameter(4)=Data.Header.Magnification;
    Data.Header.EM.Parameter(5)=Data.Header.Postmagnification.*1000; 
    Data.Header.EM.Parameter(6)=Data.Header.Exposuretime.*1000;
    Data.Header.EM.Parameter(7)=Data.Header.Objectpixelsize.*1000;
    Data.Header.EM.Parameter(9)=Data.Header.Pixelsize.*1000;
    Data.Header.EM.Parameter(10)=Data.Header.CCDArea.*1000;
    Data.Header.EM.Parameter(12)=Data.Header.Astigmatism;
    Data.Header.EM.Parameter(13)=Data.Header.AstigmatismAngle.*1000;
    Data.Header.EM.Parameter(14)=Data.Header.FocusIncrement;
    Data.Header.EM.Parameter(15)=Data.Header.CountsPerElectron.*1000;
    Data.Header.EM.Parameter(16)=Data.Header.Intensity.*1000;
    Data.Header.EM.Parameter(17)=Data.Header.EnergySlitwidth;
    Data.Header.EM.Parameter(18)=Data.Header.EnergyOffset;
    Data.Header.EM.Parameter(19)=Data.Header.Tiltangle.*1000;
    Data.Header.EM.Parameter(20)=Data.Header.Tiltaxis.*1000;
    Data.Header.EM.Parameter(24)=Data.Header.Marker_X;
    Data.Header.EM.Parameter(25)=Data.Header.Marker_Y;
    Data.Header.EM.Comment=Data.Header.Comment;
    
    
    Data.Header.Microscope = deblank(Data.Header.Microscope);
    if isequal(Data.Header.Microscope,'extern') Data.Header.EM.Parameter(8)=0; end
    if isequal(Data.Header.Microscope,'EM420') Data.Header.EM.Parameter(8)=1;end
    if isequal(Data.Header.Microscope,'CM12') Data.Header.EM.Parameter(8)=2; end;
    if isequal(Data.Header.Microscope,'CM200') Data.Header.EM.Parameter(8)=3; end;
    if isequal(Data.Header.Microscope,'CM120/Biofilter') Data.Header.EM.Parameter(8)=4; end;
    if isequal(Data.Header.Microscope,'CM300') Data.Header.EM.Parameter(8)=5; end;
    if isequal(Data.Header.Microscope,'Polara') Data.Header.EM.Parameter(8)=6; end;
    if isequal(Data.Header.Microscope,'Titan') Data.Header.EM.Parameter(8)=6; end;
    
    Data.Header.EM.Parameter(11)=Data.Header.Defocus;
    Data.Header.EM.Parameter(19)=Data.Header.Tiltangle.*1000;
    Data.Header.EM.Parameter(20)=Data.Header.Tiltaxis.*1000;

    if size(Data.Header.Comment,2) == 1 && size(Data.Header.Comment,1) > 1
        Data.Header.Comment = Data.Header.Comment';
    end
    
    if size(Data.Header.Comment,2)<80
        f=ones((80-size(Data.Header.Comment,2)),1).*32;
    	Data.Header.Comment(size(Data.Header.Comment,2)+1:80)=f;
    else
        Data.Header.Comment=Data.Header.Comment(1:80);
    end
    Data.Header.EM.Comment = int8(Data.Header.Comment);    
    %Data.Header.EM.Magic = int8([6 0 0 5]);
    tom_emwriteinc2(em_name,Data.Value,nr,Data.Header.EM.Magic,Data.Header.EM.Size,Data.Header.EM.Comment,Data.Header.EM.Parameter);
    
elseif isequal(form,'new')     
    if isstruct(Data)
        Data = Data.Value;
    end
    tom_emwriteinc2(em_name,Data,nr);
elseif isequal(form,'subregion')
    if isstruct(Data)
        Data = Data.Value;
    end
    tom_emwriteinc2(em_name,Data,nr,nrarea);
else
    error('unknown action in tom_emwritec.')
end;