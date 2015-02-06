function tom_writeemheader(em_name, header_in)
%TOM_WRITEEMHEADER writes only the header to an EM-file
%
%   tom_writeemheader(em_name,header_in)
%
%PARAMETERS
%
%  INPUT
%   em_name             name of the image
%   header_in           header structure
%  
%  OUTPUT
%
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_EMREAD, TOM_EMWRITE, TOM_EMHEADER, TOM_READEMHEADER
%
%   created by AK 08/09/06
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


% open the stream with the correct format !
fid = fopen(em_name,'r','ieee-be');
if fid==-1
    error(['Cannot open: ' em_name ' file']); 
end;
magic = fread(fid,4,'char');
fclose(fid);

% writes the header
%
% description in 'The Structure of the EM-Data Files', Herr Hegerl
% and at the bottom of this file

emtype=cellstr(['extern         '; 'EM420          '; 'CM12           '; 'CM200          '; 'CM120/Biofilter'; 'CM300          '; 'Polara         '; 'Titan          ']);
%                                               EM420=1;CM12=2;CM200=3;
%                                               CM120/Biofilter=4;CM300=5;
%                                               Polara=6;extern=0;


% set parameter from header structure

Data.Header.EM.Parameter=zeros(40,1);


Data.Header.EM.Parameter(1)=header_in.Voltage;
Data.Header.EM.Parameter(2)=header_in.Cs.*1000;
Data.Header.EM.Parameter(3)=header_in.Aperture;
Data.Header.EM.Parameter(4)=header_in.Magnification;
Data.Header.EM.Parameter(5)=header_in.Postmagnification.*1000; 
Data.Header.EM.Parameter(6)=header_in.Exposuretime.*1000;
Data.Header.EM.Parameter(7)=header_in.Objectpixelsize.*1000;
Data.Header.EM.Parameter(9)=header_in.Pixelsize.*1000;
Data.Header.EM.Parameter(10)=header_in.CCDArea.*1000;
Data.Header.EM.Parameter(12)=header_in.Astigmatism;
Data.Header.EM.Parameter(13)=header_in.AstigmatismAngle.*1000;
Data.Header.EM.Parameter(14)=header_in.FocusIncrement;
Data.Header.EM.Parameter(15)=header_in.CountsPerElectron.*1000;
Data.Header.EM.Parameter(16)=header_in.Intensity.*1000;
Data.Header.EM.Parameter(17)=header_in.EnergySlitwidth;
Data.Header.EM.Parameter(18)=header_in.EnergyOffset;

Data.Header.EM.Parameter(19)=header_in.Tiltangle.*1000;
Data.Header.EM.Parameter(20)=header_in.Tiltaxis.*1000;
Data.Header.EM.Parameter(24)=header_in.Marker_X;
Data.Header.EM.Parameter(25)=header_in.Marker_Y;

%     'Voltage',parameter(1),...
%     'Cs',parameter(2)./1000,...
%     'Aperture',parameter(3),...
%     'Magnification',parameter(4),...
%     'Exposuretime',parameter(6)./1000,...
%     'Objectpixelsize',parameter(7)./1000,...
%     'Microscope',emtype(parameter(8)+1),...
%     'Pixelsize',parameter(9)./1000,...
%     'CCDArea',parameter(10)./1000,...
%     'Defocus',parameter(11),...
%     'Astigmatism',parameter(12),...
%     'AstigmatismAngle',parameter(13)./1000,...
%     'FocusIncrement',parameter(14)./1,...
%     'CountsPerElectron',parameter(15)./1000,...
%     'Intensity',parameter(16)./1000,...
%     'EnergySlitwidth',parameter(17),...
%     'EnergyOffset',parameter(18),...
%     'Tiltangle',parameter(19)./1000,...
%     'Tiltaxis',parameter(20)./1000,...
%     'Username',num2str(fillup(1:20)),...
%     'Date',num2str(fillup(21:28)),...
%     'Magic',magic,'Size',image_size,'Comment',comment,'Parameter',parameter,'Fillup',fillup,'EM',EM);

try 
    Data.Header.Microscope = deblank(Data.Header.Microscope);
end
    

if isequal(header_in.Microscope,'extern') Data.Header.EM.Parameter(8)=0; end;
if isequal(header_in.Microscope,'EM420') Data.Header.EM.Parameter(8)=1; end;
if isequal(header_in.Microscope,'CM12') Data.Header.EM.Parameter(8)=2; end;
if isequal(header_in.Microscope,'CM200') Data.Header.EM.Parameter(8)=3; end;
if isequal(header_in.Microscope,'CM120/Biofilter') Data.Header.EM.Parameter(8)=4; end;
if isequal(header_in.Microscope,'CM300') Data.Header.EM.Parameter(8)=5; end;
if isequal(header_in.Microscope,'Polara') Data.Header.EM.Parameter(8)=6; end;
if isequal(header_in.Microscope,'Titan') Data.Header.EM.Parameter(8)=6; end;

Data.Header.EM.Parameter(11)=header_in.Defocus;
Data.Header.EM.Parameter(19)=header_in.Tiltangle.*1000;
Data.Header.EM.Parameter(20)=header_in.Tiltaxis.*1000;
% Comment
%    'Username',num2str(fillup(1:10)),...
%    'Date',num2str(fillup(11:18)),...


if isequal (computer,'PCWIN') || isequal(computer,'GLNX86') || isequal(computer,'GLNXA64') || isequal(computer,'MACI')
    fid = fopen(em_name,'r+','ieee-le'); Data.Header.Magic(1)=6;
else
    fid = fopen(em_name,'r+','ieee-be');
end;
if fid==-1 
    error(['Cannot write: ' em_name ' file']); 
end
% writes the header
%
% description in 'The Structure of the EM-Data Files', Herr Hegerl
%

fseek(fid, 0,-1);


if size(header_in.Magic,1)~=4
    header_in.Magic=[6 0 0 5]';
end
fwrite(fid,header_in.Magic,'char');

if length(header_in.Size) == 2
    header_in.Size = [header_in.Size;1];
end

fwrite(fid,header_in.Size(1:3),'int32');
if size(header_in.Comment,2)<80
    f=ones((80-size(header_in.Comment,2)),1).*32;
    header_in.Comment(size(header_in.Comment,2)+1:80)=f;
else
    header_in.Comment=header_in.Comment(1:80);
end
fwrite(fid,header_in.Comment,'char');

if size(Data.Header.EM.Parameter,1)~=40
    Data.Header.EM.Parameter=ones(40,1);
end
fwrite(fid,Data.Header.EM.Parameter,'int32');

if size(header_in.Fillup,1)~=256
    header_in.Fillup=ones(256,1);
end
fwrite(fid,header_in.Fillup,'char');


fclose(fid);





% TOM_EMREAD reads data in EM-file format
%
%    Reads an EM-Image File (V-Format) 
%	 a raw format with a 512 Byte Header.
%    If no header was provided in EMWRITE
%    the header information will be
%    abandoned in EMREAD.
%    That way data can be saved 
%    and loaded without providing a header
%    but with compatible file-format to EM.
%    The keyword 'subregion' followed by a 3D vector 
%    and a subregion 3D vector reads only a subregion
%    from a 3D volume. Only supported for float values in 3D
%    ; int2 and char are supported only for 2D images.
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
%    -256 Byte with userdata, first 20 chars username, 8 chars date (i.e.03/02/03)
%    -Raw data following with the x variable as the fastest dimension, then y and z
%
%    -The parameters are coded as follwing:
%       No.  |  Name  |  Value  |  Factor  |  Comment
%       1       U        Volt      1000       accelerating voltage
%       2       COE      ???m        1000       Cs of objective lense
%       3       APE      mrad      1000       aperture
%       4       VE       x         1          end magnification
%       5       VN       1000      1000       postmagnification of CCD (fixed value:1000!)
%       6       ET       s         1000       exposure time in seconds
%       7       OBJ      nm        1000       pixelsize in object-plane
%       8       EM       nm        1000       EM-Code:
%                                               EM420=1;CM12=2;CM200=3;
%                                               CM120/Biofilter=4;CM300=5;
%                                               Polara=6;extern=0;
%       9       CCD      nm        1000       physical pixelsize on CCD
%       10      L        nm        1000       phys_pixel_size * nr_of_pixels
%       11      DF       Angstr.   1          defocus, underfocus is neg.
%       12      FA       Angstr.   1          astigmatism
%       13      PHI      deg       1000       angle of astigmatism 
%       14      DDF      Angstr.   1          focusincr. for focus-series
%       15      CTS      -         1000       counts per primary electron, sensitivity of CCD
%       16      C2       -         1000       intensity value of C2
%       17      EW       eV        1          0 for no slit, x>0 for positive slitwidth 
%       18      EO       eV        1          energy offset from zero-loss
%       19      KW       deg       1000       tiltangle 
%       20      KR       deg       1000       axis perpend. to tiltaxis
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

