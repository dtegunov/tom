function [Data] = tom_reademheader(em_name)
%TOM_READEMHEADER reads only the header from an EM-file
%
%   [Image.Header] = tom_reademheader(em_name)
%
%PARAMETERS
%
%  INPUT
%   em_name             filename
%  
%  OUTPUT
%   Image.Header        Header information
%
%EXAMPLE
%   in=tom_reademheader('proteasome.em')
%
%REFERENCES
%
%SEE ALSO
%   TOM_EMREAD, TOM_EMWRITE, TOM_EMHEADER, TOM_READEMHEADER
%
%   created by SN 08/07/02
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

if nargin ~=1 error(['Filename not specified (e.g. tom_reademheader(''c:\Data.em'')']);  end;

% open the stream with the correct format !
fid = fopen(em_name,'r','ieee-be');
if fid==-1
    error(['Cannot open: ' em_name ' file']); 
end;
magic = fread(fid,[4],'char');
fclose(fid);

% reads the header
%
% description in 'The Structure of the EM-Data Files', Herr Hegerl
% and at the bottom of this file

emtype=cellstr(['extern         '; 'EM420          '; 'CM12           '; 'CM200          '; 'CM120/Biofilter'; 'CM300          '; 'Polara         ']);
%                                               EM420=1;CM12=2;CM200=3;
%                                               CM120/Biofilter=4;CM300=5;
%                                               Polara=6;extern=0;

% read the Header 
if isempty(magic) 
    Data=[]; %It's not a em file. Data is empty
elseif (magic(1)==3 | magic(1)==0 | magic(1)==5)
    fid = fopen(em_name,'r','ieee-be'); % for SGI or OS-9 or Mac
else
    fid = fopen(em_name,'r','ieee-le'); % for PC
end;
if ~isempty(magic)
    magic = fread(fid,[4],'char');
    image_size = fread(fid,[3],'int32');
    comment = fread(fid,[80],'char');
    parameter = fread(fid,[40],'int32');
    fillup = fread(fid,[256],'char');
    
    EM=struct('Magic',magic,'Size',image_size,'Comment',comment,'Parameter',parameter,'Fillup',fillup);
    if parameter(8)>6 
        parameter(8)=0; 
    end;
    Header=struct(...
        'Voltage',parameter(1),...
        'Cs',parameter(2)./1000,...
        'Aperture',parameter(3),...
        'Magnification',parameter(4),...
        'Postmagnification',parameter(5)./1000,...
        'Exposuretime',parameter(6)./1000,...
        'Objectpixelsize',parameter(7)./1000,...
        'Microscope',emtype(parameter(8)+1),...
        'Pixelsize',parameter(9)./1000,...
        'CCDArea',parameter(10)./1000,...
        'Defocus',parameter(11),...
        'Astigmatism',parameter(12),...
        'AstigmatismAngle',parameter(13)./1000,...
        'FocusIncrement',parameter(14),...
        'CountsPerElectron',parameter(15)./1000,...
        'Intensity',parameter(16)./1000,...
        'EnergySlitwidth',parameter(17),...
        'EnergyOffset',parameter(18),...
        'Tiltangle',parameter(19)./1000,...
        'Tiltaxis',parameter(20)./1000,...
        'Marker_X',parameter(24),...
        'Marker_Y',parameter(25),...
        'Username',num2str(fillup(1:20)),...
        'Date',num2str(fillup(21:28)),...
        'Magic',magic,'Size',image_size,'Comment',comment','Parameter',parameter,'Fillup',fillup,'EM',EM);
    Data=struct('Header',Header);
end
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

