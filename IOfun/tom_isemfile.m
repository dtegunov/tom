function [out] = tom_isemfile(em_name)
% TOM_ISEMFILE tests for a file in EM-Format
%
%[out] = tom_isemfile(em_name)
%
%    Checks if the file is in EM-file format (V-Format),
%	 a raw format with a 512 Byte Header and wether the file
%    exists or not.
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
%       2       COE      ???m        1000       Cs of objective lense
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
%
%  INPUT
%   em_name             Filename
%
%  OUTPUT
%   out                 is 1 for true, yes the file is in EM-format.
%                       is 0 for false, no the file is not in EM-format.
%                       is -1 for file doesn't exist at all.
%
%SEE ALSO
%     TOM_EMWRITE, TOM_EMHEADER, TOM_EMREADEMHEADER
%
%   created by SN 10/21/02
%   updated by AK 08/11/06 added support for files > 2G, performance enhancement
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
    out=-1;
    return;
end;
magic = fread(fid,4,'char');
fclose(fid);

if (isempty(magic))
    out=-1;    
    return;
end;

% read the Header 
if (magic(1)==3 || magic(1)==0 || magic(1)==5)
    fid = fopen(em_name,'r','ieee-be'); % for SGI or OS-9 or Mac
elseif (magic(1)==1 || magic(1)==2 || magic(1)==4 || magic(1)==6)
    fid = fopen(em_name,'r','ieee-le'); % for PC
else 
    out=0; return;
end;    
magic = fread(fid,4,'char');
image_size = fread(fid,3,'int32');
fclose(fid);

switch magic(4)
    case 1
        varsize = 1;
    case 2
        varsize = 2;
    case 5
        varsize = 4;
    case 9
        varsize = 8;
    otherwise
        out = 0;
        
        return;
end

D = dir(em_name); 

if D(1).bytes ~= image_size(1).*image_size(2).*image_size(3).*varsize+512
    out=0;
else
    out=1;
end


