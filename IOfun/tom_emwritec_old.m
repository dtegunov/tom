function tom_emwritec_old(em_name,Data,form,nr,nrarea)
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
%   i=tom_emwritec(out);
%        a fileselect-box appears and the EM-file can be selected
%
%   i=tom_emwritec(out,Data);
%       writes the data 'Data' to the file 'out', format float;
%
%   i=tom_emwritec(out,Data,'standard','short');
%       writes the data 'Data' to the file 'out', format short;
%
%   i=tom_emwritec('Newvolume',[1024 2048 512],'new');
%       creates a new file of [1024 2048 512] dimensions in x, y and z.
%
%   i=tom_emwritec('Newvolume',[1024 2048 512],'new','double');
%       creates a new file of [1024 2048 512] dimensions in x, y and z, file format is double.
%
%   i=tom_emwritec('Volume',out,'subregion',[10 20 30],[63 63 15]);
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


emtype=cellstr(['extern         '; 'EM420          '; 'CM12           '; 'CM200          '; 'CM120/Biofilter'; 'CM300          '; 'Polara         ']);
%                                               EM420=1;CM12=2;CM200=3;
%                                               CM120/Biofilter=4;CM300=5;
%                                               Polara=6;extern=0;


error(nargchk(0,5,nargin))
if nargin <2 form='standard';  nr=1;
[filename, pathname] = uigetfile({'*.em';'*.*'}, 'Save as EM-file');
if isequal(filename,0) | isequal(pathname,0) disp('No data saved.'); return; end;
em_name=[pathname filename];
end;

if isequal(form,'standard')     
%  disp('not supported'); return;
  tom_emwriteinc(em_name,single(Data));
elseif isequal(form,'new')     
    tom_emwriteinc(em_name,Data);
else
    tom_emwriteinc(em_name,single(Data),nr,nrarea);
end;