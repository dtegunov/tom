function [Data] = tom_emreadc_old(em_name,form,nr,nrarea)

% TOM_EMREADC reads data in EM-file format. C Version under development !
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
%    Syntax:
%     [Data] = tom_emreadc(em_name)
%    Input:
%    em_name:       Filename
%    Output:
%     Data:     	Structure of Image Data
%     Data.Value:   Raw data of image, or stack
%     Data.Header:  Header information
%
%    Example:
%               i=tom_emreadc;
%                   a fileselect-box appears and the EM-file can be picked
%
%               i=tom_emreadc('Proj.em');
%
%               i=tom_emreadc('HPIEMV','subregion',[102 162 1],[23 29 0]);
%                   reads a subregion starting from position (102, 162,1) of an image
%                   and gives back an image of size (23,29). Alternative to read
%                   whole volume and reduce by redi=i.Value(102:124,162:190);
%
%               i=tom_emreadc('TRIPODV','subregion',[14 17 19],[19 16 9]);
%                   reads a subregion starting from position (14,17,19) in the volume
%                   and gives back a volume of size (19,16,9). Alternative to read
%                   whole volume and reduce by redi=i.Value(14:33,17:33,19:28);
%
%               i=tom_emreadc('TRIPODV','resample',[1 2 3]);
%                   reads the file and resamples the file using every value
%                   in x direction, every second value in y direction and
%                   every third value in z direction
%
%    See Also
%     TOM_EMREAD, TOM_EMWRITE, TOM_EMHEADER, TOM_EMREADEMHEADER
%
%    09/23/02 SN
%    last change: 22/12/05 AK added resampling
%
%    Copyright (c) 2004
%    TOM toolbox for Electron Tomography
%    Max-Planck-Institute for Biochemistry
%    Dept. Molecular Structural Biology
%    82152 Martinsried, Germany
%    http://www.biochem.mpg.de/tom


emtype=cellstr(['extern         '; 'EM420          '; 'CM12           '; 'CM200          '; 'CM120/Biofilter'; 'CM300          '; 'Polara         ']);
%                                               EM420=1;CM12=2;CM200=3;
%                                               CM120/Biofilter=4;CM300=5;
%                                               Polara=6;extern=0;


error(nargchk(0,4,nargin))
if nargin <1 form='standard';  nr=1;
[filename, pathname] = uigetfile({'*.em';'*.vol';'*.norm';'*.mask';'*.*'}, 'Pick an EM-file');
if isequal(filename,0) | isequal(pathname,0) disp('No data loaded.'); return; end;
em_name=[pathname filename];
end;
if nargin <3 form='standard';  nr=1; nrarea=1; end;

if findstr(em_name,'/')
    filename=em_name(max(findstr(em_name,'/'))+1:size(em_name,2));
    pathname=em_name(1:max(findstr(em_name,'/')));
else
    filename=em_name;
    pathname=[pwd '/'];
    
end

if isequal(form,'standard')     
[Data.Value Data.Header.Magic Data.Header.Size Data.Header.Comment ...
    Data.Header.Parameter Data.Header.Fillup]=tom_emreadinc(em_name);
    Data.Header.EM=struct('Magic',Data.Header.Magic,'Size',Data.Header.Size,...
        'Comment',Data.Header.Comment,'Parameter',Data.Header.Parameter,'Fillup',Data.Header.Fillup);
    Data.Header.Voltage=Data.Header.Parameter(1);
    Data.Header.Cs=Data.Header.Parameter(2)./1000;
    Data.Header.Aperture=Data.Header.Parameter(3);
    Data.Header.Magnification=Data.Header.Parameter(4);
    Data.Header.Postmagnification=Data.Header.Parameter(5)./1000;
    Data.Header.Exposuretime=Data.Header.Parameter(6)./1000;
    Data.Header.Objectpixelsize=Data.Header.Parameter(7)./1000;
    if Data.Header.Parameter(8)<0 || Data.Header.Parameter(8)>6 Data.Header.Parameter(8)=0; end;
    Data.Header.Microscope=emtype(Data.Header.Parameter(8)+1);
    Data.Header.Pixelsize=Data.Header.Parameter(9)./1000;
    Data.Header.CCDArea=Data.Header.Parameter(10)./1000;
    Data.Header.Defocus=Data.Header.Parameter(11);
    Data.Header.Astigmatism=Data.Header.Parameter(12);
    Data.Header.AstigmatismAngle=Data.Header.Parameter(13)./1000;
    Data.Header.FocusIncrement=Data.Header.Parameter(14);
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

elseif isequal(form,'resample')
    if size(nr,1) ~= 1 | size(nr,2) ~= 3
        error('Binning must be a 3 element row vector');
    end
    if nr(1) < 1 | nr(2) < 1 | nr(3) < 1
        error('Binning factors must be positive integers');
    end
    [Data.Value Data.Header.Magic Data.Header.Size Data.Header.Comment ...
    Data.Header.Parameter Data.Header.Fillup]=tom_emreadinc_resample(em_name,nr);
    Data.Header.EM=struct('Magic',Data.Header.Magic,'Size',Data.Header.Size,...
        'Comment',Data.Header.Comment,'Parameter',Data.Header.Parameter,'Fillup',Data.Header.Fillup);
    Data.Header.Voltage=Data.Header.Parameter(1);
    Data.Header.Cs=Data.Header.Parameter(2)./1000;
    Data.Header.Aperture=Data.Header.Parameter(3);
    Data.Header.Magnification=Data.Header.Parameter(4);
    Data.Header.Postmagnification=Data.Header.Parameter(5)./1000;
    Data.Header.Exposuretime=Data.Header.Parameter(6)./1000;
    Data.Header.Objectpixelsize=Data.Header.Parameter(7)./1000;
    if Data.Header.Parameter(8)<0 || Data.Header.Parameter(8)>6 Data.Header.Parameter(8)=0; end;
    Data.Header.Microscope=emtype(Data.Header.Parameter(8)+1);
    Data.Header.Pixelsize=Data.Header.Parameter(9)./1000;
    Data.Header.CCDArea=Data.Header.Parameter(10)./1000;
    Data.Header.Defocus=Data.Header.Parameter(11);
    Data.Header.Astigmatism=Data.Header.Parameter(12);
    Data.Header.AstigmatismAngle=Data.Header.Parameter(13)./1000;
    Data.Header.FocusIncrement=Data.Header.Parameter(14);
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

else
    if nr(1) == 0 | nr(2) == 0 | nr(3) == 0
        error('Subregion start vector cannot be [0 0 0]');
    end
    
    
[Data.Value Data.Header.Magic Data.Header.Size Data.Header.Comment ...
    Data.Header.Parameter Data.Header.Fillup ]=tom_emreadinc(em_name,nr,nrarea);
    Data.Header.EM=struct('Magic',Data.Header.Magic,'Size',Data.Header.Size,...
        'Comment',Data.Header.Comment,'Parameter',Data.Header.Parameter,'Fillup',Data.Header.Fillup);
    Data.Header.Voltage=Data.Header.Parameter(1);
    Data.Header.Cs=Data.Header.Parameter(2)./1000;
    Data.Header.Aperture=Data.Header.Parameter(3);
    Data.Header.Magnification=Data.Header.Parameter(4);
    Data.Header.Postmagnification=Data.Header.Parameter(5)./1000;
    Data.Header.Exposuretime=Data.Header.Parameter(6)./1000;
    Data.Header.Objectpixelsize=Data.Header.Parameter(7);
    if Data.Header.Parameter(8)<0 || Data.Header.Parameter(8)>6 Data.Header.Parameter(8)=0; end;
    Data.Header.Microscope=emtype(Data.Header.Parameter(8)+1);
    Data.Header.Pixelsize=Data.Header.Parameter(9)./1000;
    Data.Header.CCDArea=Data.Header.Parameter(10)./1000;
    Data.Header.Defocus=Data.Header.Parameter(11);
    Data.Header.Astigmatism=Data.Header.Parameter(12);
    Data.Header.AstigmatismAngle=Data.Header.Parameter(13)./1000;
    Data.Header.FocusIncrement=Data.Header.Parameter(14);
    Data.Header.CountsPerElectron=Data.Header.Parameter(15)./1000;
    Data.Header.Intensity=Data.Header.Parameter(16)./1000;
    Data.Header.EnergySlitwidth=Data.Header.Parameter(17);
    Data.Header.EnergyOffset=Data.Header.Parameter(18);
    Data.Header.Tiltangle=Data.Header.Parameter(19)./1000;
    Data.Header.Tiltaxis=Data.Header.Parameter(20)./1000;
    Data.Header.Marker_X=Data.Header.Parameter(24);
    Data.Header.Marker_Y=Data.Header.Parameter(25);
    Data.Header.Username=num2str(Data.Header.Fillup(1:10));
    Data.Header.Date=num2str(Data.Header.Fillup(11:18));
    Data.Header.Filename=filename;
    Data.Header.Pathname=pathname;
    Data.Header.Size=size(Data.Value);if ndims(Data.Value)==2 Data.Header.Size(3)=1;end;
end;
