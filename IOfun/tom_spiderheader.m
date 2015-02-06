function Data = tom_spiderheader(in)
%TOM_SPIDERHEADER adds a  Spider header structure to a matrix
%
%   Data = tom_spiderheader(in)
%
%	Build a  Spider structure and adds a 
%   header to the in-Values with default values
%
%PARAMETERS
%
%  INPUT
%   in                  (Matrix or Volume)
%  
%  OUTPUT
%   out                 Structure in  spider format with header and in.Value
%   out.Value           Raw data of in
%   out.Header          Header information with standard values
%
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_SPIDERWRITE, TOM_SPIDERREAD, TOM_ISSPIDERFILE
%
%   created by AK 04/26/06
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

if isstruct(in)~=1
    Data = struct();
    Data.Header = struct();
    Data.Header.Spider = struct();
    Data.Value = single(in);

    %number of slices in volume (size of volume in z direction)
    Data.Header.Spider.nslice = size(in,3);

    %number of rows per slice (size of volume in x direction)
    Data.Header.Spider.nrow = size(in,1);

    %total number of records in the file (unused)
    Data.Header.Spider.irec = 0;

    %(obsolete, unused)
    Data.Header.Spider.nhistrec = 0;
    
    %file type specifier. Obsolete file types d, 8, 11, 12, 16, -1, -3, -7, and -9 are no longer supported in SPIDER.
    %iform  	(type)  	data type
    %1 	(r) 	2D image.
    %3 	(r) 	3D volume.
    %-11 	(fo) 	2D Fourier, mixed radix odd.
    %-12 	(fe) 	2D Fourier, mixed radix even.
    %-21 	(fo) 	3D Fourier, mixed radix odd.
    %-22 	(fe) 	3D Fourier, mixed radix even.
    if size(in,3) == 1
        Data.Header.Spider.iform = 1;
    else
        Data.Header.Spider.iform = 3;
    end
    
    %imami = maximum/minimum flag. Is set at 0 when the file is created, and at 1 when the maximum, minimum, average, and
    %standard deviation have been computed and stored into this header record (see following locations).
    Data.Header.Spider.imami = 0;
    
    %maximum value
    Data.Header.Spider.fmax = 0;

    %minimum value
    Data.Header.Spider.fmin = 0;

    %average value
    Data.Header.Spider.av = 0;

    %standard deviation. A value of -1.0 indicates that sig has not been computed previously.
    Data.Header.Spider.sig = -1.0;

    %(obsolete, no longer used).
    Data.Header.Spider.ihist = 0;
    
    %number of pixels per line. (size of volume in y direction)
    Data.Header.Spider.nsam = size(in,2);
    
    %number of records in file header (label).
    Data.Header.Spider.labrec = ceil(256./Data.Header.Spider.nsam);%1024 ./ (Data.Header.Spider.nsam .*4); %???? WHAT IS THIS?
    %if (mod(1024,Data.Header.Spider.nsam .*4) ~= 0)
    %    Data.Header.Spider.labrec = Data.Header.Spider.labrec + 1;
    %end
    
    %flag that tilt angles are present.
    Data.Header.Spider.iangle = 0;

    %tilt angle
    %The angle, offset & scale factor locations contained in the SPIDER header are available to communicate between 
    %different SPIDER operations. Currently they are NOT used in the code distributed with SPIDER, but some outside 
    %labs make extensive use of these positions. The angles are usually in Euler format and are given in degrees.
    Data.Header.Spider.phi = 0;

    %tilt angle
    Data.Header.Spider.theta = 0;

    %tilt angle (also called psi).
    Data.Header.Spider.gamma = 0;

    %x translation
    Data.Header.Spider.xoff = 0;

    %y translation
    Data.Header.Spider.yoff = 0;

    %z translation
    Data.Header.Spider.zoff = 0;

    %scale factor
    Data.Header.Spider.scale = 0;

    %total number of bytes in header.
    Data.Header.Spider.labbyt = Data.Header.Spider.labrec .* Data.Header.Spider.nsam .*4; 

    %record length in bytes.
    Data.Header.Spider.lenbyt = Data.Header.Spider.nsam .*4; 
    
    %This position has a value of 0 in simple 2D or 3D (non-stack) files. 
    %In an "image stack" there is one overall stack header followed by a stack of images in 
    %which each image has its own image header. (An image stack differs from a simple 3D image 
    %in that each stacked image has its own header.) A value of >0 in this position in the overall 
    %stack header indicates a stack of images. A value of <0 in this position in the overall stack 
    %header indicates an indexed stack of images and gives the maximum image number allowed in the index.
    Data.Header.Spider.istack = 0; %only non-stack files supported at the moment.

    %This position is unused now! Prior to release 9.0, a -1 at this location in an overall stack indicated 
    %a valid stack and in the stacked images, a value of 1 indicated that this image was in use (existed).
    Data.Header.Spider.NOTUSED = 0;

    %This position is only used in the overall header for a stacked image file. There, this position contains 
    %the number of the highest image currently used in the stack. This number is updated, if necessary, when an 
    %image is added or deleted from the stack.
    Data.Header.Spider.maxim = 0;

    %This position is only used in a stacked image header. There, this position contains the number of the current image or zero if the image is unused.
    Data.Header.Spider.imgnum = 0;

    %This position is only used in the overall header of indexed stacks. There, this position is the highest index currently in use.
    Data.Header.Spider.lastindx = 0;
    
    %next 2 words are unused
    Data.Header.Spider.dummy1 = 0;
    Data.Header.Spider.dummy2 = 0;
    
    %flag that additional angles are present in header. 1 = one additional rotation is present, 2 = additional rotation that preceeds the rotation that was stored in words 15..20.
    Data.Header.Spider.Kangle = 0;

    %phi1
    Data.Header.Spider.phi1 = 0;

    %theta1
    Data.Header.Spider.theta1 = 0;

    %psi1
    Data.Header.Spider.psi1 = 0;

    %phi2
    Data.Header.Spider.phi2 = 0;

    %theta2
    Data.Header.Spider.theta2 = 0;

    %psi2
    Data.Header.Spider.psi2 = 0;
    
    %reserved for Jose Maria's transforms 
    Data.Header.Spider.tforms = zeros(26,1);
    
    %creation date e.g. 27-MAY-1999 
    Data.Header.Spider.cdat = upper(date);
    
    %creation time e.g. 09:43:19 
    Data.Header.Spider.ctim = datestr(now,13);
    
    %title
    Data.Header.Spider.ctit = blanks(160);
    
else
   Data = in;
   disp('This is already a structure!'); 
end;
