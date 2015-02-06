function tom_spiderwrite(filename,Data)
%TOM_SPIDERWRITE writes out a SPIDER file
%
%   tom_spiderwrite(filename, Data)
%
%PARAMETERS
%
%  INPUT
%   filename            ...
%   Data                ...
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
%   TOM_SPIDERREAD, TOM_SPIDERHEADER, TOM_ISSPIDERFILE
%
%   created by AK 04/25/06
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

if ~ischar(filename)
    error('First argument should be a string containing a valid file name.');
end

if nargin == 1 
    error('Second argument should be a spider structure');
end

if nargin > 2
    error('This function takes only two arguments.');
end

if ~isstruct(Data)
    Data = tom_spiderheader(Data);
end

try
    fid = fopen(filename,'wb');
catch
    error(['Could not open' filename]);
end

%number of slices in volume (size of volume in z direction)
fwrite(fid,Data.Header.Spider.nslice,'float');

%number of rows per slice (size of volume in x direction)
fwrite(fid,Data.Header.Spider.nrow,'float');

%total number of records in the file (unused)
fwrite(fid,Data.Header.Spider.irec,'float');

%(obsolete, unused)
fwrite(fid,Data.Header.Spider.nhistrec,'float');

%file type specifier. Obsolete file types d, 8, 11, 12, 16, -1, -3, -7, and -9 are no longer supported in SPIDER.
%iform  	(type)  	data type
%1 	(r) 	2D image.
%3 	(r) 	3D volume.
%-11 	(fo) 	2D Fourier, mixed radix odd.
%-12 	(fe) 	2D Fourier, mixed radix even.
%-21 	(fo) 	3D Fourier, mixed radix odd.
%-22 	(fe) 	3D Fourier, mixed radix even.
fwrite(fid,Data.Header.Spider.iform,'float');

%imami = maximum/minimum flag. Is set at 0 when the file is created, and at 1 when the maximum, minimum, average, and
%standard deviation have been computed and stored into this header record (see following locations).
fwrite(fid,Data.Header.Spider.imami,'float');

%maximum value
fwrite(fid,Data.Header.Spider.fmax,'float');

%minimum value
fwrite(fid,Data.Header.Spider.fmin,'float');

%average value
fwrite(fid,Data.Header.Spider.av,'float');

%standard deviation. A value of -1.0 indicates that sig has not been computed previously.
fwrite(fid,Data.Header.Spider.sig,'float');

%(obsolete, no longer used).
fwrite(fid,Data.Header.Spider.ihist,'float');

%number of pixels per line. (size of volume in y direction)
fwrite(fid,Data.Header.Spider.nsam,'float');

%number of records in file header (label).
fwrite(fid,Data.Header.Spider.labrec,'float');

%flag that tilt angles are present.
fwrite(fid,Data.Header.Spider.iangle,'float');

%tilt angle
%The angle, offset & scale factor locations contained in the SPIDER header are available to communicate between 
%different SPIDER operations. Currently they are NOT used in the code distributed with SPIDER, but some outside 
%labs make extensive use of these positions. The angles are usually in Euler format and are given in degrees.
fwrite(fid,Data.Header.Spider.phi,'float');

%tilt angle
fwrite(fid,Data.Header.Spider.theta,'float');

%tilt angle (also called psi).
fwrite(fid,Data.Header.Spider.gamma,'float');

%x translation
fwrite(fid,Data.Header.Spider.xoff,'float');

%y translation
fwrite(fid,Data.Header.Spider.yoff,'float');

%z translation
fwrite(fid,Data.Header.Spider.zoff,'float');

%scale factor
fwrite(fid,Data.Header.Spider.scale,'float');

%total number of bytes in header.
fwrite(fid,Data.Header.Spider.labbyt,'float');

%record length in bytes.
fwrite(fid,Data.Header.Spider.lenbyt,'float');

%This position has a value of 0 in simple 2D or 3D (non-stack) files. 
%In an "image stack" there is one overall stack header followed by a stack of images in 
%which each image has its own image header. (An image stack differs from a simple 3D image 
%in that each stacked image has its own header.) A value of >0 in this position in the overall 
%stack header indicates a stack of images. A value of <0 in this position in the overall stack 
%header indicates an indexed stack of images and gives the maximum image number allowed in the index.
fwrite(fid,Data.Header.Spider.istack,'float');

%This position is unused now! Prior to release 9.0, a -1 at this location in an overall stack indicated 
%a valid stack and in the stacked images, a value of 1 indicated that this image was in use (existed).
fwrite(fid,Data.Header.Spider.NOTUSED,'float');

%This position is only used in the overall header for a stacked image file. There, this position contains 
%the number of the highest image currently used in the stack. This number is updated, if necessary, when an 
%image is added or deleted from the stack.
fwrite(fid,Data.Header.Spider.maxim,'float');

%This position is only used in a stacked image header. There, this position contains the number of the current image or zero if the image is unused.
fwrite(fid,Data.Header.Spider.imgnum,'float');

%This position is only used in the overall header of indexed stacks. There,
%this position is the highest index currently in use.
fwrite(fid,Data.Header.Spider.lastindx,'float');

%next 2 words are unused
fwrite(fid,0,'float');
fwrite(fid,0,'float');

%flag that additional angles are present in header. 1 = one additional
%rotation is present, 2 = additional rotation that preceeds the rotation that was stored in words 15..20.
fwrite(fid,Data.Header.Spider.Kangle,'float');

%phi1
fwrite(fid,Data.Header.Spider.phi1,'float');

%theta1
fwrite(fid,Data.Header.Spider.theta1,'float');

%psi1
fwrite(fid,Data.Header.Spider.psi1,'float');

%phi2
fwrite(fid,Data.Header.Spider.phi2,'float');

%theta2
fwrite(fid,Data.Header.Spider.theta2,'float');

%psi2
fwrite(fid,Data.Header.Spider.psi2,'float');

fwrite(fid,0,'float',39);

%creation date e.g. 27-MAY-1999
fwrite(fid,0,'float',134);
try
    fwrite(fid,[Data.Header.Spider.cdat ' '],'char');
catch Me
    fwrite(fid,[Data.Header.Spider.cdat' ' '],'char');
end;

%creation time e.g. 09:43:19
fwrite(fid,Data.Header.Spider.ctim,'char');

%title
fwrite(fid,sprintf('%160s',Data.Header.Spider.ctit),'char');

%finished writing header, see if the position for the values is correct.
fillup = Data.Header.Spider.labbyt-ftell(fid);
if fillup > 0
    fwrite(fid,0,'char',fillup-1);
end

%write the values
fwrite(fid,Data.Value,'float');

fclose(fid);

