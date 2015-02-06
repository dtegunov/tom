function Data = tom_readspiderheader(filename)

try
    fid = fopen(filename,'rb','ieee-le');
    testval = fread(fid,1,'float');
    fseek(fid,0,-1);
    if testval < 1 || testval > 10000
        fclose(fid);
        fid = fopen(filename,'rb','ieee-be');
    end
    
catch
    error(['Could not open ' filename]);
end


Data = struct();
Data.Header = struct();
Data.Header.Spider = struct();

%number of slices in volume (size of volume in z direction)
Data.Header.Spider.nslice = fread(fid,1,'float');

%number of rows per slice (size of volume in x direction)
Data.Header.Spider.nrow = fread(fid,1,'float');

%total number of records in the file (unused)
Data.Header.Spider.irec = fread(fid,1,'float');

%(obsolete, unused)
Data.Header.Spider.nhistrec = fread(fid,1,'float');

%file type specifier. Obsolete file types d, 8, 11, 12, 16, -1, -3, -7, and -9 are no longer supported in SPIDER.
%iform  	(type)  	data type
%1 	(r) 	2D image.
%3 	(r) 	3D volume.
%-11 	(fo) 	2D Fourier, mixed radix odd.
%-12 	(fe) 	2D Fourier, mixed radix even.
%-21 	(fo) 	3D Fourier, mixed radix odd.
%-22 	(fe) 	3D Fourier, mixed radix even.
Data.Header.Spider.iform = fread(fid,1,'float');

%imami = maximum/minimum flag. Is set at 0 when the file is created, and at 1 when the maximum, minimum, average, and
%standard deviation have been computed and stored into this header record (see following locations).
Data.Header.Spider.imami = fread(fid,1,'float');

%maximum value
Data.Header.Spider.fmax = fread(fid,1,'float');

%minimum value
Data.Header.Spider.fmin = fread(fid,1,'float');

%average value
Data.Header.Spider.av = fread(fid,1,'float');

%standard deviation. A value of -1.0 indicates that sig has not been computed previously.
Data.Header.Spider.sig = fread(fid,1,'float');

%(obsolete, no longer used).
Data.Header.Spider.ihist = fread(fid,1,'float');

%number of pixels per line. (size of volume in y direction)
Data.Header.Spider.nsam = fread(fid,1,'float');

%number of records in file header (label).
Data.Header.Spider.labrec = fread(fid,1,'float');

%flag that tilt angles are present.
Data.Header.Spider.iangle = fread(fid,1,'float');

%tilt angle
%The angle, offset & scale factor locations contained in the SPIDER header are available to communicate between 
%different SPIDER operations. Currently they are NOT used in the code distributed with SPIDER, but some outside 
%labs make extensive use of these positions. The angles are usually in Euler format and are given in degrees.
Data.Header.Spider.phi = fread(fid,1,'float');

%tilt angle
Data.Header.Spider.theta = fread(fid,1,'float');

%tilt angle (also called psi).
Data.Header.Spider.gamma = fread(fid,1,'float');

%x translation
Data.Header.Spider.xoff = fread(fid,1,'float');

%y translation
Data.Header.Spider.yoff = fread(fid,1,'float');

%z translation
Data.Header.Spider.zoff = fread(fid,1,'float');

%scale factor
Data.Header.Spider.scale = fread(fid,1,'float');

%total number of bytes in header.
Data.Header.Spider.labbyt = fread(fid,1,'float');

%record length in bytes.
Data.Header.Spider.lenbyt = fread(fid,1,'float');

%This position has a value of 0 in simple 2D or 3D (non-stack) files. 
%In an "image stack" there is one overall stack header followed by a stack of images in 
%which each image has its own image header. (An image stack differs from a simple 3D image 
%in that each stacked image has its own header.) A value of >0 in this position in the overall 
%stack header indicates a stack of images. A value of <0 in this position in the overall stack 
%header indicates an indexed stack of images and gives the maximum image number allowed in the index.
Data.Header.Spider.istack = fread(fid,1,'float');

%This position is unused now! Prior to release 9.0, a -1 at this location in an overall stack indicated 
%a valid stack and in the stacked images, a value of 1 indicated that this image was in use (existed).
Data.Header.Spider.NOTUSED = fread(fid,1,'float');

%This position is only used in the overall header for a stacked image file. There, this position contains 
%the number of the highest image currently used in the stack. This number is updated, if necessary, when an 
%image is added or deleted from the stack.
Data.Header.Spider.maxim = fread(fid,1,'float');

%This position is only used in a stacked image header. There, this position contains the number of the current image or zero if the image is unused.
Data.Header.Spider.imgnum = fread(fid,1,'float');

%This position is only used in the overall header of indexed stacks. There, this position is the highest index currently in use.
Data.Header.Spider.lastindx = fread(fid,1,'float');

%next 2 words are unused
Data.Header.Spider.dummy1 = fread(fid,1,'float');
Data.Header.Spider.dummy2 = fread(fid,1,'float');

%flag that additional angles are present in header. 1 = one additional rotation is present, 2 = additional rotation that preceeds the rotation that was stored in words 15..20.
Data.Header.Spider.Kangle = fread(fid,1,'float');

%phi1
Data.Header.Spider.phi1 = fread(fid,1,'float');

%theta1
Data.Header.Spider.theta1 = fread(fid,1,'float');

%psi1
Data.Header.Spider.psi1 = fread(fid,1,'float');

%phi2
Data.Header.Spider.phi2 = fread(fid,1,'float');

%theta2
Data.Header.Spider.theta2 = fread(fid,1,'float');

%psi2
Data.Header.Spider.psi2 = fread(fid,1,'float');

%reserved for Jose Maria's transforms 
fread(fid,14,'float');
Data.Header.Spider.tforms = fread(fid,26,'float');

%creation date e.g. 27-MAY-1999 
fread(fid,135,'float');
Data.Header.Spider.cdat = fread(fid,12,'char');

%creation time e.g. 09:43:19 
Data.Header.Spider.ctim = fread(fid,8,'char');

%title
Data.Header.Spider.ctit = fread(fid,160,'char');

%construct general header

Data.Header.Size = [Data.Header.Spider.nsam; Data.Header.Spider.nrow; Data.Header.Spider.nslice];
Data.Header.Filename=filename; % changed by SN

fclose(fid);