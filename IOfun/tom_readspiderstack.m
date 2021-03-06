function Data = tom_readspiderstack(filename,verbose)

if nargin < 2
    verbose = 0;
end

if nargin <1 
    [filename, pathname] = uigetfile({'*.spi';'*.*'}, 'Pick a spider file');
    if isequal(filename,0) || isequal(pathname,0); disp('No data loaded.'); return; end;
    filename=[pathname filename];
end;

try
    fid = fopen(filename,'rb','ieee-le');
    testval = fread(fid,1,'float');
    fseek(fid,0,'bof');
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
Data.Header.Particles = struct();

%read overall header
Data.Header.Spider = tom_readspiderheader(fid);

if Data.Header.Spider.istack == 0
    error('This spider file is not a stack of particles, use tom_spiderread.');
end

Data.Header.Size = [Data.Header.Spider.nsam; Data.Header.Spider.nrow; Data.Header.Spider.maxim];
Data.Header.Filename=filename; % changed by SN

Data.Value = zeros(Data.Header.Spider.nsam, Data.Header.Spider.nrow,Data.Header.Spider.maxim,'single');
if verbose == 1
    disp(['number of particles in stack: ' num2str(Data.Header.Spider.maxim)]);
end

%loop over all particles in stack and read them
for i=1:Data.Header.Spider.maxim
    Data.Header.Particles(i).Header = tom_readspiderheader(fid);
    tmp = fread(fid,Data.Header.Spider.nrow .* Data.Header.Spider.nsam,'float');
    Data.Value(:,:,i) = reshape(single(tmp), Data.Header.Spider.nsam, Data.Header.Spider.nrow);
    if verbose == 1 && mod(i,100) == 0
        disp([num2str(i) ' particles read']);
    end
end

if verbose == 1
    disp('Finished.');
end


fclose(fid);



function Data = tom_readspiderheader(fid)

Data = struct();

pos_before = ftell(fid);

%number of slices in volume (size of volume in z direction)
Data.nslice = fread(fid,1,'float');

%number of rows per slice (size of volume in x direction)
Data.nrow = fread(fid,1,'float');

%total number of records in the file (unused)
Data.irec = fread(fid,1,'float');

%(obsolete, unused)
Data.nhistrec = fread(fid,1,'float');

%file type specifier. Obsolete file types d, 8, 11, 12, 16, -1, -3, -7, and -9 are no longer supported in SPIDER.
%iform  	(type)  	data type
%1 	(r) 	2D image.
%3 	(r) 	3D volume.
%-11 	(fo) 	2D Fourier, mixed radix odd.
%-12 	(fe) 	2D Fourier, mixed radix even.
%-21 	(fo) 	3D Fourier, mixed radix odd.
%-22 	(fe) 	3D Fourier, mixed radix even.
Data.iform = fread(fid,1,'float');

%imami = maximum/minimum flag. Is set at 0 when the file is created, and at 1 when the maximum, minimum, average, and
%standard deviation have been computed and stored into this header record (see following locations).
Data.imami = fread(fid,1,'float');

%maximum value
Data.fmax = fread(fid,1,'float');

%minimum value
Data.fmin = fread(fid,1,'float');

%average value
Data.av = fread(fid,1,'float');

%standard deviation. A value of -1.0 indicates that sig has not been computed previously.
Data.sig = fread(fid,1,'float');

%(obsolete, no longer used).
Data.ihist = fread(fid,1,'float');

%number of pixels per line. (size of volume in y direction)
Data.nsam = fread(fid,1,'float');

%number of records in file header (label).
Data.labrec = fread(fid,1,'float');

%flag that tilt angles are present.
Data.iangle = fread(fid,1,'float');

%tilt angle
%The angle, offset & scale factor locations contained in the SPIDER header are available to communicate between 
%different SPIDER operations. Currently they are NOT used in the code distributed with SPIDER, but some outside 
%labs make extensive use of these positions. The angles are usually in Euler format and are given in degrees.
Data.phi = fread(fid,1,'float');

%tilt angle
Data.theta = fread(fid,1,'float');

%tilt angle (also called psi).
Data.gamma = fread(fid,1,'float');

%x translation
Data.xoff = fread(fid,1,'float');

%y translation
Data.yoff = fread(fid,1,'float');

%z translation
Data.zoff = fread(fid,1,'float');

%scale factor
Data.scale = fread(fid,1,'float');

%total number of bytes in header.
Data.labbyt = fread(fid,1,'float');

%record length in bytes.
Data.lenbyt = fread(fid,1,'float');

%This position has a value of 0 in simple 2D or 3D (non-stack) files. 
%In an "image stack" there is one overall stack header followed by a stack of images in 
%which each image has its own image header. (An image stack differs from a simple 3D image 
%in that each stacked image has its own header.) A value of >0 in this position in the overall 
%stack header indicates a stack of images. A value of <0 in this position in the overall stack 
%header indicates an indexed stack of images and gives the maximum image number allowed in the index.
Data.istack = fread(fid,1,'float');

%This position is unused now! Prior to release 9.0, a -1 at this location in an overall stack indicated 
%a valid stack and in the stacked images, a value of 1 indicated that this image was in use (existed).
Data.NOTUSED = fread(fid,1,'float');

%This position is only used in the overall header for a stacked image file. There, this position contains 
%the number of the highest image currently used in the stack. This number is updated, if necessary, when an 
%image is added or deleted from the stack.
Data.maxim = fread(fid,1,'float');

%This position is only used in a stacked image header. There, this position contains the number of the current image or zero if the image is unused.
Data.imgnum = fread(fid,1,'float');

%This position is only used in the overall header of indexed stacks. There, this position is the highest index currently in use.
Data.lastindx = fread(fid,1,'float');

%next 2 words are unused
Data.dummy1 = fread(fid,1,'float');
Data.dummy2 = fread(fid,1,'float');

%flag that additional angles are present in header. 1 = one additional rotation is present, 2 = additional rotation that preceeds the rotation that was stored in words 15..20.
Data.Kangle = fread(fid,1,'float');

%phi1
Data.phi1 = fread(fid,1,'float');

%theta1
Data.theta1 = fread(fid,1,'float');

%psi1
Data.psi1 = fread(fid,1,'float');

%phi2
Data.phi2 = fread(fid,1,'float');

%theta2
Data.theta2 = fread(fid,1,'float');

%psi2
Data.psi2 = fread(fid,1,'float');

%reserved for Jose Maria's transforms 
fread(fid,14,'float');
Data.tforms = fread(fid,26,'float');

%creation date e.g. 27-MAY-1999 
fread(fid,135,'float');
Data.cdat = fread(fid,12,'char');

%creation time e.g. 09:43:19 
Data.ctim = fread(fid,8,'char');

%title
Data.ctit = fread(fid,160,'char');

pos_after = ftell(fid);
headersize = pos_after-pos_before;

fseek(fid,Data.labbyt-headersize,'cof');