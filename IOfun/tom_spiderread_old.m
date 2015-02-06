function Data = tom_spiderread(filename)

%
%  Syntax:
%    Data = tom_spiderread(filename)
%
%  Input:
%    em_name:    ['PATHNAME' 'FILENAME']
%
%  Output:
%    Data:     	Structure of Image Data
%
%  Description
%	Reads an Spider-Image File 
%
%  See Also
%    spiderwrite, emread, emwrite, matread, matwrite
%
% Promotion am MPI fuer Biochemie
% Beginn: 01.07.1998
% Letzte �nderung: 01.07.1998 SN
%
% Stephan Nickell
%   Copyright (c) 2004
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute for Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

global img

% open the stream with the correct format !
fid = fopen(spider_name,'r','ieee-be');

% Word No. 
% read the Header 
nslice = fread(fid,[1],'float')
%   1.nslice = number of slices (planes) in volume (=1 for an image) The value stored in the file
%     is negative on normal VAX/VMS files. (VAX/VMS files created before 1988 had only a
%     single record for the header, irrespective of the record length which was always nsam. For
%     these early VAX/VMS files the first header had a positive value for nslice). The first
%     header position on Unix SPIDER files is always the positive value of nslice since there are
%     no short label files available on UNIX SPIDER. As of August 1996 short label files must
%     be converted to regular files on VAX/VMS using 'CP FROM SHORT' before they can be
%     used. 
nrow = fread(fid,[1],'float')
%   2.nrow = number of rows per slice. 
irec = fread(fid,[1],'float')
irec=(irec-2)./nslice; % check this !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   3.irec = total number of records in the file (unused) 
nhistrec = fread(fid,[1],'float');
%   4.nhistrec = (obsolete) 
iform = fread(fid,[1],'float');
%   5.iform = file type specifier. Currently file types d, 8, c, -, and 16 are not supported in
%     SPIDER. 
%     iform (type) data type 
%         1 (r) 2-d image. 
%         3 (r) 3-d volume. 
%         8 (8) 8 bit black and white. (Sterecon only) 
%         11 (c) 8 bit color. (Sterecon only) 
%         12 (-) 8 bit runlength black and white. (Sterecon only) 
%         16 (16) 16 bit black and white. (Sterecon only) 
%         -1 (f) 2-d Fourier transform. 
%         -3 (f) 3-d Fourier transform. 
%         -7 (q) 3-d Fourier transform, Paul's format. 
%         -9 (s) 3-d simple Fourier (Michael's format). 
%         -11 (fo) 2-d Fourier, mixed radix odd. 
%         -12 (fe) 2-d Fourier, mixed radix even. 
imami = fread(fid,[1],'float');
%   6.imami = maximum/minimum flag. Is set at 0 when the file is created, and at 1 when the
%     maximum, minimum, average, and standard deviation have been computed and stored into
%     this header record (see following locations). 
fmax = fread(fid,[1],'float');
%   7.fmax = maximum value. 
fmin = fread(fid,[1],'float');
%   8.fmin = minimum value. 
av = fread(fid,[1],'float');
%   9.av = average value. 
sig = fread(fid,[1],'float');
%  10.sig = standard deviation. A value of -1.0 Indicates that sig has not been computed
%     previously. 
ihist = fread(fid,[1],'float');
%  11.ihist = (obsolete, no longer used). 
nsam = fread(fid,[1],'float');
%  12.nsam = number of pixels per line. 
headrec = fread(fid,[1],'float');
%  13.headrec = number of records in file header. 
iangle = fread(fid,[1],'float');
%  14.iangle = flag that tilt angles have been filled. 
phi = fread(fid,[1],'float');
%  15.phi = tilt angle. 
theta = fread(fid,[1],'float');
%  16.theta = tilt angle. 
gamma = fread(fid,[1],'float');
%  17.gamma = tilt angle (also called psi). 
xoff = fread(fid,[1],'float');
%  18.xoff = x translation. 
yoff = fread(fid,[1],'float');
%  19.yoff = y translation. 
zoff = fread(fid,[1],'float');
%  20.zoff = z translation. 
scale = fread(fid,[1],'float');
%  21.scale = scale factor. 
labbyt = fread(fid,[1],'float');
%  22.labbyt = total number of bytes in header. 
lenbyt = fread(fid,[1],'float');
%  23.lenbyt = record length in bytes. 
istack = fread(fid,[1],'float');
%  24.istack = indicator that file contains a stack of images. 
inuse = fread(fid,[1],'float');
%  25.inuse = indicator that this image in stack is used (exists) 
maxim = fread(fid,[1],'float');
%  26.maxim = max. image used in stack (in stack header only) 
unused = fread(fid,[4],'float');
%  27.unused 
%  28.unused 
%  29.unused 
%  30.unused 
Kangle = fread(fid,[1],'float');
%  31.Kangle = flag that additional angles are set. 1 = one additional rotation is present, 2 =
%     additional rotation that preceeds the rotation that was stored in words 15..20. For details
%     see manual chapter voceul.man 
phi1 = fread(fid,[1],'float');
%  32.phi1 
theta1 = fread(fid,[1],'float');
%  33.theta1 
psi1 = fread(fid,[1],'float');
%  34.psi1 
phi2 = fread(fid,[1],'float');
%  35.phi2 
theta2 = fread(fid,[1],'float');
%  36.theta2 
psi2 = fread(fid,[1],'float');
%  37.psi2 
fillup = fread(fid,[218],'float');
%38.- 50..76 reserved for Jose Maria's transforms 

%212-214 == cdat = character * 10 containing creation date
%215-216 -- ctim = character * 8 containing creation time
%217-256 -- ctit = character * 160 containing title 
%
% see also http://www.wadsworth.org/spider_doc/spider/docs/image_doc.html
%

Data_read=zeros(irec,nrow,nslice);
D=zeros(irec,nrow,nslice);

% h = waitbar(0,'Load Spider-Data...');

clear Data;

for lauf=1:nslice
%	waitbar(lauf./nslice);

Data_read(:,:,lauf) = fread(fid,[nrow,irec],'float')';

end;

comment=[];
parameter=[];
img=0;
marks=[];
four=[nslice];
angle=0;
notilt=0;

Header=struct('Magic',5,'Size',[nslice nrow irec],'Comment',comment,'Parameter',parameter,'Fillup',fillup);

Data=struct('Value',Data_read,'Header',Header,'Angle',angle,'NoTilt',notilt,'Range_min',fmin,'Range_max',...
	fmax,'Handle',img,'Marker',marks,'FFT',four);

clear Data_read;

%close(h);

fclose(fid);
