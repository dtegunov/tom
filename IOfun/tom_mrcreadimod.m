function [Data] = tom_mrcreadimod(mrc_name,format)
%TOM_MRCREADIMOD reads data in MRC-file format imod style
%
%   [Data] = tom_mrcreadimod(mrc_name,format)
%
%PARAMETERS
%
%  INPUT
%   filelabel           ...
%   threshold           ...
%   label               ...
%   color               ...
%   transformmatrix     ...
%   iconposition        ...
%   host                ...
%  
%  OUTPUT
%   data		...
%
%DESCRIPTION
% Info copied from: http://bio3d.colorado.edu/imod/doc/mrc_format.txt
%
% The MRC file format used by IMOD.
% 
% The MRC header. length 1024 bytes
% 
% SIZE DATA    NAME	Description
% 
% 4    int     nx;	Number of Columns 
% 4    int     ny;        Number of Rows
% 4    int     nz;        Number of Sections.
% 
% 4    int     mode;      Types of pixel in image.  Values used by IMOD:
% 			 0 = unsigned bytes, 
% 			 1 = signed short integers (16 bits), 
% 			 2 = float,  
% 			 3 = short * 2, (used for complex data)
% 			 4 = float * 2, (used for complex data)
% 			16 = unsigned char * 3 (used for rgb data)
% 
% 4    int     nxstart;     Starting point of sub image.  
% 4    int     nystart;
% 4    int     nzstart;
% 
% 4    int     mx;         Grid size in X, Y, and Z       
% 4    int     my;
% 4    int     mz;
% 
% 4    float   xlen;       Cell size; pixel spacing = xlen/mx   
% 4    float   ylen;
% 4    float   zlen;
% 
% 4    float   alpha;      cell angles 
% 4    float   beta;
% 4    float   gamma;
% 
% 			 Ignored by imod.
% 4    int     mapc;       map coloumn 1=x,2=y,3=z.       
% 4    int     mapr;       map row     1=x,2=y,3=z.       
% 4    int     maps;       map section 1=x,2=y,3=z.       
% 
%                          These need to be set for proper scaling of
% 	 		 non byte data.
% 4    float   amin;       Minimum pixel value.
% 4    float   amax;       Maximum pixel value.
% 4    float   amean;      Mean pixel value.
% 
% 2    short   ispg;       image type 
% 2    short   nsymbt;     space group number 
% 4    int     next;       number of bytes in extended header 
% 2    short   creatid;    Creator ID 
% 30   ---     extra data (not used)
% 
%                          These two values specify the structure of data in the
%                          extended header; their meaning depend on whether the
%                          extended header has the Agard format, a series of
%                          4-byte integers then real numbers, or has data 
%                          produced by SerialEM, a series of short integers
% 2    short   nint;       Number of integers per section (Agard format) or
%                          number of bytes per section (SerialEM format)
% 2    short   nreal;      Number of reals per section (Agard format) or
%                          flags for which types of short data (SerialEM format):
%                          1 = tilt angle * 100  (2 bytes)
%                          2 = piece coordinates for montage  (6 bytes)
%                          4 = Stage position * 25    (4 bytes)
%                          8 = Magnification / 100 (2 bytes)
%                          16 = Intensity * 25000  (2 bytes)
%                          32, 128, 512: Reserved for 4-byte items
%                          64, 256, 1024: Reserved for 2-byte items
%                          If the number of bytes implied by these flags does
%                          not add up to the value in nint, then nint and nreal
%                          are interpreted as ints and reals per section
% 
% 28   ---     extra data (not used)
%                       Explanation of type of data.
% 2    short   idtype;  ( 0 = mono, 1 = tilt, 2 = tilts, 3 = lina, 4 = lins)
% 2    short   lens;
% 2    short   nd1;	for idtype = 1, nd1 = axis (1, 2, or 3)     
% 2    short   nd2;
% 2    short   vd1;                       vd1 = 100. * tilt increment
% 2    short   vd2;                       vd2 = 100. * starting angle
% 
%         		Used to rotate model to match new rotated image.
% 24   float   tiltangles[6];  0,1,2 = original:  3,4,5 = current 
% 
% OLD-STYLE MRC HEADER - IMOD 2.6.19 and below:
% 2    short   nwave;     # of wavelengths and values
% 2    short   wave1;
% 2    short   wave2;
% 2    short   wave3;
% 2    short   wave4;
% 2    short   wave5;
% 
% 4    float   zorg;      Origin of image.  Used to auto translate model
% 4    float   xorg;      to match a new image that has been translated.
% 4    float   yorg;
% 
% NEW-STYLE MRC image2000 HEADER - IMOD 2.6.20 and above:
% 4    float   xorg;      Origin of image.  Used to auto translate model
% 4    float   yorg;      to match a new image that has been translated.
% 4    float   zorg;
% 
% 4    char    cmap;      Contains "MAP "
% 4    char    stamp;     First byte has 17 for big- or 68 for little-endian
% 4    float   rms;
% 
% ALL HEADERS:
% 4    int     nlabl;  	Number of labels with useful data.
% 800  char[10][80]    	10 labels of 80 charactors.
% ------------------------------------------------------------------------
% 
% Total size of header is 1024 bytes plus the size of the extended header.
%  
% Image data follows with the origin in the lower left corner,
% looking down on the volume.
% 
% The size of the image is nx * ny * nz * (mode data size).
%
%EXAMPLE
%   ... = tom_mrcreadimod(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%    tom_emread, tom_mrc2emseries, tom_mrc2emstack
%
%   created by SN 09/25/02
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


error(nargchk(0,2,nargin))
if nargin <1 
    [filename, pathname] = uigetfile({'*.mrc';'*.*'}, 'Pick an MRC-file');
    if isequal(filename,0) | isequal(pathname,0) 
        disp('No data loaded.'); return; 
    end;
    em_name=[pathname filename];
    format='le';
end;
if nargin==1 %default format: little endian (le)
    format='le';
end

if isequal(format,'le');
    fid = fopen(mrc_name,'r','ieee-le');
else
    fid = fopen(mrc_name,'r','ieee-be');
end;

if fid==-1
    error(['Cannot open: ' mrc_name ' file']); 
end;

%Header=struct(...
%    'Nx','','Ny','','Nz','',...
%    'Mode','',...
%    'Nxstart','','Nystart','','Nzstart','',...
%    'Mx','','My','','Mz','',...
%    'Xlen','','Ylen','','Zlen','',...
%    'Alpha','','Beta','','Gamma','',...
%    'Mapc','','Mapr','','Maps','',...
%    'Amin','','Amax','','Amean','',...
%    'Ispg','',...
%    'Nsymt','',...
%    'Next','',...
%    'Creatid','',...
%    'Extra_data1','',...
%    'Nint','','Nreal','',...
%    'Extra_data2','',...
%    'Idtype','',...
%    'Lens','',...
%    'Nd1','','Nd2','',...
%    'Vd1','','Vd2','',...
%    'Tiltangles','',...
%    'Xorg','','Yorg','','Zorg','',...
%    'Cmap','','Stamp','','Rms','',...
%    'Nlabl','',...
%    'Comment','');

Header.nx = fread(fid,[1],'long');
Header.ny = fread(fid,[1],'long');
Header.nz = fread(fid,[1],'long');

Header.mode = fread(fid,[1],'long');

Header.nxstart = fread(fid,[1],'long');
Header.nystart = fread(fid,[1],'long');
Header.nzstart = fread(fid,[1],'long');

Header.mx = fread(fid,[1],'long');
Header.my = fread(fid,[1],'long');
Header.mz = fread(fid,[1],'long');

Header.xlen = fread(fid,[1],'float');
Header.ylen = fread(fid,[1],'float');
Header.zlen = fread(fid,[1],'float');

Header.alpha = fread(fid,[1],'float');
Header.beta = fread(fid,[1],'float');
Header.gamma = fread(fid,[1],'float');

Header.mapc = fread(fid,[1],'long');
Header.mapr = fread(fid,[1],'long');
Header.maps = fread(fid,[1],'long');

Header.amin = fread(fid,[1],'float');
Header.amax = fread(fid,[1],'float');
Header.amean = fread(fid,[1],'float');

Header.ispg = fread(fid,[1],'short');
Header.nsymbt = fread(fid,[1],'short');
Header.next = fread(fid,[1],'long');
Header.creatid = fread(fid,[1],'short');
Header.extra_data1=fread(fid,[30],'char');
Header.nint=fread(fid,[1],'short');
Header.nreal=fread(fid,[1],'short');
Header.extra_data2=fread(fid,[28],'char');
Header.idtype=fread(fid,[1],'short');
Header.lens=fread(fid,[1],'short');
Header.nd1=fread(fid,[1],'short');
Header.nd2 = fread(fid,[1],'short');
Header.vd1 = fread(fid,[1],'short');
Header.vd2 = fread(fid,[1],'short');
for i=1:6
    Header.tiltangles(i)=fread(fid,[1],'float');
end
%new_style
    Header.xorg=fread(fid,[1],'float');
    Header.yorg=fread(fid,[1],'float');
    Header.zorg=fread(fid,[1],'float');
    Header.cmap=fread(fid,[1],'char');
    Header.stamp=fread(fid,[1],'char');
    Header.rms=fread(fid,[1],'float');

Header.nlabl=fread(fid,[1],'long');
Header.label=fread(fid,[800],'char');
Header.ExtentedHeader=fread(fid,[Header.next./4],'float');

fseek(fid,0,'eof'); %go to the end of file
fprintf('loading');
for i=1:Header.nz
    fprintf('.');
    if Header.mode==0
        beval=Header.nx*Header.ny*Header.nz;
        fseek(fid,-beval,0); %go to the beginning of the values
        Data_read(:,:,i) = fread(fid,[Header.nx,Header.ny],'int8');        
    elseif ( (Header.mode==1 || Header.mode==6) )
        beval=Header.nx*Header.ny*Header.nz*2;
        fseek(fid,-beval,0); %go to the beginning of the values
        Data_read(:,:,i) = fread(fid,[Header.nx,Header.ny],'int16');        
    elseif (Header.mode==2 || Header.mode==6)
       beval=Header.nx*Header.ny*Header.nz*4;
       fseek(fid,-beval,0); %go to the beginning of the values
        Data_read(:,:,i) = fread(fid,[Header.nx,Header.ny],'float');        
    else
        error(['Sorry, i cannot read this as an MRC-File !!!']);
        Data_read=[];
    end;
    Data_read2(:,:,i)=flipdim(Data_read(:,:,i),2);
end;
disp('done');
fclose(fid);
%MRC=struct('Comment',Header.label,'Parameter',Parameter);
%Header=struct('Size',[Header.nx Header.ny Header.nz]','MRC',MRC);

Data=struct('Value',Data_read2,'Header',Header);

clear Data_read;
