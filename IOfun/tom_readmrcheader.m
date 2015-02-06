function [Data] = tom_readmrcheader(varargin)
%TOM_READMRCHAEDER reads the header information of an MRC-file
%
%   [Data] = tom_readmrcheader(varargin)
%
%   out=tom_readmrcheader
%   out=tom_readmrcheader('source_file','le')
%
%PARAMETERS
%
%  INPUT
%   'source_file'       ...
%   le                  ...
%  
%  OUTPUT
%   out                 ...
%
%EXAMPLE
%   i=tom_readmrcheader;
%   A fileselect-box appears and the MRC-file can be picked. Open the file for little-endian only (PC format)
%
%   i=tom_readmrcheader('c:\test\mrs_001.mrc','le'); 
%   Open MRC file mrs_001.mrc in little-endian(PC) format
%
%REFERENCES
%
%SEE ALSO
%   TOM_MRCREAD, TOM_MRC2EM, TOM_ISMRCFILE, TOM_READEMHEADER
%
%   created by SN 09/25/02
%   updated by WDN 05/23/05
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


[comp_typ,maxsize,endian] = computer;
switch endian
    case 'L'
        sysfor='ieee-le';
    case 'B'
        sysfor='ieee-be';
end

if nargin <1 
    [filename, pathname] = uigetfile({'*.mrc';'*.*'}, 'Pick an MRC-file');
    if isequal(filename,0) | isequal(pathname,0) 
        disp('No data loaded.'); return; 
    end;
    mrc_name=[pathname filename];
end;
if nargin==1
    mrc_name=varargin{1};
end

fid = fopen(mrc_name,'r',sysfor);
if fid==-1
    error(['Cannot open: ' mrc_name ' file']); 
end;
MRC.nx = fread(fid,[1],'int');        %integer: 4 bytes
MRC.ny = fread(fid,[1],'int');        %integer: 4 bytes
MRC.nz = fread(fid,[1],'int');        %integer: 4 bytes
MRC.mode = fread(fid,[1],'int');      %integer: 4 bytes
MRC.nxstart= fread(fid,[1],'int');    %integer: 4 bytes
MRC.nystart= fread(fid,[1],'int');    %integer: 4 bytes
MRC.nzstart= fread(fid,[1],'int');    %integer: 4 bytes
MRC.mx= fread(fid,[1],'int');         %integer: 4 bytes
MRC.my= fread(fid,[1],'int');         %integer: 4 bytes
MRC.mz= fread(fid,[1],'int');         %integer: 4 bytes
MRC.xlen= fread(fid,[1],'float');     %float: 4 bytes
MRC.ylen= fread(fid,[1],'float');     %float: 4 bytes
MRC.zlen= fread(fid,[1],'float');     %float: 4 bytes
MRC.alpha= fread(fid,[1],'float');    %float: 4 bytes
MRC.beta= fread(fid,[1],'float');     %float: 4 bytes
MRC.gamma= fread(fid,[1],'float');    %float: 4 bytes
MRC.mapc= fread(fid,[1],'int');       %integer: 4 bytes
MRC.mapr= fread(fid,[1],'int');       %integer: 4 bytes
MRC.maps= fread(fid,[1],'int');       %integer: 4 bytes
MRC.amin= fread(fid,[1],'float');     %float: 4 bytes
MRC.amax= fread(fid,[1],'float');     %float: 4 bytes
MRC.amean= fread(fid,[1],'float');    %float: 4 bytes
MRC.ispg= fread(fid,[1],'short');     %integer: 2 bytes
MRC.nsymbt = fread(fid,[1],'short');  %integer: 2 bytes
MRC.next = fread(fid,[1],'int');      %integer: 4 bytes
MRC.creatid = fread(fid,[1],'short'); %integer: 2 bytes
MRC.unused1 = fread(fid,[30]);        %not used: 30 bytes
MRC.nint = fread(fid,[1],'short');    %integer: 2 bytes
MRC.nreal = fread(fid,[1],'short');   %integer: 2 bytes
MRC.unused2 = fread(fid,[28]);        %not used: 28 bytes
MRC.idtype= fread(fid,[1],'short');   %integer: 2 bytes
MRC.lens=fread(fid,[1],'short');      %integer: 2 bytes
MRC.nd1=fread(fid,[1],'short');       %integer: 2 bytes
MRC.nd2 = fread(fid,[1],'short');     %integer: 2 bytes
MRC.vd1 = fread(fid,[1],'short');     %integer: 2 bytes
MRC.vd2 = fread(fid,[1],'short');     %integer: 2 bytes
for i=1:6                             %24 bytes in total
    MRC.tiltangles(i)=fread(fid,[1],'float');%float: 4 bytes
end
MRC.xorg = fread(fid,[1],'float');    %float: 4 bytes
MRC.yorg = fread(fid,[1],'float');    %float: 4 bytes
MRC.zorg = fread(fid,[1],'float');    %float: 4 bytes
MRC.cmap = fread(fid,[4],'char');     %Character: 4 bytes
MRC.stamp = fread(fid,[4],'char');    %Character: 4 bytes
MRC.rms=fread(fid,[1],'float');       %float: 4 bytes
MRC.nlabl = fread(fid,[1],'int');     %integer: 4 bytes
MRC.labl = fread(fid,[800],'char');   %Character: 800 bytes

Extended.magnification(1)=0;
Extended.exp_time(1)=0;
Extended.pixelsize(1)=0;
Extended.defocus(1)=0;
Extended.a_tilt(1:MRC.nz)=0;
Extended.tiltaxis(1)=0;
if MRC.next~=0%Extended Header
    nbh=MRC.next./128;%128=lengh of FEI extended header
    if nbh==1024%FEI extended Header
        for lauf=1:nbh
            Extended.a_tilt(lauf)= fread(fid,[1],'float');        %float: 4 bytes
            Extended.b_tilt(lauf)= fread(fid,[1],'float');        %float: 4 bytes
            Extended.x_stage(lauf)= fread(fid,[1],'float');       %float: 4 bytes
            Extended.y_stage(lauf)=fread(fid,[1],'float');        %float: 4 bytes
            Extended.z_stage(lauf)=fread(fid,[1],'float');        %float: 4 bytes
            Extended.x_shift(lauf)=fread(fid,[1],'float');        %float: 4 bytes
            Extended.y_shift(lauf)=fread(fid,[1],'float');        %float: 4 bytes
            Extended.defocus(lauf)=fread(fid,[1],'float');        %float: 4 bytes
            Extended.exp_time(lauf)=fread(fid,[1],'float');       %float: 4 bytes
            Extended.mean_int(lauf)=fread(fid,[1],'float');       %float: 4 bytes
            Extended.tiltaxis(lauf)=fread(fid,[1],'float');       %float: 4 bytes
            Extended.pixelsize(lauf)=fread(fid,[1],'float');      %float: 4 bytes
            Extended.magnification(lauf)=fread(fid,[1],'float');  %float: 4 bytes
            fseek(fid,128-52,0);
            %position = ftell(fid)
        end
    else %IMOD extended Header
        fseek(fid,MRC.next,'cof');%go to end end of extended Header
    end
end
fclose(fid);
%Header=struct('Size',[MRC.nx MRC.ny MRC.nz]','MRC',MRC);
Header=struct(...
    'Voltage',0,...
    'Cs',0,...
    'Aperture',0,...
    'Magnification',Extended.magnification(1),...
    'Postmagnification',0,...
    'Exposuretime',Extended.exp_time(1),...
    'Objectpixelsize',Extended.pixelsize(1).*1e9,...
    'Microscope',0,...
    'Pixelsize',0,...
    'CCDArea',0,...
    'Defocus',Extended.defocus(1),...
    'Astigmatism',0,...
    'AstigmatismAngle',0,...
    'FocusIncrement',0,...
    'CountsPerElectron',0,...
    'Intensity',0,...
    'EnergySlitwidth',0,...
    'EnergyOffset',0,... 
    'Tiltangle',Extended.a_tilt(1:MRC.nz),...
    'Tiltaxis',Extended.tiltaxis(1),...
    'Username',num2str(zeros(20,1)),...
    'Date',num2str(zeros(8)),...
    'Size',[MRC.nx,MRC.ny,MRC.nz],...
    'Comment',num2str(zeros(80,1)),...
    'Parameter',num2str(zeros(40,1)),...
    'Fillup',num2str(zeros(256,1)),...
    'Filename',mrc_name,...
    'Marker_X',0,...
    'Marker_Y',0,...
    'MRC',MRC);
Data=struct('Header',Header);

