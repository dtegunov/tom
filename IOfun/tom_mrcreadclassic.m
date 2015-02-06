function [Data] = tom_mrcreadclassic(mrc_name,format)
%TOM_MRCREADCLASSIC reads data in MRC-file format
%
%   [Data] = tom_mrcreadclassic(mrc_name,format)
%
%PARAMETERS
%
%  INPUT
%   mrc_name            filename
%   format              machine-byte-order  set format 'le' for little endian files (PC,Linux) 
%                       or 'be' for big endian (SGI,Mac)
%  
%  OUTPUT
%   Data                Structure of Image Data
%   Data.Value          Raw data of image, or stack
%   Data.Header         Header information
%
%  Reads an MRC-Image File
%  a raw format with a 1024 Byte header.
%  Structure of MRC-data files:
%  Header with a length = 1024 bytes, organized as 56 LONG words followed
%  by space for 10 80 byte text labels.
%     
%      1      NX              # of Columns    (fastest changing in map)
%      2      NY              # of Rows
%      3      NZ              # of Sections   (slowest changing in map)
%      4      MODE            Data type
%                               0 = Image     stored as Integer*1
%                               1 = Image     stored as Integer*2
%                               2 = Image     stored as Reals
%                               3 = Transform stored as Complex Integer*2
%                               4 = Transform stored as Complex Reals
%     
%      5      NXSTART         Number of first COLUMN  in map (Default = 0)
%      6      NYSTART         Number of first ROW     in map      "
%      7      NZSTART         Number of first SECTION in map      "
%      8      MX              Number of intervals along X
%      9      MY              Number of intervals along Y
%     10      MZ              Number of intervals along Z
%     11      X length        Cell Dimensions (Angstroms)
%     12      Y length                     "
%     13      Z length                     "
%     14      Alpha           Cell Angles     (Degrees)
%     15      Beta                         "
%     16      Gamma                        "
%     17      MAPC            Which axis corresponds to Columns  (1,2,3 for X,Y,Z)
%     18      MAPR            Which axis corresponds to Rows     (1,2,3 for X,Y,Z)
%     19      MAPS            Which axis corresponds to Sections (1,2,3 for X,Y,Z)
%     20      AMIN            Minimum density value
%     21      AMAX            Maximum density value
%     22      AMEAN           Mean    density value    (Average)
%     23      ISPG            Space group number       (0 for images)
%     24      NSYMBT          Number of bytes used for storing symmetry operators
%     25      EXTRA           Extra, user defined storage space. 29 words max.
%     .          "
%     .          "
%     .          "   (all set to zero by default)
%     .          "
%     53         "
%     54      XORIGIN         X origin
%     55      YORIGIN         Y origin
%     56      NLABL           Number of labels being used
%     57-256  LABEL(20,10)    10  80 character text labels (ie. A4 format)
%     
%     Symmetry records follow - if any - stored as text as in International
%     Tables, operators separated by * and grouped into 'lines' of 80
%     characters (ie. symmetry operators do not cross the ends of the
%     80-character 'lines' and the 'lines' do not terminate in a *).
%     
%                     Data records follow.
%
%     In general, the X-Y origin is taken as 0,0 being in the
%     lower-left corner of the image AND the first data point in the file
%     (normally corresponding to array element 1,1).
%
%                        Top of Image             1,NY     Array    NX,NY
%                ^ !                                ^ !
%                ! !                                ! !
%                ! !                                ! !
%                y !                                y !
%                  !_____________   x-->              !_______________  x-->
%                0,0                                 1,1             NX,1
%
%
%     Complex images are not supported!
%
%EXAMPLE
%    i=tom_mrcread;
%    a fileselect-box appears and the EM-file can be picked
%    i=data=tom_mrcreadclassic('test.mrc','be');
%    figure; tom_imagesc(i.Value); 
%
%REFERENCES
%
%SEE ALSO
%   tom_emread, tom_mrc2emseries, tom_mrc2emstack
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

if isequal(format,'le');
    fid = fopen(mrc_name,'r','ieee-le');
else
    fid = fopen(mrc_name,'r','ieee-be');
end;

if fid==-1
    error(['Cannot open: ' mrc_name ' file']); 
end;
NX = fread(fid,[1],'long');
NY = fread(fid,[1],'long');
NZ = fread(fid,[1],'long');
MODE = fread(fid,[1],'long');
PARAMETER=fread(fid,[52],'long');
COMMENT=fread(fid,[800],'char');
Data_read=zeros(NX,NY,NZ);
%Extended=fread(fid,[(132096-1024)./4],'float');

for i=1:NZ
    
    if MODE==0
        Data_read(:,:,i) = fread(fid,[NX,NY],'int8');
        
    elseif MODE==1
        Data_read(:,:,i) = fread(fid,[NX,NY],'int16');
        
    elseif MODE==2
        Data_read(:,:,i) = fread(fid,[NX,NY],'float');
        
    else
        error(['Sorry, i cannot read this as an MRC-File !!!']);
        Data_read=[];
    end;
    Data_read2(:,:,i)=flipdim(Data_read(:,:,i),2);
end;

fclose(fid);
MRC=struct('Comment',COMMENT,'Parameter',PARAMETER);
Header=struct('Size',[NX NY NZ]','MRC',MRC);
Data=struct('Value',Data_read2,'Header',Header);

clear Data_read;
