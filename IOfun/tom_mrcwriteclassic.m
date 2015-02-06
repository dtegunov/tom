function tom_mrcwriteclassic(varargin)
%TOM_MRCWRITECLASSIC writes data in an MRC-file format
%
%   tom_mrcwriteclassic(varargin)
%
%	 Writes an MRC-Image File a raw format with
%    a 1024 Byte header. If input is not a structure 
%    a default header is created
%
%	 tom_mrcwriteclassic(em_name,Data,format)
%
%PARAMETERS
%
%  INPUT
%   em_name             ['PATHNAME' 'FILENAME'] of the output file
%   Data                Structure of Image Data
%   Data.Value          Raw data of image, or stack
%   Data.Header         Header information
%   format              'le' for little, 'be' for big endian format
%  
%  OUTPUT
%
%EXAMPLE
%   a fileselect-box appears and the data can be saved
%   with the selected filename as little endian
%   load clown;
%   tom_mrcwriteclassic(X,'le');
%
%   Save X in test.mrc as big endian format
%   load clown
%   tom_mrcwriteclassic('test.mrc',X,'be');
%
%REFERENCES
%
%SEE ALSO
%   TOM_MRCREAD, TOM_EMWRITE
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



if nargin<=1
    error(['Data not specified (e.g. tom_mrcwrite(out)']);
elseif nargin==2
    Data=varargin{1};
    format=varargin(2);    
    %Data=em_name;
    [filename, pathname] = uiputfile({'*.mrc';'*.*'}, 'Save as MRC-file');
    if isequal(filename,0) | isequal(pathname,0) disp('Data not saved.'); return; end;
    em_name=[pathname filename];
    if isempty(findstr('.mrc',em_name))
        em_name=[em_name '.mrc'];
    end
elseif nargin==3
    em_name=varargin{1};
    Data=varargin{2};
    format=varargin(3);    
end;

%if nargin <1 error(['Data not specified (e.g. tom_mrcwrite(out)']);  end;

if isstruct(Data)
    if isa(Data.Value,'double') | isa(Data.Value,'single')
        MODE=2;
    end;
    if isa(Data.Value,'int16') 
        MODE=1;
    end;
    if isa(Data.Value,'int8') 
        MODE=0;
    end;
    if size(Data.Header.Comment,1)<800
        fillup=char(zeros(800-size(Data.Header.Comment,1),1));
    Data.Header.Comment=[Data.Header.Comment' fillup']';
    end;
    if size(Data.Header.Parameter,1)<52
        fillup=zeros(52-size(Data.Header.Parameter,1),1);
    Data.Header.Parameter=[Data.Header.Parameter' fillup']';
    end;
    if ndims(Data.Value)==1 Data.Header.Size=[size(Data.Value,1) 1 1]; end;
    if ndims(Data.Value)==2 Data.Header.Size=[size(Data.Value,1) size(Data.Value,2) 1]; end;
    if ndims(Data.Value)==3 Data.Header.Size=[size(Data.Value,1) size(Data.Value,2) size(Data.Value,3)]; end;
    Data.Header.Parameter(14)=Data.Header.Tiltangle(1);
end;

if isstruct(Data)~=1
    if isa(Data,'double') | isa(Data,'single')
            MODE=2;
    end;
    if isa(Data,'int16') 
          MODE=1;
    end;
    if isa(Data,'int8') 
            MODE=0;
    end;
    comment=char(zeros(800,1));
    parameter=zeros(52,1);
    if ndims(Data)==1 image_size=[size(Data,1) 1 1]; end;
    if ndims(Data)==2 image_size=[size(Data,1) size(Data,2) 1]; end;
    if ndims(Data)==3 image_size=[size(Data,1) size(Data,2) size(Data,3)]; end;
    Header=struct('Size',image_size','Comment',comment,'Parameter',parameter);
    Data=struct('Value',Data,'Header',Header);
    disp('Default MRC-header was created.');
end;

% open the stream always with the big endian format !
if isequal(format,'le');
    fid = fopen(em_name,'w','ieee-le');
else
    fid = fopen(em_name,'w','ieee-be');
end;

if fid==-1
    error(['Cannot open: ' em_name ' file']); 
end;
    fwrite(fid,Data.Header.Size(1),'long');
    fwrite(fid,Data.Header.Size(2),'long');
    fwrite(fid,Data.Header.Size(3),'long');
    fwrite(fid,MODE,'long');
    fwrite(fid,Data.Header.Parameter(1:52),'long');
    fwrite(fid,Data.Header.Comment(1:800),'char');
    xdim = Data.Header.Size(1);
    ydim = Data.Header.Size(2);
    zdim = Data.Header.Size(3);

for lauf=1:zdim
	
	Data_write=Data.Value(1:xdim,1:ydim,lauf);

	if MODE==0
		fwrite(fid,Data_write,'int8');
	elseif MODE==1
		fwrite(fid,Data_write,'int16');
	elseif MODE==2
		fwrite(fid,Data_write,'float');
	else
		disp('Sorry, i cannot write this as an EM-File !!!');
	end;
end;


fclose(fid);
