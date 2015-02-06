function [Data_read] = tom_rawread(em_name, form, endian, dims, byte, header)
%TOM_RAWREAD reads data in RAW-file format
%
%   [Data_read] = tom_rawread(em_name,form,endian,dims,byte,header)
%
%PARAMETERS
%
%  INPUT
%   em_name             filename
%   form                'char', 'short, 'int16', 'int32', 'float', 'double'
%   endian              'le' for little-endian(PC), 'be' for big-endian(SGI,MAC)
%   dims                dimensions
%   byte                block length (currently not used)
%   header              header skippes in bytes
%  
%  OUTPUT
%   Data_read			Image or 3D data
%
%EXAMPLE
%   i=tom_rawread;
%       calls tom_rawreadgui the user interface and reads the
%       file on the variable i.
%
%   i=tom_rawread('recon_data_final.dat','float','le',[200 200 100],0,512);
%       skips 512 bytes header and reads data in float format, big endian type
%       with a dimension of 200x200x100 voxels.
%
%   to import an Amria segmentation mask save it as .raw (3D) and read it:
%   i=tom_rawread('labels_final.raw','char','le',[512 512 64],0,0);
%
%REFERENCES
%
%SEE ALSO
%   TOM_EMWRITE, TOM_EMHEADER, TOM_READEMHEADER, TOM_RAWREADGUI
%
%   created by SN 10/23/02
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

if nargin==0 Data_read=tom_rawreadgui;
else

% open the stream with the correct format !
if isequal(endian,'be')
    fid = fopen(em_name,'r','ieee-be');
else
    fid = fopen(em_name,'r','ieee-le');
end;

if fid==-1
    error(['Cannot open: ' em_name ' file']); 
end;
fseek(fid,header,-1);
xdim = dims(1);
ydim = dims(2);
zdim = dims(3);

% which byte format?
% to adapt to EM, transpose
% internal loop is faster ! (hopefully!)
Data_read=zeros(xdim,ydim,zdim);

    if isequal('char',form)
    	for lauf=1:zdim
    		Data_read(:,:,lauf) = fread(fid,[xdim,ydim],'char');
    	end;
    elseif isequal('short',form)
    	for lauf=1:zdim
    		Data_read(:,:,lauf) = fread(fid,[xdim,ydim],'short');
    	end;
    elseif isequal('uint8',form)
    	for lauf=1:zdim
    		Data_read(:,:,lauf) = fread(fid,[xdim,ydim],'uint8');
    	end;
    elseif isequal('long',form)
    	for lauf=1:zdim
		    Data_read(:,:,lauf) = fread(fid,[xdim,ydim],'long');
	    end;
    elseif isequal('uint16',form)
    	for lauf=1:zdim
		    Data_read(:,:,lauf) = fread(fid,[xdim,ydim],'uint16');
	    end;
    elseif isequal('int16',form)
    	for lauf=1:zdim
		    Data_read(:,:,lauf) = fread(fid,[xdim,ydim],'int16');
	    end;
    elseif isequal('ushort',form)
    	for lauf=1:zdim
		    Data_read(:,:,lauf) = fread(fid,[xdim,ydim],'ushort');
	    end;
    elseif isequal('int32',form)
    	for lauf=1:zdim
		    Data_read(:,:,lauf) = fread(fid,[xdim,ydim],'int32');
	    end;
    elseif isequal('uint32',form)
    	for lauf=1:zdim
		    Data_read(:,:,lauf) = fread(fid,[xdim,ydim],'uint32');
	    end;
    elseif isequal('float',form)
	    for lauf=1:zdim
    		Data_read(:,:,lauf) = fread(fid,[xdim,ydim],'float');
        end;
    elseif isequal('double',form)
	    for lauf=1:zdim
    		Data_read(:,:,lauf) = fread(fid,[xdim,ydim],'double');
        end;
    else
    disp('Sorry, i cannot read this as an RAW-File !!!');
    Data_read=[];
    end;
fclose(fid);
end;

