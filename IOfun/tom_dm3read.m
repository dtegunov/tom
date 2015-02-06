function [m sx sy]=tom_dm3read(filename,logfile)
% This file was adapted from a routine by Fred Sigworth ReadDM3.m by SN,
% 2010.
%
% function [m sx sy]=tom_dm3read(filename,logfile)
% Read a Digital Micrograph file and return the first image in its native
% data format (e.g. int16), along with its pixel scale information
% (typically in nm).  If a logfile name is specified, a description of the
% entire tree will be written to the console and also to the file.

% F. Sigworth, July 2009
% Code was based on the description by Greg Jefferis at
% http://rsb.info.nih.gov/ij/plugins/DM3Format.gj.html
% 
% This function has been written assuming that it will run on a
% little-endian (e.g. Intel) machine, reading a file written by a
% little-endan machine. Otherwise the use of the swapbytes function will
% have to be made conditional on the local machine type and the byteorder
% variable. Also, swapbytes may have to be added to the GetData function.
% 
% Note that the code here allows any
% fields from the tree to be extracted. Here is where we define the fields
% that we will extract.  We grab the data value at the first match of each
% of these tags.  Here numerals represent unnamed fields.  To see what all
% the tags are, specify a logfile to receive the hundreds of lines of
% information!
%
% Copyright (c) 2010, Fred Sigworth
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

celltags={'ImageList 2 ImageData Calibrations Dimension 1 Scale'
    'ImageList 2 ImageData Calibrations Dimension 2 Scale'
    'ImageList 2 ImageData Dimensions 1'
    'ImageList 2 ImageData Dimensions 2'
    'ImageList 2 ImageData Data'};

found=zeros(size(celltags));
output=cell(size(celltags));

% Set up the log file.  We use my mprintf function to allow printing to the
% console as well, and suppressing printing when the handle is zero.
if nargin>1
    flog=fopen(logfile,'w');
    hs=[1 flog];  % log file handles
else
    flog=0;
    hs=[0];
end;
tabstring='| ';
level=0;
maxprint=4;
OutputOn=1;

% Read the whole file into memory as a byte array.
fb=fopen(filename,'r');
d=fread(fb,inf,'*uint8');  % read the whole file as bytes
fclose(fb);

p=int32(1);  % byte pointer--also a global variable
Tags=cell(1,10); % Keeps track of the sequence of tags, for matching with the tag strings.


% Pick up the header
version=GetLong;
if version ~=3
    error(['tom_dm3read: Wrong file type.  Version = ' num2str(version)]);
end;
nbytes=GetLong;

% Handle little- and big-endian files and machines
dle=GetLong;  % boolean to tell whether data is little endian
[str,maxsize,endian] = computer;
mle= (endian=='L');  % machine is little endian: we'll have to swap bytes in reading the tree.
dswap=(dle~=mle);  % swap byte-order when reading data

% Traverse the tree
GetTagGroup;

% Close the logfile
if flog>0
    fclose(flog);
end;

% Extract the output parameters
sx=output{1};
sy=output{2};
xdim=output{3};
ydim=output{4};
m=reshape(output{5},xdim,ydim);

% end of main function

% ---- here are all the local functions, called recursively ----

    function GetTagGroup
        sorted=GetByte;
        open=GetByte;
        NumTags=GetLong;
        for i=1:NumTags
            GetTagEntry(i);
        end;
    end

    function GetTagEntry(MemberIndex)
        level=level+1;
        PutNew;
        PutTabs;
        isdata=GetByte;
        labelsize=GetInt;
        labelstring=GetString(labelsize);
        PutStr('-');
        PutStr([labelstring ':']);
        if numel(labelstring)<1
            labelstring=num2str(MemberIndex);
        end;
        Tags{level}=labelstring;
        if isdata==21
            GetTagType
        elseif isdata==20
            GetTagGroup
        else
            error(['Unknown TagEntry type ' num2str(isdata)]);
        end;
        Tags{level}=[];
        level=level-1;
    end

    function GetTagType
        dum=GetLong;
        if dum ~= 623191333
            disp(['Illegal TagType value ' num2str(dum)]);
        end
        deflen=GetLong;
        EncType=GetLong;
        x=GetData(EncType);
        index=CheckTags;
        if index>0
            output{index}=x;
        end;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This is the function that slows everything down, which checks the
    % entire tag list for each tagged type.  It is inefficient, but simple.
    function r=CheckTags
        for i=1:numel(celltags)
            ok=~found(i);
            c=celltags{i};
            j=1;
            while ok && (numel(c)>0)
                [s c]=strtok(c);
                ok=strcmp(s,'*')||strcmp(s,Tags{j});
                j=j+1;
            end;
            if ok
                r=i;
                return
            end;
        end;
        r=0;
    end


    function x=GetData(ftype,num)
        if nargin<2
            num=1;
        end;
        x=[];
        %         disp(['GetData ' num2str(ftype)]);
        switch ftype
            case 2  % short
                x=typecast(d(p:p+num*2-1),'int16');
                p=p+2*num;
            case 3  % long
                x=typecast(d(p:p+num*4-1),'int32');
                p=p+4*num;
            case 4  % ushort
                x=typecast(d(p:p+num*2-1),'uint16');
                p=p+2*num;
            case 5  % ulong
                x=typecast(d(p:p+num*4-1),'uint32');
                p=p+4*num;
            case 6  % float
                x=typecast(d(p:p+num*4-1),'single');
                p=p+4*num;
            case 7  % double
                x=typecast(d(p:p+num*8-1),'double');
                p=p+8*num;
            case 8  % boolean
                x=d(p:p+num-1);
                p=p+num;
            case 9  % char
                x=char(d(p:p+num-1));
                p=p+num;
            case 10  % octet
                x=(d(p:p+num-1));
                p=p+num;
            case 15  % Struct
                PutStr('struct');
                StructNameLength=GetLong;
                NumFields=GetLong;
                x=[];
                for i=1:NumFields
                    FieldNameLength(i)=GetLong;
                    FieldType(i)=GetLong;
                end;
                StructName=GetString(StructNameLength);
                PutStr(StructName);
                PutNew; PutTabs;
                for i=1:NumFields
                    %                     FieldNameLen=FieldNameLength(i);
                    FieldName=GetString(FieldNameLength(i));
                    FieldTy=FieldType(i);
                    PutStr(FieldName);
                    x(i)=GetData(FieldType(i));
                    PutNew; PutTabs;
                end;
            case 18 % string
                length=GetLong;
                x=char(d(p:p+length-1)');
                PutVal(x); PutNew;
                p=p+length;
                
            case 20  % Array
                ArrayType=GetLong;
                if ArrayType==15  % Struct is special case
                    StructNameLength=GetLong;
                    NumFields=GetLong;
                    x=[];
                    for i=1:NumFields
                        FieldNameLength(i)=GetLong;
                        FieldType(i)=GetLong;
                    end;
                end;
                ArrayLength=GetLong;
                
                if ArrayType ~=4
                    PutStr('array of');
                    PutVal(ArrayLength);
                    PutStr(' --type'); PutVal(ArrayType);
                end;
                
                if ArrayType==15
                    PutStr('structs');
                    PutNew;
                    for j=1:ArrayLength
                        OutputOn=j<=maxprint;
                        for i=1:NumFields
                            FieldNameLen=FieldNameLength(i);
                            FieldName=GetString(FieldNameLength(i));
                            FieldTy=FieldType(i);
                            PutTabs;
                            PutStr(FieldName);
                            x(i)=GetData(FieldType(i));
                            PutNew;
                        end;
                        OutputOn=1;
                    end;
                elseif ArrayType==4
                    OutputOn=0;
                    for j=1:ArrayLength
                        x(j)=GetData(ArrayType);
                    end;
                    OutputOn=1;
                    PutVal(char(x'));
                else
                    % Might be long data
                    if (ArrayLength > 1000)  % try to handle a long array
                        OutputOn=0;
                        x=GetData(ArrayType,ArrayLength);
                        OutputOn=1;
                    else
                        PutNew;
                        for j=1:ArrayLength
                            OutputOn=j<=maxprint;
                            PutTabs;
                            x(j)=GetData(ArrayType);
                            PutNew;
                        end;
                        OutputOn=1;
                    end; % long data
                end;
            otherwise
                x=0;
                disp(['Unrecognized data type ' num2str(ftype)]);
        end; % switch
        if (ftype < 15) && OutputOn
            PutVal(x);
        end;
    end % GetData


    function PutStr(s)
        if OutputOn
            mprintf(hs,'%s ',s);
        end;
    end

    function PutTabs
        if OutputOn
            for i=1:level-1
                mprintf(hs,'%s',tabstring);
            end;
        end;
    end

    function PutVal(x)
        if OutputOn
            if isa(x,'char')
                mprintf(hs,'%s ',x);
            else
                mprintf(hs,'%d ',x);
            end;
        end;
    end

    function PutNew
        if OutputOn
            mprintf(hs,'\n');
        end;
    end

    function s=GetString(len)
        len=int32(len);
        if len<1
            s=[];
        else
            s=char(d(p:p+len-1)');
            p=p+len;
        end;
    end

    function x=GetLong
        
        x=typecast(d(p:p+3),'int32');
        x=swapbytes(x);
        p=p+4;
        
    end

    function x=GetWord
        x=typecast(d(p:p+1),'int16');
        x=swapbytes(x);
        p=p+2;
    end

    function x=GetInt
        x=typecast(d(p:p+1),'int16');
        x=swapbytes(x);
        p=p+2;
    end
    function x=GetByte
        x=d(p);
        p=p+1;
    end

    function mprintf(handles,varargin)
    % function mprintf(handles,varargin) % copy of my utility function to
    % make ReadDM3 self-contained.
    % Write the same formatted text to multiple files.  handles is an array of
    % file handles.  The function fprintf is called multiple times, once for
    % each handle number.  Handles of 0 are ignored.
    for i=1:numel(handles)
        if handles(i)>0
            fprintf(handles(i),varargin{:});
        end;
    end
    end % mprintf
end