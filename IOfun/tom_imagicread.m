function out = tom_imagicread(filename)
%TOM_IMAGICREAD reads imagic data files
%
%   out = tom_imagicread(filename)
%
%PARAMETERS
%
%  INPUT
%   filename     pathname and filename of the header file or image file
%
%  OUTPUT
%   out          image structure
%
%REFERENCES
%
%   http://www.imagescience.de/formats/formats.htm
%
%SEE ALSO
%   TOM_IMAGICWRITE, TOM_IMAGIC_CREATEHEADER
%
%   created by AK 08/12/06
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


if nargin < 1
    error('No filename given');
end

if length(filename) > 4
    if isequal(filename(end-3:end),'.hed') || isequal(filename(end-3:end),'.img')
        filename = filename(1:end-4);
    end
end

filename_body = [filename '.img'];
filename_header = [filename '.hed'];

fid_body = fopen(filename_body,'rb','ieee-le');
fid_header = fopen(filename_header,'rb','ieee-le');
if fid_body==-1
    error(['Cannot open: ' filename_body ' file for reading.']); 
end

if fid_header==-1
    error(['Cannot open: ' filename_header ' file for reading.']); 
end

%find out how many headers are stored in header file
fseek(fid_header,0,'eof');
header_size = ftell(fid_header);
fseek(fid_header,0,'bof');
if mod(header_size,1024) ~= 0
    error('Header file has invalid length');
end

numimages = header_size ./ 1024;

%read in every image
for i=1:numimages
    
    %read in header
    header = fread(fid_header,14,'int32');
    out.Header.IMAGIC(i).IMN = header(1);
    out.Header.IMAGIC(i).IFOL = header(2);
    out.Header.IMAGIC(i).IERROR = header(3);
    out.Header.IMAGIC(i).NHFR = header(4);
    out.Header.IMAGIC(i).NMONTH = header(5);
    out.Header.IMAGIC(i).NDAY = header(6);
    out.Header.IMAGIC(i).NYEAR = header(7);
    out.Header.IMAGIC(i).NHOUR = header(8);
    out.Header.IMAGIC(i).NMINUT = header(9);
    out.Header.IMAGIC(i).NSEC = header(10);
    out.Header.IMAGIC(i).NPIX2 = header(11);
    out.Header.IMAGIC(i).NPIXEL = header(12);
    out.Header.IMAGIC(i).IXLP = header(13);
    out.Header.IMAGIC(i).IYLP = header(14);
    header = fread(fid_header,4,'char');    
    out.Header.IMAGIC(i).TYPE = char(header');
    header = fread(fid_header,2,'int32');
    out.Header.IMAGIC(i).IXOLD = header(1);
    out.Header.IMAGIC(i).IYOLD = header(2);
    header = fread(fid_header,6,'float');
    out.Header.IMAGIC(i).AVDENS = header(1);
    out.Header.IMAGIC(i).SIGMA = header(2);
    out.Header.IMAGIC(i).VARIAN = header(3);
    out.Header.IMAGIC(i).OLDAVD = header(4);
    out.Header.IMAGIC(i).DENSMAX = header(5);
    out.Header.IMAGIC(i).DENSMIN = header(6);
    header = fread(fid_header,1,'int32');
    out.Header.IMAGIC(i).COMPLEX = header(1);    
    header = fread(fid_header,5,'float');
    out.Header.IMAGIC(i).CXLENGTH = header(1);
    out.Header.IMAGIC(i).CYLENGTH = header(2);
    out.Header.IMAGIC(i).CZLENGTH = header(3);
    out.Header.IMAGIC(i).CALPHA = header(4);
    out.Header.IMAGIC(i).CBETA = header(5);    
    header = fread(fid_header,80,'char');
    out.Header.IMAGIC(i).NAME = char(header');    
    header = fread(fid_header,1,'float');
    out.Header.IMAGIC(i).CGAMMA = header(1);    
    header = fread(fid_header,14,'int32');
    out.Header.IMAGIC(i).MAPC = header(1);
    out.Header.IMAGIC(i).MAPR = header(2);
    out.Header.IMAGIC(i).MAPS = header(3);
    out.Header.IMAGIC(i).ISPG = header(4);
    out.Header.IMAGIC(i).NXSTART = header(5);
    out.Header.IMAGIC(i).NYSTART = header(6);
    out.Header.IMAGIC(i).NZSTART = header(7);
    out.Header.IMAGIC(i).NXINTV = header(8);
    out.Header.IMAGIC(i).NYINTV = header(9);
    out.Header.IMAGIC(i).NZINTV = header(10);
    out.Header.IMAGIC(i).IZLP = header(11);
    out.Header.IMAGIC(i).I4LP = header(12);
    out.Header.IMAGIC(i).I5LP = header(13);
    out.Header.IMAGIC(i).I6LP = header(14);
    header = fread(fid_header,3,'float');
    out.Header.IMAGIC(i).ALPHA = header(1);
    out.Header.IMAGIC(i).BETA = header(2);
    out.Header.IMAGIC(i).GAMMA = header(3);    
    header = fread(fid_header,32,'int32');
    out.Header.IMAGIC(i).IMAVERS = header(1);
    out.Header.IMAGIC(i).REALTYPE = header(2);
    out.Header.IMAGIC(i).RONLY = header(32);
    header = fread(fid_header,2,'float');    
    out.Header.IMAGIC(i).ANGLE = header(1);
    out.Header.IMAGIC(i).RCP = header(2);
    header = fread(fid_header,2,'int32');    
    out.Header.IMAGIC(i).IXPEAK = header(1);
    out.Header.IMAGIC(i).IYPEAK = header(2);
    header = fread(fid_header,3,'float');    
    out.Header.IMAGIC(i).CCC = header(1);
    out.Header.IMAGIC(i).ERRAR = header(2);
    out.Header.IMAGIC(i).ERR3D = header(3);
    header = fread(fid_header,1,'int32');    
    out.Header.IMAGIC(i).REF = header(1);
    header = fread(fid_header,13,'float');    
    out.Header.IMAGIC(i).CLASSNO = header(1);    
    out.Header.IMAGIC(i).LOCOLD = header(2);
    out.Header.IMAGIC(i).OLDAVD = header(3);
    out.Header.IMAGIC(i).OLDSIGMA = header(4);    
    out.Header.IMAGIC(i).XSHIFT = header(5);
    out.Header.IMAGIC(i).YSHIFT = header(6);
    out.Header.IMAGIC(i).NUMCLS = header(7);
    out.Header.IMAGIC(i).OVQVAL = header(8);
    out.Header.IMAGIC(i).EANGLE = header(9);
    out.Header.IMAGIC(i).EXSHIFT = header(10);
    out.Header.IMAGIC(i).EYSHIFT = header(11);
    out.Header.IMAGIC(i).CMTOTVAR = header(12);
    out.Header.IMAGIC(i).INFORMAT = header(13);
    header = fread(fid_header,2,'int32');    
    out.Header.IMAGIC(i).NUMEIGEN = header(1);    
    out.Header.IMAGIC(i).NIACTIVE = header(2);    
    header = fread(fid_header,8,'float');    
    out.Header.IMAGIC(i).RESOLX = header(1);    
    out.Header.IMAGIC(i).RESOLY = header(2);    
    out.Header.IMAGIC(i).RESOLZ = header(3);    
    out.Header.IMAGIC(i).ALPHA2 = header(4);    
    out.Header.IMAGIC(i).BETA2 = header(5);    
    out.Header.IMAGIC(i).GAMMA2 = header(6);    
    out.Header.IMAGIC(i).NMETRIC = header(7);        
    out.Header.IMAGIC(i).ACTMSA = header(8);        
    fseek(fid_header,69*4,'cof');
    header = fread(fid_header,228,'char');
    out.Header.IMAGIC(i).HISTORY = char(header');        
    
    out.Header.Size = [out.Header.IMAGIC(i).IXLP, out.Header.IMAGIC(i).IYLP, numimages];
    out.Header.Objectpixelsize = out.Header.IMAGIC(i).RESOLX;
    %read in image
    if i==1
        switch out.Header.IMAGIC(i).TYPE
            case 'REAL'
                type = 'float';
                t2 = 'single';
            case 'INTG'
                type = 'int16';
                t2 = 'int16';
            case 'PACK'
                type = 'uint8';
                t2 = 'uint8';
            case 'COMP'
                type = 'double';
                t2 = 'double';
            otherwise
                error(['Format ' out.Header.IMAGIC(i).TYPE ' not supported.']);
        end
    end
    
    %if out.Header.Size(3) < 2
        if i==1
            out.Value = zeros(out.Header.Size(1),out.Header.Size(2),numimages,t2);
        end
        out.Value(:,:,i) = reshape(fread(fid_body,out.Header.Size(1).*out.Header.Size(2),type),out.Header.Size(1),out.Header.Size(2));
    %else
    %    if i==1
    %        out.Value = zeros(out.Header.Size(1),out.Header.Size(2),out.Header.Size(3),numimages);
    %    end
    %    out.Value(:,:,:,i) = fread(fid_body,out.Header.Size(1).*out.Header.Size(2).*out.Header.Size(3),type);
    %end
    
end

fclose(fid_body);
fclose(fid_header);

