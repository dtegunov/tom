function tom_imagicwrite(Data,type,filename_body,filename_header)
%TOM_IMAGICWRITE writes data as an imagic file which can be used in eman
%
%   tom_imagicwrite(Data,type,filename_body,filename_header)
%
%PARAMETERS
%
%  INPUT
%   Data                input matrix
%   type                must be '2dstack'
%   filename_body       pathname and filename of the image file
%   filename_header     pathname and filename of the header file
%
%  OUTPUT
%
%
%REFERENCES
%
%   http://www.imagescience.de/formats/formats.htm
%
%SEE ALSO
%   TOM_IMAGICREAD, TOM_IMAGIC_CREATEHEADER
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

fid_body = fopen(filename_body,'wb','ieee-le');
fid_header = fopen(filename_header,'wb','ieee-le');
if fid_body==-1
    error(['Cannot open: ' filename_body ' file for writing.']); 
end

if fid_header==-1
    error(['Cannot open: ' filename_header ' file for writing.']); 
end


if isequal(lower(type),'2dstack')
    
    out = tom_imagic_createheader(Data,'2dstack');
    
    
    for i=1:size(out.Header.IMAGIC,2)
        %write header
        fwrite(fid_header,out.Header.IMAGIC(i).IMN,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).IFOL,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).IERROR,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).NHFR,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).NMONTH,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).NDAY,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).NYEAR,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).NHOUR,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).NMINUT,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).NSEC,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).NPIX2,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).NPIXEL,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).IXLP,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).IYLP,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).TYPE,'char');        
        fwrite(fid_header,out.Header.IMAGIC(i).IXOLD,'int32');        
        fwrite(fid_header,out.Header.IMAGIC(i).IYOLD,'int32');                
        fwrite(fid_header,out.Header.IMAGIC(i).AVDENS,'float');        
        fwrite(fid_header,out.Header.IMAGIC(i).SIGMA,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).VARIAN,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).OLDAVD,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).DENSMAX,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).DENSMIN,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).COMPLEX,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).CXLENGTH,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).CYLENGTH,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).CZLENGTH,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).CALPHA,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).CBETA,'float');
        fwrite(fid_header,sprintf('%80s',out.Header.IMAGIC(i).NAME),'char');
        fwrite(fid_header,out.Header.IMAGIC(i).CGAMMA,'float');        
        fwrite(fid_header,out.Header.IMAGIC(i).MAPC,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).MAPR,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).MAPS,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).ISPG,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).NXSTART,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).NYSTART,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).NZSTART,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).NXINTV,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).NYINTV,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).NZINTV,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).IZLP,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).I4LP,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).I5LP,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).I6LP,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).ALPHA,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).BETA,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).GAMMA,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).IMAVERS,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).REALTYPE,'int32');
        fwrite(fid_header,zeros(1,29),'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).RONLY,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).ANGLE,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).RCP,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).IXPEAK,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).IYPEAK,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).CCC,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).ERRAR,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).ERR3D,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).REF,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).CLASSNO,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).LOCOLD,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).OLDAVD,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).OLDSIGMA,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).XSHIFT,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).YSHIFT,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).NUMCLS,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).OVQVAL,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).EANGLE,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).EXSHIFT,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).EYSHIFT,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).CMTOTVAR,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).INFORMAT,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).NUMEIGEN,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).NIACTIVE,'int32');
        fwrite(fid_header,out.Header.IMAGIC(i).RESOLX,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).RESOLY,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).RESOLZ,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).ALPHA2,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).BETA2,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).GAMMA2,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).NMETRIC,'float');
        fwrite(fid_header,out.Header.IMAGIC(i).ACTMSA,'float');
        fwrite(fid_header,zeros(1,69),'float');
        fwrite(fid_header,sprintf('%228s',out.Header.IMAGIC(i).HISTORY),'char');

        
        %write body
        if i==1
            switch out.Header.IMAGIC(i).TYPE
                case 'REAL'
                    type = 'float';
                case 'INTG'
                    type = 'int16';
                case 'PACK'
                    type = 'uint8';
                case 'COMP'
                    type = 'double';
                otherwise
                    error(['Format ' out.Header.IMAGIC(i).TYPE ' not supported.']);
            end
        end
        fwrite(fid_body,out.Value(:,:,i),type);
        
    end
    
    
    
else
    error('Operation not supported.');
end


fclose(fid_body);
fclose(fid_header);



