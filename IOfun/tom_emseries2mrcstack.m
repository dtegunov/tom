function tlt=tom_emseries2mrcstack(source_em_files,extension,from_nr,to_nr, new_mrc_stack)
%TOM_EMSERIES2MRCSTACK converts an EM series to MRC stack (FEI style)
%
%   tom_emseries2mrcstack(source_em_files,extension,from_nr,to_nr, new_mrc_stack)
%
%PARAMETERS
%
%  INPUT
%   source_em_files     source of name of EM files (e.g., 'MyFile_')
%   extension           extension of EM files (e.g., '.em')
%   from_nr             from number (=first index of EM files)
%   to_nr               to number (=last index of EM files)
%   new_mrc_stack       name of MRC stack
%  
%  OUTPUT
%
%EXAMPLE
%   tom_emseries2mrcstack('ROS_020205_','.em',1,47, 'ROS_020205a.mrc')
%
%REFERENCES
%
%SEE ALSO
%   tom_mrc2emstack, tom_emread, tom_emwrite
%
%   created by SN 04/19/05
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

fid = fopen(new_mrc_stack,'w','ieee-le');

if fid==-1
    error(['Cannot open: ' new_mrc_stack ' file']);
end;

em=tom_emread([source_em_files num2str(from_nr) extension]);
tlt = zeros(to_nr-from_nr+1,1);

fwrite(fid,em.Header.Size(1),'int32');
fwrite(fid,em.Header.Size(2),'int32');
fwrite(fid,to_nr-from_nr+1,'int32');
if em.Header.Magic(4)==5
    fwrite(fid,2,'int32');
end;
if em.Header.Magic(4)==2
    fwrite(fid,1,'int32');
end;

fwrite(fid,0,'int32');
fwrite(fid,0,'int32');
fwrite(fid,0,'int32');
fwrite(fid,em.Header.Size(1),'int32');
fwrite(fid,em.Header.Size(2),'int32');
fwrite(fid,to_nr-from_nr+1,'int32');
fwrite(fid,em.Header.Size(1),'float');
fwrite(fid,em.Header.Size(2),'float');
fwrite(fid,to_nr-from_nr+1,'float');
fwrite(fid,90,'float');
fwrite(fid,90,'float');
fwrite(fid,90,'float');
fwrite(fid,1,'int32');
fwrite(fid,2,'int32');
fwrite(fid,3,'int32');
fwrite(fid,min(min(em.Value)),'float');
fwrite(fid,max(max(em.Value)),'float');
fwrite(fid,mean(mean(em.Value)),'float');
fwrite(fid,0,'short');
fwrite(fid,0,'short');
fwrite(fid,[1024.*128],'long');
fwrite(fid,[0 0],'int32');      %   : integer;
fwrite(fid,0,'int16');
fwrite(fid,0,'int16');
fwrite(fid,0,'int16');
fwrite(fid,0,'int16');        %   : smallint;
fwrite(fid,zeros(9,1),'float');  % : array[0..8] of single;
fwrite(fid,0,'float');
fwrite(fid,0,'float');
fwrite(fid,0,'float');       %   : single;
fwrite(fid,0,'int32');     %    : integer;
fwrite(fid,zeros(800,1),'char');     %   : array[0..9] of p80 end of std MRC format, rest is FEI special;
%  TMrcExtendedHeader = record

% some extra bytes, debugged SN ???
fwrite(fid,zeros(1024-964,1),'char'); % fill up

for lauf=from_nr:to_nr
    em=tom_emread([source_em_files num2str(lauf) extension]);
    
    fwrite(fid,em.Header.Tiltangle,'float');
    fwrite(fid,em.Header.Tiltaxis,'float');
    tlt(lauf)=em.Header.Tiltangle;
    fwrite(fid,0,'float');
    fwrite(fid,0,'float');
    fwrite(fid,0,'float');
    fwrite(fid,0,'float');
    fwrite(fid,0,'float');
    fwrite(fid,em.Header.Defocus.*10e-11,'float');
    fwrite(fid,em.Header.Exposuretime,'float');
    fwrite(fid,mean2(em.Value),'float');
    fwrite(fid,em.Header.Tiltaxis,'float');
    %    em.Header.Objectpixelsize=1e-9; % 1 nm by default
    fwrite(fid,em.Header.Objectpixelsize.*1e-10,'float');
    fwrite(fid,em.Header.Magnification,'float');
    
    fwrite(fid,zeros(32-13,1),'float'); % fill up
    
    
    %    fwrite(fid,ones((1024-52)./4,1),'float');
    %    fseek(fid,128*1024-52,0);
    %    fseek(fid,128-52,0);
end;
fwrite(fid,zeros(128.*1024-((to_nr-from_nr+1).*128),1),'char'); % fill up

%fseek(fid,128*(1024-Parameter.nz),0);
for lauf=from_nr:to_nr
    em=tom_emread([source_em_files num2str(lauf) extension]);
    if em.Header.Magic(4)==2
        fwrite(fid,em.Value,'int16');
    elseif em.Header.Magic(4)==5
        fwrite(fid,em.Value,'float');
    else
        error('only int16 and float values!');
    end;
    
    %    fseek(fid,128-52,0);
end;

fclose(fid);

clear Data_read;

