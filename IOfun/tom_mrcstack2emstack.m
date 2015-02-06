function tom_mrcstack2emstack(source_mrc_file, new_em_stack_name, extension, format)
%TOM_MRCSTACK2EMSTACK converts a MRC stack into a EM stack
%
%   tom_mrcstack2emstack(source_mrc_file,new_em_stack_name,extension,format)
%
% tom_mrcstack2emseries converts a MRC(FEI style) stack into an Em-Stack
% by reading the header and then the images
%
%PARAMETERS
%
%  INPUT
%   source_mrc_file     mrc input file
%   new_em_stack_name   name of the output file
%   extension           extension of the output files
%   format              machine-byte-order  set format 'le' for little endian files (PC,Linux) 
%                       or 'be' for big endian (SGI,Mac)
%  
%  OUTPUT
%
%EXAMPLE
%   tom_mrcstack2emstack('A_03.mrc','Vol','em','le');
%
%REFERENCES
%
%SEE ALSO
%   tom_mrc2emseries, tom_emread, tom_emwrite
%
%   created by SN 08/01/04
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
    fid = fopen(source_mrc_file,'r','ieee-le');
else
    fid = fopen(source_mrc_file,'r','ieee-be');
end;

if fid==-1
    error(['Cannot open: ' source_mrc_file ' file']); 
end;
Parameter.nx = fread(fid,[1],'int32');
Parameter.ny = fread(fid,[1],'int32');
Parameter.nz = fread(fid,[1],'int32');
Parameter.mode = fread(fid,[1],'int32')
    Parameter.nxstart= fread(fid,[1],'int32');
    Parameter.nystart= fread(fid,[1],'int32');
    Parameter.nzstart= fread(fid,[1],'int32');
    Parameter.mx= fread(fid,[1],'int32');
    Parameter.my= fread(fid,[1],'int32');
    Parameter.mz= fread(fid,[1],'int32');               % : integer;
    Parameter.xlen= fread(fid,[1],'float');
    Parameter.ylen= fread(fid,[1],'float');
    Parameter.zlen= fread(fid,[1],'float');
    Parameter.alpha= fread(fid,[1],'float');
    Parameter.beta= fread(fid,[1],'float');
    Parameter.gamma= fread(fid,[1],'float');           % : single;
    Parameter.mapc= fread(fid,[1],'int32'); 
    Parameter.mapr= fread(fid,[1],'int32'); 
    Parameter.maps= fread(fid,[1],'int32');         %  : integer;
    Parameter.amin= fread(fid,[1],'float'); 
    Parameter.amax= fread(fid,[1],'float'); 
    Parameter.amean= fread(fid,[1],'float');        %  : single;
    Parameter.ispg= fread(fid,[1],'int32'); 
    Parameter.nsymbt = fread(fid,[1],'int32');      % : integer;
    Parameter.unk = fread(fid,[16],'int32');        %   : array[0..15] of integer; { this is unused junk }
    Parameter.idtype= fread(fid,[1],'int32');      %   : integer;
    Parameter.nd1= fread(fid,[1],'int16'); 
    Parameter.nd2= fread(fid,[1],'int16'); 
    Parameter.vd1= fread(fid,[1],'int16'); 
    Parameter.vd2= fread(fid,[1],'int16');        %   : smallint;
    Parameter.tiltangles= fread(fid,[9],'float');  % : array[0..8] of single;
    Parameter.zorg= fread(fid,[1],'float'); 
    Parameter.xorg= fread(fid,[1],'float'); 
    Parameter.yorg = fread(fid,[1],'float');       %   : single;
    Parameter.nlabl = fread(fid,[1],'int32');     %    : integer;
    Parameter.labl = fread(fid,[800],'char');     %   : array[0..9] of p80 end of std MRC format, rest is FEI special;
%  TMrcExtendedHeader = record

Data_read=zeros(Parameter.nx,Parameter.ny,1);

for lauf=1:Parameter.nz
    Extended(lauf).a_tilt= fread(fid,[1],'float');
    Extended(lauf).b_tilt= fread(fid,[1],'float');
    Extended(lauf).x_stage= fread(fid,[1],'float');
    Extended(lauf).y_stage=fread(fid,[1],'float');
    Extended(lauf).z_stage=fread(fid,[1],'float');
    Extended(lauf).x_shift=fread(fid,[1],'float');
    Extended(lauf).y_shift=fread(fid,[1],'float');
    Extended(lauf).defocus=fread(fid,[1],'float');
    Extended(lauf).exp_time=fread(fid,[1],'float');
    Extended(lauf).mean_int=fread(fid,[1],'float');
    Extended(lauf).tiltaxis=fread(fid,[1],'float');
    Extended(lauf).pixelsize=fread(fid,[1],'float').*10^9;
    Extended(lauf).magnification=fread(fid,[1],'float'); % : single
    
    fseek(fid,128-52,0);
end;

fseek(fid,128*(1024-Parameter.nz),0);

% create one time a large EM stack file
tom_emwritec([new_em_stack_name '.' extension],[Parameter.nx Parameter.ny Parameter.nz],'new');

for lauf2=1:Parameter.nz
    lauf2
    if Parameter.mode==0
        Data_read(:,:,1) = (fread(fid,[Parameter.nx,Parameter.ny],'int8'));
    elseif Parameter.mode==1
        Data_read(:,:,1) = (fread(fid,[Parameter.nx,Parameter.ny],'int16'));
    elseif Parameter.mode==2
        Data_read(:,:,1) = (fread(fid,[Parameter.nx,Parameter.ny],'float'));
    else
        error(['Sorry, i cannot read this as an MRC-File !!!']);
        Data_read=[];
    end;
    lauf2
        %write slice by slice in EM format
        tom_emwritec([new_em_stack_name '.' extension], Data_read, 'subregion',[1 1 lauf2] ,[Parameter.nx-1 Parameter.ny-1 1]);
    lauf2
end
fclose(fid);
clear Data_read;

