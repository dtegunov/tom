function tom_vol2chimera(vol,vol2,vol3)
% tom_vol2chimera opens a 3D volume in Chimera, UCSF.
% 
% Input:
%           vol:                        input 3D volume.
%           vol2:                       second input 3D volume.
% Output:
%           opens the 3D volume(s) in chimera 
%
% tom_vol2chimera creates a temporary EM-file in the 'tempdir' of Matlab.
%
% Example:
% tom_vol2chimera(vol_1.Value);
%
% tom_vol2chimera(vol_1.Value,vol_2.Value);
%
% References:
% http://www.cgl.ucsf.edu/chimera/
%
%SEE ALSO
%   tom_chimera_volumes
%
%   created by SN 07/05/10
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2010
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

[a b]=unix('whoami');
b=b(1:end-1);


if (isnumeric(vol))
    filename=[tempdir b 'tom_tmp2chimera.em'];
    tom_emwritec(filename,vol);
    type1=' tom_em:';
end;

if (iscell(vol))    
    type1=' ';
    filename=[tempdir b 'tom_tmp2chimera1.bild'];
    write_sel(filename,vol);
end;

if (isstruct(vol))    
    type1=' ';
    filename=[tempdir b 'tom_tmp2chimera1.pdb'];
    tom_pdbwrite(filename,vol);
end;


if nargin>1
    filename2=[tempdir b 'tom_tmp2chimera2.em'];
    if (isnumeric(vol2))
        tom_emwritec(filename2,vol2);
        type2=' tom_em:';
    end;
    if (iscell(vol2))
        type2=' ';
        filename2=[tempdir b 'tom_tmp2chimera2.bild'];
        write_sel(filename2,vol2);
    end;
    if (isstruct(vol2))
        type2=' ';
        filename2=[tempdir b 'tom_tmp2chimera2.pdb'];
        tom_pdbwrite(filename2,vol2);
    end;
end;

if nargin==3
    filename3=[tempdir b 'tom_tmp2chimera3.em'];
    if (isnumeric(vol3))
        tom_emwritec(filename3,vol3);
        type3=' tom_em:';
    end;
    if (iscell(vol3))
        type3=' ';
        filename3=[tempdir b 'tom_tmp2chimera3.bild'];
        write_sel(filename3,vol3);
    end;
    if (isstruct(vol3))
        type3=' ';
        filename3=[tempdir b 'tom_tmp2chimera3.pdb'];
        tom_pdbwrite(filename3,vol3);
    end;
    
end;


if isunix
    [a chimera_path]=unix(['which chimera']);
    if isequal(a,0)
        if nargin==1
            unix(['chimera ' type1 '' filename ' &']);
            disp(['chimera ' type1 '' filename ' &']);
        end;
        if nargin==2
            unix(['chimera ' type1 '' filename '' type2 '' filename2 ' &']);
            disp(['chimera ' type1 '' filename '' type2 '' filename2 ' &']);
        end;
        if nargin==3
            unix(['chimera ' type1 '' filename ' ' type2 '' filename2 '' type3 '' filename3 ' &']);
            disp(['chimera ' type1 '' filename ' ' type2 '' filename2 '' type3 '' filename3 ' &']);
        end;
        
        disp(['open volume in chimera, path=' chimera_path]);
    else
        disp('cannot find chimera. Add it to the standard path');
    end;
    
end;
if ispc
    try
    chimera_ver=winqueryreg('HKEY_LOCAL_MACHINE','SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Chimera_is1','DisplayName');
    chimera_path=winqueryreg('HKEY_LOCAL_MACHINE','SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Chimera_is1','Inno Setup: App Path');
    chimera_path=['"' chimera_path '\bin\chimera.exe"'];
    disp(['open volume in ' chimera_ver ' at ' chimera_path]);
    catch
        chimera_path=['"' 'C:\Program Files (x86)\Chimera 1.5.3\bin\chimera.exe"'];
    end
    %    system(['"C:\Program Files\Chimera_1_4_1b\bin\chimera.exe"' ' --title ' filename ' tom_em:' filename ' &']);
    if nargin==2
        system([chimera_path ' --title ' filename ' tom_em:' filename ' tom_em:' filename2 ' &']);
    else
        system([chimera_path ' --title ' filename ' tom_em:' filename ' &']);
    end;
end;

function write_sel(filename,tmp_sel)


 fid=fopen(filename,'wt');
 for i=1:length(tmp_sel)
    fprintf(fid,'%s \n',tmp_sel{i});
 end;
 fclose(fid);



