function tom_av3_copy_subvols(directory, targetdir, filename, fileextension)
%TOM_A3_COPY_SUBVOLS was used to copy particle volumes from Willy's subdirs, ... 
%not maintained, do not use
%
%   tom_av3_copy_subvols(directory, targetdir, filename, fileextension)
%
%PARAMETERS
%
%  INPUT
%   directory           ...
%   targetdir           ...
%   filename            ...
%   fileextension       ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av3_copy_subvols(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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
%


Align = struct();
count = 1;

if strcmp(fileextension(1), '.') ~= 1
    fileextension = strcat('.', fileextension);
end;

subdirs = dir(directory);

for i = 3:size(subdirs,1)
    
    if subdirs(i).isdir == 1

        files = dir(strcat(directory, '/', subdirs(i).name, '/sorted/particles/'));
        disp('');
        disp('------------------------------------------------------');
        disp(subdirs(i).name);
        disp('------------------------------------------------------');
        for j = 3:size(files,1)
            
            if files(j).isdir == 0
                disp(files(j).name);
                unix(['cp ', directory, '/', subdirs(i).name , '/sorted/particles/', files(j).name, ' /', targetdir, '/', filename, '_', num2str(count), fileextension]);
                Align(count).ID = strcat(filename, '_', num2str(i), fileextension); 
                Align(count).Source = strcat(subdirs(i).name);
                Align(count).Shift.X = 0;
                Align(count).Shift.Y = 0;
                Align(count).Shift.Z = 0;
                Align(count).Angle.Phi = 0;
                Align(count).Angle.Psi = 0;
                Align(count).Angle.Theta = 0;
                Align(count).CCF = 0;
                Align(count).Class = 0;
                Align(count).ProjectionClass = 0;

%                 picklist = tom_emread(strcat(directory, '/', subdirs(i).name, '/sorted/picklist.em'));
%                 Align(count).OriginalCoordinates.X = picklist.Value(6,j-2);
%                 Align(count).OriginalCoordinates.Y = picklist.Value(7,j-2);                
%                 Align(count).OriginalCoordinates.Z = picklist.Value(8,j-2);

                count2 = 1;
                fid = 0;
                ang = [];
                while fid ~= -1
                    fid = fopen(strcat(directory, '/', subdirs(i).name, '/sorted/AWR_notbin_', num2str(count2),'.em'));
                    if fid ~= -1
                        fclose(fid);
                        header = tom_reademheader(strcat(directory, '/', subdirs(i).name, '/sorted/AWR_notbin_',num2str(count2), '.em'));
                        ang = [ang, header.Header.Tiltangle];
                    end;
                    count2 = count2 + 1;
                end

                Align(count).Tiltseries.Max = max(ang);
                Align(count).Tiltseries.Min = min(ang);
                Align(count).Tiltseries.Angles = ang;
                
                count = count + 1;

            end;
            
        end;
        
    end;
    
end;
save(strcat('/',targetdir,'/alignlist.mat'), 'Align');