function Align = tom_av3_create_alignlist(filenames, Align, waitbarflag)
%TOM_AV3_CREATE_ALIGNLIST creates or appends alignment list
%
%   Align = tom_av3_create_alignlist(filenames, Align, waitbarflag)
%
%PARAMETERS
%
%  INPUT
%   filenames           2D char array with the full pathnames of the input volumes,
%                        can be generated with strvcat
%   Align               alignment structure to append (use [] to create a new
%                        structure)
%   waitbarflag         1 = waitbar gui
%                       0 = text messages
%  
%  OUTPUT
%   Align               ...
%
%EXAMPLE
%   ... = tom_av3_create_alignlist(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV3_PASTE4OSCARGUI, TOM_AV3_OSCARGUI, TOM_AV3_PASTE4OSCAR,
%   TOM_AV3_EXTRACT_ANGLESSHIFTS
%
%   created by AK 10/10/05
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


if nargin < 3
    waitbarflag = 0;
end

if waitbarflag == 1
    h = waitbar(0,'Creating alignment list');
end

for j = 1:size(filenames,1)
    if ~exist(strtrim(filenames(j,:)),'file')
        error(['File ' strtrim(filenames(j,:)) ' does not exist!']);
    end;
end;

if isempty(Align)
    run = 1;
else 
    run = size(Align,1)+1;
end

for i = 1:size(filenames,1)
 
    header = tom_reademheader(strtrim(filenames(i,:)));
    %Evaluate comment, format: Position.x Position.y Position.z offset binning TomogramFilename
    comment = char(header.Header.Comment');
    token = zeros(1,5);
    try
    if size(strtrim(comment),1) > 1
        for j = 1:5
                [tok, comment] = strtok(comment);
                token(j) = str2num(tok);
        end
    end

    end;
    
    
    Align(run,i).Filename = strcat(strtrim(filenames(i,:)));
    Align(run,i).Tomogram.Filename = strtrim(comment);
    Align(run,i).Tomogram.Header = header.Header;
    Align(run,i).Tomogram.Position.X = token(1); %Position of particle in Tomogram (values are unbinned)
    Align(run,i).Tomogram.Position.Y = token(2);
    Align(run,i).Tomogram.Position.Z = token(3);
    Align(run,i).Tomogram.Regfile = '';
    Align(run,i).Tomogram.Offset = token(4);     %Offset from Tomogram
    Align(run,i).Tomogram.Binning = token(5);    %Binning of Tomogram
    Align(run,i).Tomogram.AngleMin = header.Header.Tiltangle; 
    Align(run,i).Tomogram.AngleMax = header.Header.Tiltaxis;
    Align(run,i).Shift.X = 0; %Shift of particle, will be filled by tom_av3_extract_anglesshifts
    Align(run,i).Shift.Y = 0;
    Align(run,i).Shift.Z = 0;
    Align(run,i).Angle.Phi = 0; %Rotational angles of particle, will be filled by tom_av3_extract_anglesshifts
    Align(run,i).Angle.Psi = 0;
    Align(run,i).Angle.Theta = 0;
    Align(run,i).Angle.Rotmatrix = []; %Rotation matrix filled up with function tom_align_sum, not needed otherwise
    Align(run,i).CCC = 0; % cross correlation coefficient of particle, will be filled by tom_av3_extract_anglesshifts
    Align(run,i).Class = 0;
    Align(run,i).ProjectionClass = 0;
    Align(run,i).NormFlag = 0; %is particle phase normalized?
    Align(run,i).Filter = [0 0]; %is particle filtered with bandpass?

    
    if waitbarflag == 1
        waitbar(i./size(filenames,1),h,[num2str(i), ' of ', num2str(size(filenames,1)), ' files done.']);
    end
    
end




if waitbarflag == 1
    close(h);
end