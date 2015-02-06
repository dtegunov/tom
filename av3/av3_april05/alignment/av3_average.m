function average = av3_average(motivelist, filename, r, flag, class)
%AV3_AVERAGE  performs averaging of 3D-particles from entire TOMOGRAM
%
%USAGE
%   average = av3_average(motivelist, filename, r, flag, class);
%
%   AV3_AVERAGE averages subtomograms (coordinates from MOTL file). For 
%   correct weighting, best use AV3_AVERAGE_WEI.
%
%PARAMETERS
%   MOTIVELIST    motivelist (motl); for format see AV3_TRANS-ROT-ALIG
%   FILENAME      string; filename of tomogram
%                      -> files have to be EM-format
%   R             Radius of particles
%   FLAG          - 'var'    normalize volumes to variance before averaging
%                 - 'normal' not normalized
%
%OUTPUT
%   AVERAGE       averaged particle
%
%SEE ALSO
%   AV3_COLLECTPARTICLES, AV3_AVERAGE
%
%   FF 09/24/02 
% last change 03/31/05 FF - changed docu
eps = 300;   % use as parameter for creating mask out of ps
error(nargchk(3,5,nargin));
if (nargin < 5) 
    class = 0;
end;
if (nargin < 4)
    flag = 'normal';
end;
npart=size(motivelist,2);
ipart = 0;
volheader = tom_reademheader(filename);
for i=1:npart,
    if (nargin < 4) | (motivelist(20,i) == class )
        ifile = motivelist(4,i); %corresponding filenumber
        %file=strcat(filename, '_', num2str(ifile), '.em');
        x=motivelist(8,i);y=motivelist(9,i);z=motivelist(10,i);
        if ((x > r+1) & (y>r+1) & (z>r+1) & ...
                (x<volheader.Header.Size(1)-r-1) & (y<volheader.Header.Size(2)-r-1) ...
                & (z<volheader.Header.Size(3)-r-1)) 
            part = tom_emread(filename,'subregion',[x-r y-r z-r], [2*r 2*r 2*r]);
            phi = motivelist(17,i); psi = motivelist(18,i); theta=motivelist(19,i);
            %add translation later! to be done!
            part = double(tom_rotate(part.Value, [-psi, -phi, -theta]));
            % normalize if flag is set
            if (nargin > 3)
                if strmatch(flag,'var')
                    [meanp dummy dummy rms ] = tom_dev(part);
                    part = (part-meanp)/rms;
                end;
            end;
            % average
            if i == 1
                average = part*0; %initialize
            end;
            average = part + average;
            ipart = ipart + 1; % number of particles averaged
        end;
    end;
end;
disp(['  ' num2str(ipart) '  particles averaged '  ])
file=strcat(filename, '_average.em');
tom_emwrite(file, average);