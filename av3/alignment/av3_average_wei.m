function [average, average_wei] = av3_average_wei(motivelist, filename, r, flag, class,semiangle)
% AV3_AVERAGE_WEI performs WEIGTED averaging of 3D-particles from entire
%           tomograms 
%
% USAGE
%   [average, average_wei] = av3_average_wei(motivelist, filename, r, flag, class,semiangle);
%
%   The principle is the following: A MOTIVELIST contains features of
%   interest and proposed orientations. The locations of the particles are
%   stored in columns 8-10. The corresponding tomogram is stored as
%   FILENAME. The features are boxed out around a radius of interest R,
%   rotated into the suggested orientation (MOTIVELIST column 17-19) and
%   averaged. The FLAG set to 'var' subtracts the mean value from every 
%   individual tomogram and normalizes it according to its variance prior 
%   to averaging. The CLASS of the features is specified in column 20 of
%   the MOTIVELIST. The SEMIANGLE has to be chosen according to the
%   tilt-series and is used to create a proper weighting function for the
%   final average.
%
% PARAMETERS
%   MOTIVELIST    motivelist (motl); for format see AV3_COLLECTPARTICLES 
%   FILENAME      string; filename of tomogram
%                      -> files have to be EM-format
%   R             Radius of particles
%   FLAG          - 'var'    normalize volumes to variance before averaging
%                 - 'normal' not normalized
%   CLASS         class of particles (column 20 in MOTL)
%   SEMIANGLE     semiangle of missing wedge - e.g. 30, if tilt range
%                   -60 deg to 60 deg
%
% OUTPUT
%   AVERAGE       averaged particle
%   AVERAGE_WEI   weighted average, according to missing wedge of data
%
% SEE ALSO
%   AV3_COLLECTPARTICLES, AV3_AVERAGE
%
%   FF 08/18/03
%   last change
%   FF 02/03/04
error(nargchk(6,6,nargin));
npart=size(motivelist,2);
ipart = 0;
volheader = tom_reademheader(filename);
for i=1:npart,
    if (motivelist(20,i) == class )
        %ifile = motivelist(4,i); %corresponding filenumber
        %file=strcat(filename, '_', num2str(ifile), '.em');
        x=motivelist(8,i);y=motivelist(9,i);z=motivelist(10,i);
        if ((x > r+2) & (y>r+2) & (z>r+2) & ...
                (x<volheader.Header.Size(1)-r-1) & (y<volheader.Header.Size(2)-r-1) ...
                & (z<volheader.Header.Size(3)-r-1)) 
            part = tom_emread(filename,'subregion',[x-r-1 y-r-1 z-r-1], [2*r+1 2*r+1 2*r+1]);
            phi = motivelist(17,i); psi = motivelist(18,i); theta=motivelist(19,i);
            %add translation later! to be done!
            part = double(tom_rotate(part.Value,[ -psi, -phi, -theta]));
            % normalize if flag is set
            if strmatch(flag,'var')
                [meanp dummy dummy rms ] = tom_dev(part,'noinfo');
                part = (part-meanp)/rms;
            end;
            % average and calculate weighting
            if ipart == 0
                wedge=tom_wedge(part,semiangle);
                wei = 0*part;%initialize psf and average
                average = part*0; 
            end;
            average = part + average;
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi,-phi,-theta])),0.5,1,'z'),0,0.5);
            wei = wei + tmpwei;% calculate pointspread function in fourier space
            ipart = ipart + 1; % number of particles averaged
        end;
        disp(['processed line ' num2str(i) ' of MOTL']);
        if (ipart-floor(ipart/10)*10 < 1)
            %tom_dspcub(average);
        end;
    end;
end;
disp(['  ' num2str(ipart) '  particles averaged '  ])
file=strcat(filename, '_average.em');
tom_emwrite(file, average);
lowp = floor(size(average,1)/2)-3;
wei = 1./wei;rind = find(wei > 100000);wei(rind) = 0;% take care for inf
average_wei = real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(average)).*wei,lowp))));
file=strcat(filename, '_average_wei.em');
tom_emwrite(file, average_wei);
figure;tom_dspcub(average);
figure;tom_dspcub(average_wei);
