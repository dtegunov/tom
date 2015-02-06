function motivelist = av3_collectparticles(volume, corrvolume, anglevol, nparticles, radius, boxdim, ...
    filename, phi_end, psi_end, theta_end, angular_incr, phi_start, psi_start, theta_start)
%AV3_COLLECTPARTICLES 
%   motivelist = av3_collectparticles(volume, corrvolume, anglevol, nparticles, radius, boxdim, ...
%       filename, phi_end, psi_end, theta_end, angular_incr, phi_end, psi_end, theta_end)
%
%   Routine cuts out subtomograms of a tomogram (=VOLUME). Peaks of a correlation volume (=CORRVOLUME) 
%   are determined - after determination of the peak the values in CORRVOLUME are set to zero within RADIUS. 
%   In VOLUME boxes of dimension BOXDIM corresponding to each peak in the tomograms are cut out 
%   and rotated as proposed by ANGLEVOL. NPARTICLES boxes are cut out and stored in the files
%   FILENAME_{counter}.em . If FILENAME == 'xxx' then no files are stored,
%   only the motivelist is created.
%
%INPUT
%   volume         : volume containing particles
%   corrvolume     : volume of cross correlation
%   anglevol       : volume containg angles as indexes
%   nparticles     : number of desired peaks
%   radius         : radius of particles (ccf is set zero within radius after peak determination)
%   boxdim         : dimension of boxes to cut out
%   filename[char] : name of files stored as <file>_num.em
%   phi_end        : parameter from cross correlation
%   psi_end        : s.a.
%   theta_end      : s.a.
%   angular_incr   : increment (same for ALL angles)
%
%OUTPUT
%   motivelist     : matrix containing parameters of particles:
%
%   The following parameters are stored in the matrix MOTIVELIST of dimension (20, NPARTICLES):
%   column 
%      1         : Cross-Correlation Coefficient
%      2         : x-coordinate in full tomogram
%      3         : y-coordinate in full tomogram
%      4         : peak number
%      8         : x-coordinate in full tomogram
%      9         : y-coordinate in full tomogram
%      10        : z-coordinate in full tomogram
%      14        : x-shift in subvolume
%      15        : y-shift in subvolume
%      16        : z-shift in subvolume
%      17        : Phi
%      18        : Psi
%      19        : Theta 
%      20        : class
%
%SEE ALSO
%   oscar, AV3_AVERAGE, TOM_CHOOSER
%
%   started 08/09/02 FF
%   last revision 08/14/02 FF

error(nargchk(11,14,nargin));
if (nargin < 12)
    phi_start = 0;
end;
if (nargin < 13)
    psi_start = 0;
end;
if (nargin < 14)
    theta_start = 0;
end;

motivelist = zeros(20, nparticles);
boxdim2 = floor(boxdim/2);
dlimx = boxdim2; ulimx = size(corrvolume,1) - boxdim2;    %set parameters for borders
dlimy = boxdim2; ulimy = size(corrvolume,2) - boxdim2;
dlimz = boxdim2; ulimz = size(corrvolume,3) - boxdim2;
i = 1;
for lauf = 1:nparticles,
    [r CCC corrvolume] = tom_peak(corrvolume, radius); 
    x = r(1); y=r(2); z=r(3);
    %   x = 18; y = 25; z = 30; CCC = 0.5;
    %   check whether x,y,z not too close to border
    if (x > dlimx) & (x < ulimx) & (y > dlimy) & (y < ulimy) &(z > dlimz) & (z < ulimz) 
        disp([num2str(i) '. Peak:......... = ' num2str(x) ' ' num2str(y) ' ' num2str(z) ])
        motivelist(2,i) = x; motivelist(8,i) = x;
        motivelist(3,i) = y; motivelist(9,i) = y;
        motivelist(10,i) = z;
        motivelist(1,i) = CCC;
        [phi, psi, theta] = av3_index2angle(anglevol(x,y,z), phi_end, psi_end, theta_end, angular_incr, ...
            phi_start, psi_start, theta_start);
        disp([num2str(i) '. Peak angles   = ' num2str(phi) ' ' num2str(psi) ' ' num2str(theta) ])
        motivelist(17,i) = phi;
        motivelist(18,i) = psi;
        motivelist(19,i) = theta;
        if  ~ strmatch(filename,'xxx')
            box = volume(x-boxdim2:x-boxdim2+boxdim-1,y-boxdim2:y-boxdim2+boxdim-1,z-boxdim2:z-boxdim2+boxdim-1);
            rotbox = tom_rotate(box, [-psi, -phi, -theta]);
            rotbox = tom_emheader(rotbox);
            name = [filename '_' num2str(i) '.em'];
            tom_emwrite( name, rotbox);
        end;
        i = i + 1;
    end;
end;
motivelist = motivelist(:,1:i-1);
temp = tom_emheader(motivelist);
name = ['MOTL' filename '_' '.em'];
tom_emwrite(name, temp);
