function motivelist = av3_createmotl(corrvolume, anglevol, nparticles, radius, ...
    phi_end, psi_end, theta_end, angular_incr, phi_start, psi_start, theta_start)
%AV3_COLLECTPARTICLES 
%   motivelist = av3_createmotl(corrvolume, 'anglevol', nparticles, radius,  ...
%       phi_end, psi_end, theta_end, angular_incr, phi_end, psi_end, theta_end)
%
%   Routine determines Peaks of a correlation volume (=CORRVOLUME) - after
%   determination of the peak the values in CORRVOLUME are set to zero
%   within RADIUS. 
%   A MOTIVELIST is created conatiang the cross-correlation values and best
%   orientations of the potential particles as given by the file ANGLEVOL.
%   In order to interpret this file, the angles PHI_END, PSI_END, and
%   THETA_END and the ANGULAR_INCrment have to be known. In case the
%   starting values xxx_START are non-zero, they have to be given, too.
%
%INPUT
%   corrvolume     : volume of cross correlation
%   'anglevol'     : [character] filename of volume containg angles as
%   indexes
%   nparticles     : number of desired peaks
%   radius         : radius of particles (ccf is set zero within radius
%                         after peak determination)
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
%   oscar, AV3_AVERAGE, AV3_COLLECTPARTICLES, TOM_CHOOSER
%
%   started 08/09/02 FF
%   last revision 11/27/02 FF

error(nargchk(8,11,nargin));
if (nargin < 9)
    phi_start = 0;
end;
if (nargin < 10)
    psi_start = 0;
end;
if (nargin < 11)
    theta_start = 0;
end;

motivelist = zeros(20, nparticles);
dlimx = radius+1; ulimx = size(corrvolume,1) - radius-1;    %set parameters for borders
dlimy = radius+1; ulimy = size(corrvolume,2) - radius-1;
dlimz = radius+1; ulimz = size(corrvolume,3) - radius-1;
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
        angleindex = tom_emreadc(anglevol, 'subregion', [x y z],[0 0 0]);
        [phi, psi, theta] = av3_index2angle(double(angleindex.Value), phi_end, psi_end, theta_end, angular_incr, ...
            phi_start, psi_start, theta_start);
        disp([num2str(i) '. Peak angles   = ' num2str(phi) ' ' num2str(psi) ' ' num2str(theta) ])
        motivelist(17,i) = phi;
        motivelist(18,i) = psi;
        motivelist(19,i) = theta;
        motivelist(4,i) = i;
        i = i + 1;
    end;
end;
motivelist = motivelist(:,1:i-1);

