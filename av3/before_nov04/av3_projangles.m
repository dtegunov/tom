function anglelist = av3_projangles(motl,alignmentfile,class,z0)
% AV3_PROJANGLES calculates projection angles of individual projections 
%   of particles located in an entire tomogram. 
%   
%   anglelist = av3_projangles(motl,alignmentfile,class,irefmark,r0)
%
%   The routine is needed for
%   exact weighted reconstruction of averaged particles. Attention: so far,
%   only for default centering of backprojection, generalization to be
%   done. 
%   
%
%   
% INPUT
%   motl            : motivelist (array) 
%   alignmentfile   : file [char] 
%   class           : class of particles (default: 0)
%   irefmark        : Index of reference marker
%   z0              : z-origin from rec (corresponding to motl) - default
%                       1025
%
% OUTPUT
%   anglelist:
%       array of dim  npart x ntilt x 6
%       1: index particle (ipart) -> from MOTL - but only belonging to class!
%       2: index of tilt (itilt) -> from alignmentfile - corresponds to
%               running index of micrograph
%       3: phi, psi, theta, gamma_4baproj, x, y -> phi, theta projection
%               vector; psi rot of micrograph; gamma_4baproj is the
%               rotation angle that has to be applied on the individual
%               micrographs; x and y are centers of particles on individual
%               micrographs.
%   
%   see J.Walz, PhD-Thesis p. 170 f.
%   
%   04/14/03 FF

error(nargchk(2,4,nargin));
alig = tom_emread(alignmentfile);
alig = alig.Value;
[alig beta_tilt] = tom_alignment3d(alig); % beta_tilt = image rotation 
deltax =alig(5,:,1);
deltay = alig(6,:,1);
if nargin < 3
    class=0;
end;
if nargin < 4
    z0 = 1025;
    x0=1025;
    y0=1025;
end;
%   define Rotation Matrix after Projection (in EM)
R = [tom_cos(-beta_tilt) -tom_sin(-beta_tilt); tom_sin(-beta_tilt) tom_cos(-beta_tilt)];
irun = 0;
theta_p=alig(1,:,1); % tiltangle
psi_p(1:size(alig,2))=0; % psi - normally assumed to be zero
for ipart =1: size(motl,2)
    iclass = motl(20,ipart);
    if (iclass == class)
        irun = irun +1;
        theta_r = motl(19,ipart);% read Euler Angles from MOTL!
        psi_r = motl(18,ipart);
        phi_r = motl(17,ipart);
        for itilt = 1:size(alig,2)
            %   1st: determine Projection angle theta_pprime and phi_prime
            %   o M_pprime total projection matrix for individual particle in
            %       micrograph
            %   o M_r Rotation Matrix according to MOTL in tomogram
            %   o M_p Projection matrix according to tomogram
            %   o resulting total projection:
            %       M_pprime = M_p * M_r   
            %   note: M_p consists of Proj P_p and Rotation R_p
            %       M_pprime = R_p * P_p * M_r   
            cos_theta_pprime = tom_sin(theta_p(itilt))*tom_sin(phi_p(itilt))*tom_sin(theta_r)*tom_sin(psi_r)-...
                tom_sin(theta_p(itilt))*tom_cos(phi_p(itilt))*tom_sin(theta_r)*tom_cos(psi_r) + tom_cos(phi_p(itilt))*tom_cos(theta_r);
            nomin = tom_sin(theta_p(itilt))*tom_sin(phi_p(itilt)) * ( tom_cos(psi_r)*tom_cos(phi_r)-tom_cos(theta_r)*tom_sin(psi_r)*tom_sin(phi_r) ) + ...
                tom_sin(theta_p(itilt))*tom_cos(phi_p(itilt)) * ( tom_cos(psi_r)*tom_cos(phi_r)+tom_cos(theta_r)*tom_cos(psi_r)*tom_sin(psi_r) );
            denum = -tom_sin(theta_p(itilt))*tom_sin(phi_p(itilt)) * ( tom_cos(psi_r)*tom_sin(phi_r)+tom_cos(theta_r)*tom_sin(psi_r)*tom_cos(phi_r) ) - ...
                tom_sin(theta_p(itilt))*tom_cos(phi_p(itilt)) * ( tom_sin(psi_r)*tom_sin(phi_r)-tom_cos(theta_r)*tom_cos(psi_r)*tom_cos(phi_r) );
            tan_phi_prime=nomin/denum;
            nomin = tom_cos(theta_p(itilt))*tom_sin(phi_p(itilt))*tom_sin(theta_r(itilt))*tom_sin(psi_r(itilt)) - ...
                tom_cos(theta_p(itilt))*tom_cos(phi_p(itilt))*tom_sin(theta_r)*tom_cos(psi_r) - tom_sin(theta_p(itilt))*tom_cos(theta_r);
            denum = tom_cos(phi_p(itilt))*tom_sin(theta_r)*tom_sin(psi_r) + tom_sin(phi_p(itilt))*tom_sin(theta_r)*tom_cos(psi_r);
            tan_psi_prime = nomin/denum;
                        anglelist(irun,itilt,2) = atan(tan_phi_prime)*180/pi;
            anglelist(irun,itilt,3) = atan(tan_psi_prime)*180/pi;
            anglelist(irun,itilt,1) = acos(cos_theta_pprime)*180/pi;
            %   2nd: positions of particles on micrographs
            x = motl(14,ipart);
            y = motl(15,ipart);
            z = motl(16,ipart);
            %   proj before rotation in EM
            x_pprime = deltax(itilt) + (x - x0)*tom_cos(theta_p(itilt)) - (z - z0)*tom_sin(theta_p(itilt)) ;
            y_pprime = deltay(itilt) +  (y-y0);
            %   take into account rot in EM
            xy_onproj = R*[x_pprime-orig  y_pprime-orig] + orig;
            gamma_4baproj = theta_p(itilt)+anglelist(irun,itilt,3)+90;
            anglelist(irun,itilt,4) = gamma_4baproj;
            anglelist(irun,itilt,5) = x_pprime;
            anglelist(irun,itilt,6) = y_pprime;
        end;
    end;
end;