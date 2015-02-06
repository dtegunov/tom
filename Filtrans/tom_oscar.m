function [ccfmax, phi_max, psi_max, theta_max] = tom_oscar(vol,template,phi_start, phi_end, phi_incr, ...
        psi_start, psi_end, psi_incr, the_start, the_end, the_incr)
%TOM_OSCAR performs 3D cross correlation of template and volume.
%
%[ccfmax, phi_max, psi_max, theta_max] = tom_oscar(vol,templ,phi_start, phi_end, phi_incr, psi_start, psi_end, psi_incr, the_start, the_end, the_incr)
%
%PARAMETERS
%
%  INPUT
%   vol                 volume (e.g. sub-tomogram)
%   templ               template to be used for cross correlation
%   phi_start           start value of phi (in deg)
%   phi_end             end value of phi (in deg)
%   phi_incr            angular increment of phi (in deg)
%   psi_start           start value of psi (in deg)
%   psi_end             end value of psi (in deg)
%   psi_incr            angular increment of psi (in deg)
%   the_start           start value of start (in deg)
%   the_end             end value of start (in deg)
%   the_incr            angular increment of start (in deg)
%  
%  OUTPUT
%   ccfmax              CCF map with maximum correlation (respective to angles)
%   phi_max             map with corresponding phi
%   psi_max             map with corresponding psi
%   theta_max           map with corresponding the
%
%   TOM_OSCAR performs 3D cross correlation of template and volume.
%   Exhaustive six-dimensional search is performed - takes long...
%   Only maxima of CCF wih respect to the orientation are stored.
%   NO normalization is implemented yet!
%
%EXAMPLE
%   vol = tom_emread('Ribosome_frank_1.7nm.em');
%   [ccfmax, phi_max, psi_max, theta_max] = tom_oscar(vol.Value,vol.Value,0, 10, 5, 0, 10, 5, 0, 10, 5);
%
%REFERENCES
%
%SEE ALSO
%   TOM_ORCD, TOM_CORR, TOM_CCC
%
%   created by FF 06/26/03
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


[nvolx nvoly nvolz] = size(vol);
[ntemplx ntemply ntemplz] = size(template);
offsetx = nvolx/2+1-ntemplx/2;
offsety = nvoly/2+1-ntemply/2;
offsetz = nvolz/2+1-ntemplz/2;
%fourier transform volume and conjugate
vol = conj(fftn(vol));
%initialize outputs
ccfmax = -1000*ones(nvolx, nvoly, nvolz);
phi_max = int16(-ones(nvolx, nvoly, nvolz));
psi_max = int16(-ones(nvolx, nvoly, nvolz));
theta_max = int16(-ones(nvolx, nvoly, nvolz));
for phi=phi_start:phi_incr:phi_end
    for psi=psi_start:psi_incr:psi_end
        for the=the_start:the_incr:the_end
            temp = zeros(nvolx, nvoly, nvolz);
            temp(offsetx:offsetx+ntemplx-1,offsety:offsetx+ntemply-1,offsetz:offsetz+ntemplz-1)= ...
                double(tom_rotate(template,[phi,psi,the]));
            temp = real(ifftshift(ifftn(fftn(temp).*vol)));%perform convolution
            indx = find(temp > ccfmax);
            ccfmax(indx) = temp(indx);
            phi_max(indx) = int16(phi);
            psi_max(indx) = int16(psi);
            theta_max(indx) = int16(the);
        end;
    end;
    disp([num2str(((phi-phi_start)/phi_incr+1)/((phi_end-phi_start)*phi_incr+1)*100) '%']);
end;
            
