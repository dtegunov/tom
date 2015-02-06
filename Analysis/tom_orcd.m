function ccf = tom_orcd(vol,template,phi_start, phi_end, phi_incr, ...
        psi_start, psi_end, psi_incr, the_start, the_end, the_incr, flag)  
% TOM_ORCD performs three-dim. orientational search
%
%   ccc =tom_orcd(vol,template,phi_start, phi_end, phi_incr, psi_start, ...
%       psi_end, psi_incr, the_start, the_end, the_incr, flag) 
%
%   The second array is rotated with respect to the first one, thereby
%   scanning the given angular range for the best correlation. For the
%   definition of angles see TOM_ROTATE3D. The resulting set of cross
%   correlation coefficients  is  stored  at CCF and printed, thereby
%   representing phi as rows, psi as columns and theta as layers.  
%
%  INPUT
%   vol         volume
%   template    template is rotated with respect to volume
%   phi_start   start value of phi
%   phi_end     end value of phi
%   phi_incr    increment of phi
%   psi_start   start value of psi
%   psi_end     end value of psi
%   psi_incr    increment of psi
%   the_start   start value of theta
%   the_end     end value of theta
%   the_incr    increment of theta
%   flag        flag - if set to 'norm' normalized CCCs will be calculated
%
%  OUTPUT
%   ccc         cross-correlation coefficients
%
%EXAMPLE
%   vol = tom_emread('Ribosome_frank_1.7nm.em');
%   ccc =tom_orcd(vol.Value,vol.Value,0, 330, 30, 0, ...
%       330, 30, 0, 180, 30, 'norm') 
%
%REFERENCES
%
%   created by FF 06/24/03
%   updated by ...
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


if nargin < 12
    vol= vol-mean(mean(mean(vol)));
    vol = vol/sqrt(sum(sum(sum(vol.*vol))));
end;
iphi=1;ipsi=1;ithe=1;
for phi=phi_start:phi_incr:phi_end
    for psi=psi_start:psi_incr:psi_end
        for the=the_start:the_incr:the_end
            if nargin == 12
                ccf(iphi,ipsi,ithe) = tom_ccc(vol,double(tom_rotate(template,[phi,psi,the])),'norm');
            else
                ccf(iphi,ipsi,ithe) = tom_ccc(vol,double(tom_rotate(template,[phi,psi,the])));
            end;
            ithe=ithe+1;
        end;
        ithe=1;
        ipsi=ipsi+1;
    end;
    ipsi=1;
    iphi=iphi+1;
end;

function ccc = corr3(a,b)
%a= a-mean(mean(mean(a)));  % a already normalized before
b= b-mean(mean(mean(b)));
ccc=sum(sum(sum(a.*b)))/sqrt(sum(sum(sum(b.*b))));
