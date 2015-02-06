function [shiftvec, angles, ccc, ccf_max, phivol, psivol, thevol] = ...
    av3_dock_phase_only(ref,particle,hipass,lowpass,angincr,startangles,endangles,flag)
% AV3_DOCK_PHASE_ONLY aligns 3D reference (e.g. from X-ray) to particle (from EM)
%
%   [shiftvec, angles, ccc, ccf_max, phi_opt, psi_opt, the_opt] = ...
%       av3_dock_phase_only(ref,particle,hipass,lowpass,angincr, ...
%       startangles,endangles,flag)
%
%   The function works just like AV3_DOCK, apart from the fact that the
%   phase-only filter is used instead of the CCF.
%
%  PARAMETERS
%   REF             reference - e.g. generated from X-ray structure
%   PARTICLE        particle - e.g. EM map to dock REF in
%   LOWPASS         lowpass for data prior to X-corr
%   HIPASS          hipass for data prior to X-corr
%   ANGINCR         angular increment (in deg)
%   STARTANGLES     start angles - vector [phi psi theta] in deg
%   ENDANGLES       end angles - vector [phi psi theta] in deg
%   
%   SHIFTVEC        vector to shift ref
%   ANGLES          angles to turn ref - vector [phi psi theta] in deg
%   CCC             maximum X-corr coefficient
%   CCF_MAX         X-correlation function max volume
%   PHI_OPT         optimal phi volume
%   PSI_OPT         optimal psi volume
%   THE_OPT         optimal the volume
%
%
%  SEE ALSO
%   AV3_DOCK, AV3_DOCK-NONORM
%
%   12/12/03 FF
%
%
cent= [floor(size(ref,1)/2)+1 floor(size(ref,2)/2)+1 floor(size(ref,3)/2)+1];
scf = sqrt(size(ref,1)*size(ref,2)*size(ref,3));
npixels = size(ref,1)*size(ref,2)*size(ref,3);
%[mref xx1 xx2 mstd] = tom_dev(ref,'noinfo');
%ref = (ref - mref)./mstd;
tshift = 0;
error(nargchk(5,8,nargin))
if (nargin < 8)
    flag = 'normal';
end;
if (nargin < 7)
    phiend = 360 - angincr;
    psiend = 360 - angincr;
    theend = 180;
else
    phiend = endangles(1);
    psiend = endangles(2);
    theend = endangles(3);
end;
if (nargin < 6)
    phistart = 0;
    psistart = 0;
    thestart = 0;
else
    phistart = startangles(1);
    psistart = startangles(2);
    thestart = startangles(3);
end;
nsteps = (phiend-phistart+angincr)/angincr *(psiend-psistart+angincr)/angincr * (theend-thestart+angincr)/angincr;
ccf_max = -1000*ones(size(ref,1),size(ref,2),size(ref,3));
phivol = -ones(size(ref,1),size(ref,2),size(ref,3));
psivol = -ones(size(ref,1),size(ref,2),size(ref,3));
thevol = -ones(size(ref,1),size(ref,2),size(ref,3));
fpart = tom_fourier(particle);
if (nargin == 8)
    if (strmatch(flag,'mutual','exact') == 1);
        disp('using mutual correlation function')
        fpart = fpart./sqrt(abs(fpart));
        fref = tom_fourier(ref);
        ref = real(tom_ifourier(fref./sqrt(abs(fref))));
    elseif (strmatch(flag,'phase-only','exact') == 1);
        disp('using phase-only correlation function')
        fpart = fpart./(abs(fpart));
        fref = tom_fourier(ref);
        ref = real(tom_ifourier(fref./(abs(fref))));
    end;
end;
%apply bandpass - check if necessary...
fpart= ifftshift(tom_spheremask(fftshift(fpart),lowpass,1) - tom_spheremask(fftshift(fpart),hipass,1));
icount = 0;
flag = 0;
disp([' ' num2str(flag*5) ' % done']);
for phi = phistart:angincr:phiend
    for psi = psistart:angincr:psiend
        for the = thestart:angincr:theend
            rref = double(tom_rotate(ref,[phi,psi,the]));
            %apply wedge and bandpass %calculate rms - IMPORTANT! - missing wedge!!!
            fref = fftshift(tom_fourier(rref));
            fref = ifftshift(tom_spheremask(fref,lowpass,1) - tom_spheremask(fref,hipass,1));
            ccf = real(fftshift(tom_ifourier(fpart.*conj(fref))))/npixels;
            ind = find(ccf > ccf_max);
            ccf_max(ind) = ccf(ind);
            phivol(ind) = phi;
            psivol(ind) = psi;
            thevol(ind) = the;
            icount = icount+1;
            tmp = double(int16(icount*20/nsteps));
            if tmp > flag
                flag = tmp;
                disp([' ' num2str(flag*5) ' % done']);
            end;
        end;
    end;
end;
[pos ccc] = peak_det_2(ccf_max);
phi_opt=phivol(round(pos(1)),round(pos(2)),round(pos(3)));
psi_opt=psivol(round(pos(1)),round(pos(2)),round(pos(3)));
the_opt=thevol(round(pos(1)),round(pos(2)),round(pos(3)));
shiftvec = pos-cent;
if size(ccf)>100
    %tom_dspcub(tom_bin(ccf_max));
else
    %tom_dspcub(ccf);
end;
angles = [phi_opt psi_opt the_opt];

