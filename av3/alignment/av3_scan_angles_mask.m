function av3_scan_angles_mask(refilename, motlfilename, particlefilename, startindx, iterations,angincr, ...
    angiter,semiangle,mask,iclass,hipass,lowpass,nfold,rcorr)
% AV3_SCANANGLES aligns 3D subtomograms to reference
%
%   av3_scan_angles_mask(refilename, motlfilename, particlefilename, startindx, iterations,angincr,...
%   angiter,semiangle,mask,iclass,hipass,lowpass,nfold,rcorr)
%   Filenames are expected as:
%       'filename'_#no.em
%
%   In this routine, the reference is rotated into the proposed orientation
%   (given in the motl). A X-correlation with the particle under scrutiny is
%   performed for this orientation and the ANGITER*ANGINCR neighbouring
%   orienations. 
%   A clever sampling on the unit sphere is done for psi and theta: rings around 
%   (theta_old, psi_old) are drawn. 
%   The routine takes into account the missing wedge - the semi angle of
%   the wedge SEMIANGLE has to be given. ROI is the radius if interest
%   which is used for masking the particles and RSMOOTH is the smoothing of
%   this sphere (additional to ROI). Only particles of the class ICLASS
%   (column 20 of MOTL) are correlated. For determination of CCF-peaks
%   spline interpolation is used - sub-pixel accuracy.
%
%  PARAMETERS
%   refilename          filename of reference(s) - 'refilename'_#no.em
%   motlfilename        filename of corresponding motl - 'motlfilename'_#no.em
%   particlefilename    filename of 3D-particles to be aligned and averaged
%                           'particlefilename'_#no.em
%   startindx           start index - Index of first reference AND motl
%   iterations          number of iterations to be performed -after each
%                           iteration motl and average are stored
%   angincr             angular increment
%   angiter             iterations for each angle phi, psi, theta
%   semiangle           semi angle of missing wedge
%   mask                mask - make sure dims are all right!
%   iclass              class of particles
%   hipass              hipass - for X-corr
%   lowpass             lowpass - for X-corr
%   nfold               symmetry of particle
%   rcorr               radius for masking X-Corr function
%
%   09/17/03 FF
%    last change 11/05/04 FF

npixels = sum(sum(sum(mask)));
ind = startindx;
name = [refilename '_' num2str(ind) '.em'];
ref = tom_emread(name);
wedge=tom_wedge(ref.Value,semiangle);
for ind = startindx:startindx+iterations
    name = [refilename '_' num2str(ind) '.em'];
    ref = tom_emread(name);
    disp(['read in file ' name]);
    ref = ref.Value;
    ref = tom_symref(ref,nfold);
    [mref xx1 xx2 mstd] = tom_dev(ref,'noinfo');
    ref = (ref - mref)./mstd;
    ref = tom_limit(ref,-3,4,'z'); % throw away the gold
    [mref xx1 xx2 mstd] = tom_dev(ref,'noinfo');
    ref = (ref - mref)./mstd;
    average=ref*0;
    wei = zeros(size(average,1),size(average,2),size(average,3));%weighting function
    %mask = tom_spheremask(ones(size(average,1),size(average,2),size(average,3)),roi,rsmooth);
    name = [motlfilename '_' num2str(ind) '.em'];
    motl = tom_emread(name);
    motl = motl.Value;
    for indpart = 1:size(motl,2)
        if (motl(20,indpart) == iclass)
            tshift = 0;
            phi_old=motl(17,indpart);
            psi_old=motl(18,indpart);
            the_old=motl(19,indpart);
            % read shift BEFORE rot
            xshift = motl(14,indpart);
            yshift = motl(15,indpart);
            zshift = motl(16,indpart);ccc = -1;
            ifile = motl(4,indpart);
            name = [particlefilename '_' num2str(ifile) '.em'];
            particle = tom_emread(name);particle = particle.Value;
            [meanv, maxv, minv, stdv] = tom_dev(particle,'noinfo');
            particle = (particle-meanv)/stdv;
            particle = tom_limit(particle,-3,4,'z'); % throw away the gold
            particle4av = particle;
            % do not shift particle but mask !
            rshift = [xshift yshift zshift];
            rmask = double(tom_rotate(mask,[phi_old,psi_old,the_old]));
            % rshift = tom_pointrotate(rshift,[phi_old,psi_old,the_old]);
            % %changed 09/04/03 FF - now compatible with Reiner!
            shiftmask = tom_shift(rmask,rshift);
            particle=shiftmask.*particle;
            particle= particle - shiftmask.*(sum(sum(sum(particle)))/npixels);%subtract mean in sphere
            fpart=fftshift(tom_fourier(particle));
            %apply bandpass
            fpart= ifftshift(tom_spheremask(fpart,lowpass,3) - tom_spheremask(fpart,hipass,2));
            %normalize
            fpart(1,1,1)=0;
            fpart = fpart/sqrt((sum(sum(sum(fpart.*conj(fpart))))));
            for phi = phi_old-angiter*angincr:angincr:phi_old+angiter*angincr
                for ithe =  0:ceil(angiter/2)
                    if ithe == 0
                        npsi=1;
                        dpsi=360;
                    else
                        %sampling for psi and the on unit sphere in rings
                        dpsi=angincr/sin(ithe*angincr/180*pi);
                        npsi = ceil(360/dpsi);
                    end;
                    for ipsi = 0:(npsi-1)
                        r = [ 0 0 1];
                        r = tom_pointrotate(r,0,ipsi*dpsi,ithe*angincr);
                        r = tom_pointrotate(r,0,psi_old,the_old);
                        the = 180/pi*atan2(sqrt(r(1).^2+r(2).^2),r(3));
                        psi = 180/pi*atan2(r(2),r(1)) + 90;
                        rpart=double(tom_rotate(ref,[phi,psi,the]));
                        rpart=tom_ifourier(ifftshift(tom_spheremask(wedge.*fftshift(tom_fourier(rpart)))));
                        rmask=double(tom_rotate(mask,[phi,psi,the]));
                        rpart = rpart.*rmask;%mask with smoothened edges
                        rpart = rpart - rmask.*(sum(sum(sum(rmask.*rpart)))/npixels); %subtract mean in mask
                        fref=fftshift(tom_fourier(rpart));
                        %apply bandpass
                        fref=ifftshift(tom_spheremask(fref,lowpass,3) - tom_spheremask(fref,hipass,2));
                        fref(1,1,1)=0;
                        %calculate rms - IMPORTANT! - missing wedge!!!
                        fref = fref/sqrt((sum(sum(sum(fref.*conj(fref))))));
                        ccf = npixels*tom_spheremask(real(fftshift(tom_ifourier(fpart.*conj(fref)))),rcorr,1);%changed 09/17/03
                        [pos ccctmp] = peak_det_2(real(ccf));
                        if ccctmp > ccc
                            ccc = ccctmp;
                            phi_opt=phi;
                            psi_opt=psi;
                            the_opt=the;
                            tshift = pos-33;
                            tom_dspcub(ccf)
                        end;
                    end;
                end;
            end;
            motl(17,indpart)=phi_opt;
            motl(18,indpart)=psi_opt;
            motl(19,indpart)=the_opt;
            motl(11,indpart) = tshift(1);
            motl(12,indpart) = tshift(2);
            motl(13,indpart) = tshift(3);
            %rshift = tom_pointrotate(rshift,-psi_opt,-phi_opt,-the_opt);
            motl(14,indpart) = tshift(1);%changed 09/04/03
            motl(15,indpart) = tshift(2);
            motl(16,indpart) = rshift(3);
            motl(1,indpart)=ccc;
            % take care: particle4av is NOT pre-shifted
            average = average + double(tom_rotate(tom_shift(particle4av,-tshift),[-psi_opt,-phi_opt,-the_opt]));
            tom_dspcub((average));
            %weighting - avoid interpolation artefacts
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi_opt,-phi_opt,-the_opt])),0.5,1,'z'),0,0.5);
            wei = wei + tmpwei;
            disp(['Particle no ' num2str(indpart) ' , Iteration no ' num2str(ind)]);
        end; %endif 
    end;% end particle loop
    name = [motlfilename '_' num2str(ind+1) '.em'];
    tom_emwrite(name,motl);
    % do weighting
    lowp = floor(size(average,1)/2)-3;%lowpass
    wei = 1./wei;rind = find(wei > 100000);wei(rind) = 0;% take care for inf
    average = real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(average)).*wei,lowp))));
    name = [refilename '_' num2str(ind+1) '.em'];
    tom_emwrite(name,average);
    disp(['wrote reference ' name]);
end; % end iteration loop
