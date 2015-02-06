function av3_scan_phi(refilename, motlfilename, particlefilename, startindx, iterations,angincr, ...
    angiter,mask,hipass,lowpass,nfold,threshold,wedgelist,iclass, ibin)
% AV3_SCAN_PHI aligns 3D subtomograms to reference
%
%   av3_scan_phi(refilename, motlfilename, particlefilename, startindx, iterations,angincr,...
%   angiter,mask,hipass,lowpass,nfold,threshold,wedgelist,iclass, ibin)
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
%   iclass              class of particles - fixed to 1 or 2
%   hipass              hipass - for X-corr
%   lowpass             lowpass - for X-corr
%   nfold               symmetry of particle
%   threshold           Threshold*mean(ccc) is cutoff for averaging - e.g.
%                           choose 0.5
%   wedgelist           array containing the tilt range - 1st: no of
%                           tomogram, 2nd: minimum tilt, 3rd: max tilt
%   iclass              class of particles - default:0
%   ibin                binning - default 0
%
%
% Format of MOTL:
%    The following parameters are stored in the matrix MOTIVELIST of dimension (20, NPARTICLES):
%    column 
%       1         : Cross-Correlation Coefficient
%       2         : x-coordinate in full tomogram
%       3         : y-coordinate in full tomogram
%       4         : particle number
%       5         : running number of tomogram - used for wedgelist
%       6         : index of feature in tomogram (optional)
%       8         : x-coordinate in full tomogram
%       9         : y-coordinate in full tomogram
%       10        : z-coordinate in full tomogram
%       14        : x-shift in subvolume
%       15        : y-shift in subvolume
%       16        : z-shift in subvolume
%       17        : Phi (in deg)
%       18        : Psi
%       19        : Theta 
%       20        : class no
%
%   09/18/03 FF
%last change 03/07/04 - binning added

error(nargchk(13,15,nargin))
if nargin < 14 
    iclass = 0;
end;
if nargin < 15
    ibin = 0;
end;
if ibin > 0
    mask = tom_bin(mask,ibin);
end;
npixels = sum(sum(sum(mask)));
cent= [floor(size(mask,1)/2)+1 floor(size(mask,2)/2)+1 floor(size(mask,3)/2)+1];
scf = size(mask,1)*size(mask,2)*size(mask,3);
ind = startindx;
name = [refilename '_' num2str(ind) '.em'];
ref = tom_emread(name);
if ibin > 0
    ref = tom_bin(ref.Value,ibin);
end;
%wedge=tom_wedge(ref.Value,semiangle);
for ind = startindx:startindx+iterations-1
    name = [refilename '_' num2str(ind) '.em'];
    ref = tom_emread(name);
    disp(['read in file ' name]);
    ref = ref.Value;
    if ibin > 0
        ref = tom_bin(ref,ibin);
    end;
    ref = tom_symref(ref,nfold);
    [mref xx1 xx2 mstd] = tom_dev(ref,'noinfo');
    ref = (ref - mref)./mstd;
    %ref = tom_limit(ref,-3,4,'z'); % throw away the gold
    [mref xx1 xx2 mstd] = tom_dev(ref,'noinfo');
    ref = (ref - mref)./mstd;
    average=ref*0;
    wei = zeros(size(average,1),size(average,2),size(average,3));%weighting function
    %mask = tom_spheremask(ones(size(average,1),size(average,2),size(average,3)),roi,rsmooth);
    name = [motlfilename '_' num2str(ind) '.em'];
    motl = tom_emread(name);
    motl = motl.Value;
    if ibin > 0
        motl(11:16) = motl(11:16)/(2^ibin);%take binning into account
    end;
    indx = find (motl(1,:) > 0); meanv = mean(motl(1,indx));
    indx = find (motl(1,:) > threshold*meanv);
    itomo_old = 0;
    for indpart = 1:size(motl,2)
        if ((motl(20,indpart) == 1) | (motl(20,indpart) == 2) | (motl(20,indpart) == iclass))
            itomo = motl(5,indpart);
            if itomo_old ~= itomo %wedge stuff - exact weighting according to list
                xx = find(wedgelist(1,:) == itomo);
                minangle= wedgelist(2,xx);
                maxangle= wedgelist(3,xx);
                wedge = av3_wedge(ref,minangle,maxangle);
                itomo_old = itomo;
            end;
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
            if ibin > 0
                particle = tom_bin(particle,ibin);
            end;
            particle = tom_limit(particle,-3,4,'z'); % throw away the gold
            particle4av = particle;
            % do not shift particle but mask !
            rshift(1) = motl(11,indpart);
            rshift(2) = motl(12,indpart);
            rshift(3) = motl(13,indpart);
            %rshift = [xshift yshift zshift];
            rmask = double(tom_rotate(mask,[phi_old,psi_old,the_old]));
            %rshift = tom_pointrotate(rshift,phi_old,psi_old,the_old);
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
                psi = psi_old;
                the = the_old;
                rpart=double(tom_rotate(ref,[phi,psi,the]));
                rpart=tom_ifourier(ifftshift(tom_spheremask(wedge.*fftshift(tom_fourier(rpart)))));%changed back F
                rmask=double(tom_rotate(mask,[phi,psi,the]));
                rpart = rpart.*rmask;%mask with smoothened edges
                rpart = rpart - rmask.*(sum(sum(sum(rmask.*rpart)))/npixels); %subtract mean in mask
                fref=fftshift(tom_fourier(rpart));
                %apply bandpass
                fref=ifftshift(tom_spheremask(fref,lowpass,3) - tom_spheremask(fref,hipass,2));
                fref(1,1,1)=0;
                %calculate rms - IMPORTANT! - missing wedge!!!
                % changed back FF
                % fref = fref/sqrt((sum(sum(sum(fftshift(fref.*conj(fref)).*wedge)))));
                fref = fref/sqrt((sum(sum(sum(fref.*conj(fref))))));
                ccf = npixels*tom_spheremask(real(fftshift(tom_ifourier(fpart.*conj(fref)))),size(average,1)/5,size(average,1)/16);
                [pos ccctmp] = peak_det_2(real(ccf));
                if ccctmp > ccc
                    ccc = ccctmp;
                    phi_opt=phi;
                    psi_opt=psi;
                    the_opt=the;
                    tshift = pos-cent;
                    if size(ccf)>100
                        tom_dspcub(tom_bin(ccf));drawnow;
                    else
                        tom_dspcub(ccf);drawnow;
                    end;
                end;
            end;
            motl(17,indpart)=phi_opt;
            motl(18,indpart)=psi_opt;
            motl(19,indpart)=the_opt;
            motl(11,indpart) = tshift(1)*(2^ibin);
            motl(12,indpart) = tshift(2)*(2^ibin);
            motl(13,indpart) = tshift(3)*(2^ibin);
            rshift = tom_pointrotate(tshift,-psi_opt,-phi_opt,-the_opt);
            motl(14,indpart) = rshift(1)*(2^ibin);
            motl(15,indpart) = rshift(2)*(2^ibin);
            motl(16,indpart) = rshift(3)*(2^ibin);
            motl(1,indpart)=ccc;
            % take care: particle4av is NOT pre-shifted
            if (ccc > threshold*meanv) %kick off bad particles
                average = average + double(tom_rotate(tom_shift(particle4av,-tshift),[-psi_opt,-phi_opt,-the_opt]));
                if size(average,1)>100
                    tom_dspcub((tom_bin(average)));drawnow;
                else
                    tom_dspcub((average));drawnow;
                end;
                %weighting - avoid interpolation artefacts
                tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi_opt,-phi_opt,-the_opt])),0.5,1,'z'),0,0.5);
                wei = wei + tmpwei;
                motl(20,indpart)=1;%good particles -> class one
            else
                motl(20,indpart)=2;%bad CCF: kick into class 2
            end;
            disp(['Particle no ' num2str(indpart) ' , Iteration no ' num2str(ind)]);
            if ( rem(indpart,10) == 0 )
                name = [motlfilename '_tmp_' num2str(ind+1) '.em'];
                tom_emwrite(name,motl);
            end;
        end; %endif 
    end;% end particle loop
    name = [motlfilename '_' num2str(ind+1) '.em'];
    tom_emwrite(name,motl);
    % do weighting
    lowp = floor(size(average,1)/2)-3;%lowpass
    wei = 1./wei;rind = find(wei > 100000);wei(rind) = 0;% take care for inf
    average = real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(average)).*wei,lowp))));
    if ibin > 0
        average = tom_zoom(average,ibin);
    end;
    name = [refilename '_' num2str(ind+1) '.em'];
    tom_emwrite(name,average);
    disp(['wrote reference ' name]);
end; % end iteration loop
