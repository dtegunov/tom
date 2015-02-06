function av3_spalig_phi2(iref, refilename, motlfilename, particlefilename, startindx, mask, ...
    hipass, lowpass, nfold, threshold, wedgelist,iclass, ibin)  
% AV3_SPALIG_PHI aligns 3D subtomograms to reference
%
%   av3_spalig_phi(iref, refilename, motlfilename, particlefilename, startindx, mask, ...
%       hipass, lowpass, nfold, threshold, wedgelist,iclass, ibin)  
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
%last change 11/05/04

error(nargchk(11,13,nargin))
if nargin < 12
    iclass = 0;
end;
if nargin < 13
    ibin = 0;
end;
if ibin > 0
    mask = tom_bin(mask,ibin);
end;
npixels = sum(sum(sum(mask)));
cent= [floor(size(mask,1)/2)+1 floor(size(mask,2)/2)+1 floor(size(mask,3)/2)+1];
centphi = 2*floor(size(mask,1)/2)+1;
scf = size(mask,1)*size(mask,2)*size(mask,3);
ind = startindx;
%wedge=tom_wedge(ref.Value,semiangle);
name = [motlfilename '_' num2str(ind) '.em'];
motl = tom_emread(name);
motl = motl.Value;
if ibin > 0
    motl(11:16) = motl(11:16)/(2^ibin);%take binning into account
end;
indref = find(motl(4,:) == iref);% determine 
name = [particlefilename '_' num2str(iref) '.em'];
ref = tom_emread(name);
disp(['read in file ' name]);
ref = ref.Value;
if ibin > 0
    ref = tom_bin(ref,ibin);
end;
tshift(1) = motl(11,indref)/(2^ibin);tshift(2) = motl(12,indref)/(2^ibin);tshift(3) = motl(13,indref)/(2^ibin);
phi_opt=motl(17,indref);psi_opt=motl(18,indref);the_opt=motl(19,indref);
ref = double(tom_rotate(tom_shift(ref,-tshift),[-psi_opt,-phi_opt,-the_opt]));
ref = tom_symref(ref,nfold);
[mref xx1 xx2 mstd] = tom_dev(ref,'noinfo');
ref = (ref - mref)./mstd;
%ref = tom_limit(ref,-3,4,'z'); % throw away the gold
[mref xx1 xx2 mstd] = tom_dev(ref,'noinfo');
ref = (ref - mref)./mstd;
average=ref;
wei = zeros(size(average,1),size(average,2),size(average,3));%weighting function
%mask =
%tom_spheremask(ones(size(average,1),size(average,2),size(average,3)),roi,rsmooth);
indx = find (motl(1,:) > 0); meanv = mean(motl(1,indx));
indx = find (motl(1,:) > threshold*meanv);
itomo_old = 0;
for indpart = 1:size(motl,2)
    if ( ((motl(20,indpart) == 1) | (motl(20,indpart) == 2) | (motl(20,indpart) == iclass)) && (indpart~=indref) )
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
        tshift(1) = motl(11,indpart);tshift(2) = motl(12,indpart);tshift(3) = motl(13,indpart);
        ifile = motl(4,indpart);
        name = [particlefilename '_' num2str(ifile) '.em'];
        particle = tom_emread(name);particle = particle.Value;
        if ibin > 0
            particle = tom_bin(particle,ibin);
        end;
        [mref xx1 xx2 mstd] = tom_dev(particle,'noinfo');
        particle = (particle - mref)./mstd;
        particle = tom_limit(particle,-3,4,'z'); % throw away the gold
        particle4av = particle;
        %apply bandpass
        fpart=fftshift(tom_fourier(particle));
        fpart= ifftshift(tom_spheremask(fpart,lowpass,3) - tom_spheremask(fpart,hipass,2));
        particle = real(tom_ifourier(fpart));
        particle = double(tom_rotate(tom_shift(particle,-tshift),[-psi_old,-phi_old,-the_old]));
        particle = mask.*particle;
        particle= particle - mask.*(sum(sum(sum(particle)))/npixels);%subtract mean in sphere
        pol = tom_cart2cyl(particle);[rpol phii zi] = ndgrid(1:size(pol,1),1:size(pol,2),1:size(pol,3));
        rpol = rpol.*pol;clear phii zi;
        fpol = fft(rpol,[],2);
        % prepare ref
        %rpart=tom_ifourier(ifftshift(tom_spheremask(wedge.*fftshift(tom_fourier(rpart)))));
        fav = fftshift(tom_fourier(average));
        fav = ifftshift(tom_spheremask(fav,lowpass,3) - tom_spheremask(fav,hipass,2));
        ref = real(tom_ifourier(fav));
        ref = mask.*ref;
        ref = ref - mask.*(sum(sum(sum(ref)))/npixels);%subtract mean in sphere
        refpol = tom_cart2cyl(ref);[rrefpol phii zi] = ndgrid(1:size(pol,1),1:size(pol,2),1:size(pol,3));
        rrefpol = rrefpol.*refpol; clear phii zi;
        fref = fft(rrefpol,[],2);
        %compute polar CCF
        ccf = fftshift(real( ifft(fpol.*conj(fref),[],2) ));
        phiccf = sum(sum(ccf,1),3);
        [phi_opt ccc] = peak_det_2(phiccf);
        deltaphi = (phi_opt(2)-centphi)*360/size(phiccf,2);
        disp(num2str(deltaphi) )
        phi_opt = (phi_old + deltaphi);
        psi_opt = psi_old;
        the_opt = the_old;
        motl(17,indpart) = phi_opt;
        %motl(18,indpart) = psi_opt;
        %motl(19,indpart) = the_opt;
        tshift(1) = motl(11,indpart);
        tshift(2) = motl(12,indpart);
        tshift(3) = motl(13,indpart);
        rshift = tom_pointrotate(tshift,-psi_opt,-phi_opt,-the_opt);
        motl(14,indpart) = rshift(1);
        motl(15,indpart) = rshift(2);
        motl(16,indpart) = rshift(3);
        %motl(1,indpart)=ccc;
        % take care: particle4av is NOT pre-shifted
        if (ccc > threshold*meanv) %kick off bad particles
            average = average + double(tom_rotate(tom_shift(particle4av,-tshift),[-psi_opt,-phi_opt,-the_opt]));
            if size(average,1)>100
                tom_dspcub((tom_bin(average)),1);drawnow;
            else
                tom_dspcub((average),1);drawnow;
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

