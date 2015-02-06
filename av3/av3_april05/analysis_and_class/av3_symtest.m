function [ccfsum, phiccf] = av3_symtest(motlfilename, particlefilename, startindx, mask, ...
    hipass, lowpass, threshold, iclass, ibin)  
% AV3_SYMTEST perfroms rotational symmetry test
%
%   [ccfsum, phiccf] = av3_symtest(motlfilename, particlefilename, startindx, mask, ...
%               hipass, lowpass, threshold, iclass, ibin)  
%
%
% PARAMETERS
%  INPUT
%   iref                particle number of 1st reference - no of particle
%                           file
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
%  OUTPUT
%   ccfsum              rotational auto-correlation function summed over all
%                           particles
%   phiccf              stack containg the individual ARCFs
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
%last change 10/25/04 (bug ibin resolved)

error(nargchk(7,9,nargin))
if nargin < 8
    iclass = 0;
end;
if nargin < 9
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
indx = find (motl(1,:) > 0); meanv = mean(motl(1,indx));
indx = find (motl(1,:) > threshold*meanv);
ir =0;
for ii = 1:size(motl,2)
    indpart = ii;
    if ( ((motl(20,indpart) == 1) | (motl(20,indpart) == 2) | ...
            (motl(20,indpart) == iclass)) )
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
        %particle = tom_limit(particle,-3,4,'z'); % throw away the gold
        %apply bandpass
        fpart=fftshift(tom_fourier(particle));
        fpart= ifftshift(tom_spheremask(fpart,lowpass,3) - tom_spheremask(fpart,hipass,2));
        particle = real(tom_ifourier(fpart));
        particle = double(tom_rotate(tom_shift(particle,-tshift),[-psi_old,-phi_old,-the_old]));
        particle = mask.*particle;
        particle= particle - mask.*(sum(sum(sum(particle)))/npixels);%subtract mean in sphere
        pol = tom_cart2cyl(particle);
        if ir == 0
            [r phii zi] = ndgrid(1:size(pol,1),1:size(pol,2),1:size(pol,3));clear phii zi;
            phiccf = zeros(size(motl,2),size(pol,2));
            ccfsum = zeros(1,size(pol,2));
            phivals = ((1:size(pol,2))-1)*360/(size(pol,2));
            ir =1;
        end;
        rpol = pol;
        fpol = fft(rpol,[],2);
        ccf = (real( ifft(fpol.*conj(fpol),[],2) )).*r;
        phiccf(indpart,:) = sum(sum(ccf,1),3);%figure(2);plot(phivals,phiccf(indpart,:));title('polar CCF');
        ccfsum = ccfsum + phiccf(indpart,:);%figure(1);plot(phivals,ccfsum)
        disp(['Particle no ' num2str(indpart) ]);
    end; %endif 
end;% end particle loop
