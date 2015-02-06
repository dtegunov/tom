
function av3_scan_fast(refilename, motlfilename, particlefilename, startindx, iterations,angincr, ...
    angiter,mask,hipass,lowpass,nfold,threshold,wedgelist,iclass, ibin, dispflag, limit)
% av3_scan_fast aligns 3D subtomograms to reference
%
%   av3_scan_fast (refilename, motlfilename, particlefilename, startindx, iterations,angincr,...
%       angiter,mask,hipass,lowpass,nfold,threshold,wedgelist,iclass, ibin, dispflag, limit)
%   Filenames are expected as:
%       'filename'_#no.em
%
%   In this routine, the reference is rotated into the proposed orientation
%   (given in the motl). A X-correlation with the particle under scrutiny is
%   performed for this orientation and the ANGITER*ANGINCR neighbouring
%   orienations. 
%   A clever sampling on the unit sphere is done for psi and theta: rings around 
%   (theta_old, psi_old) are drawn. 
%   The routine takes into account the missing wedge - the extreme tilt
%   angles have to be specified in WEDGELIST. The specified MASK is used 
%   for masking the individual particles and the reference prior to 
%   X-correlation. Preferrably the mask should have smoothened edges. 
%   Only particles of the class ICLASS (column 20 of MOTL) are aligned.
%   However, take care; classes '1' and '2' are reserved in the following
%   way: Particles of these classes are ALWAYS aligned. Particles of class
%   '1' are expected to contribute to the average, particles of class '2'
%   had X-correlation coefficients below the specified THRESHOLD in the
%   previous iteration (if THRESHOLD >0). A clarifying example: Given, you 
%   start with a MOTL containing particles of classes 5 and 6 and a
%   THRESHOLD=0. You chose to align solely the particles of ICLASS=5. After 
%   finishing one cycle you will end up with particles of the CLASSES 1 and
%   6. The particles of class 6 were not aligned, the particles of class 5
%   were aligned and all contributed to the resulting average. If you now
%   raise the THRESHOLD to say 0.9 you will end up with particles of
%   classes 1, 2, and 6 after the next cycle: The particles of classes 1
%   and 2 were aligned, but only the particles of class 1 contributed to
%   the average.
%   For determination of CCF-peaks spline interpolation is used leading to
%   sub-pixel accuracy.
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
%   mask                mask - make sure dims are all right!
%   hipass              hipass - for X-corr
%   lowpass             lowpass - for X-corr
%   nfold               symmetry of particle (rotational symmetry along z 
%                           (e.g. 3 for 3-fold rotational symmetry)
%   threshold           Threshold*mean(ccc) is cutoff for averaging - e.g.
%                           choose 0.5 (mean of ccc from PREVIOUS iteration!)
%   wedgelist           array containing the tilt range - 1st column: no of
%                           tomogram, 2nd: minimum tilt, 3rd: max tilt
%   iclass              class of particles - default:0
%   ibin                binning - default 0
%   dispflag            flag for display - if set to 'nodisp' then no
%                       display, e.g. for batch mode - otherwise display.
%                       'Dspcub' views of average and X-correlation
%                       function are shown. - default: 'disp'
%   limit               optional : to limit the histogram of single particles
%                       and the reference, give as a factor of how many times
%                       std to cut, cut off intensities will be set to mean (lim 'z')
%                       default (omitted) : no normmalisation, no cutting
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
%       11        : x-shift in subvolume - AFTER rotation of template
%       12        : y-shift in subvolume - AFTER rotation of template
%       13        : z-shift in subvolume - AFTER rotation of template
%     ( 14        : x-shift in subvolume - BEFORE rotation of template )
%     ( 15        : y-shift in subvolume - BEFORE rotation of template )
%     ( 16        : z-shift in subvolume - BEFORE rotation of template )
%       17        : Phi (in deg)
%       18        : Psi
%       19        : Theta 
%       20        : class no
%
%  SEE ALSO
%   AV3_PHIALIG, AV3_CREATEMOTL, molmatch.exe, omnimatch.exe
%
%   09/18/03 FF
% modified 03/21/05 FF - corrected bug in binning
%
% last change 14/11/05, should be almost 3x faster,
% adapted from Vince (don't waste time on masking),
% limit parameter introduced,
% bugfix in threshold calculation,
% tested KG & MB

% start new masking stuff, 14/11/05
mask_template = ones( size(mask,1),size(mask,2),size(mask,3));
mask_lowpass_3 = tom_spheremask( mask_template, lowpass,3);
mask_hipass_2 = tom_spheremask( mask_template,hipass,2);
mask_x = mask_lowpass_3-mask_hipass_2;
clear mask_lowpass_3;
clear mask_hipass_2;
mask_y = tom_spheremask(mask_template,size(mask,1)/5,size(mask,1)/16);
% do some funny masking inline :
mcenter=[floor(size(mask,1)/2)+1, floor(size(mask,2)/2)+1, floor(size(mask,3)/2)+1];
[x,y,z]=ndgrid(0:size(mask,1)-1,0:size(mask,2)-1,0:size(mask,3)-1);
radius = floor((min(min(size(mask,1),size(mask,2)),size(mask,3))-1)/2) ;
x = sqrt((x+1-mcenter(1)).^2+(y+1-mcenter(2)).^2+(z+1-mcenter(3)).^2);
ind = find(x>=radius);
clear x y z mcenter;
mask_z = ones(size(mask,1), size(mask,2), size(mask,3));
mask_z(ind) = 0;
% end new masking stuff, 14/11/05

error(nargchk(13,17,nargin))
if nargin < 14 
    iclass = 0;
end;
if nargin < 15
    ibin = 0;
end;
if nargin < 16
    dispflag = 'disp';
end;
if ibin > 0
    mask = tom_bin(mask,ibin);
    mask_x = tom_bin(mask_x,ibin);mask_y = tom_bin(mask_y,ibin);mask_z = tom_bin(mask_z,ibin); % added 14/11/05, MB
    % account for binning in rotation center
    rcent = floor(size(mask)./2)-0.25;% so far only 'experimentally' checked for even dim
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
if nargin < 15, 
    disp(['Starting cycle of ' num2str(iterations) ' iterations ...']);
elseif nargin <= 16,
    disp(['Starting cycle of ' num2str(iterations) ' iterations with binning of ' num2str(ibin)]);
elseif nargin > 16,
    disp(['Starting cycle of ' num2str(iterations) ' iterations with binning of ' num2str(ibin) ' and limit of ' num2str(limit) 'x rms']);
end;
    %wedge=tom_wedge(ref.Value,semiangle);
for ind = startindx:startindx+iterations-1
    name = [refilename '_' num2str(ind) '.em'];
    ref = tom_emread(name);
    disp(['read in file ' name]);
    ref = ref.Value;average=ref*0;
    wei = zeros(size(ref,1),size(ref,2),size(ref,3));%weighting function
    if ibin > 0
        ref = tom_bin(ref,ibin);
    end;
    ref = tom_symref(ref,nfold);
    [mref xx1 xx2 mstd] = tom_dev(ref,'noinfo');
    ref = (ref - mref)./mstd;
    %ref = tom_limit(ref,-limit,limit,'z'); % throw away the gold, TEMPORARILY TAKEN OUT MB 31/11/05
    [mref xx1 xx2 mstd] = tom_dev(ref,'noinfo');
    ref = (ref - mref)./mstd;
    %mask = tom_spheremask(ones(size(average,1),size(average,2),size(average,3)),roi,rsmooth);
    name = [motlfilename '_' num2str(ind) '.em'];
    motl = tom_emread(name);
    motl = motl.Value;
    if ibin > 0
        motl(11:16) = motl(11:16)/(2^ibin);%take binning into account
    end;
    indx = find ((motl(20,:) ==1 ) | (motl(20,:) == 2) | (motl(20,:) == iclass) ); meanv = mean(motl(1,indx));
    disp(['mean ccc of ' num2str(size(indx,2)) ' particles in classes 1, 2 and ' num2str(iclass) ' is: ' num2str(meanv)]);
    indx = find (motl(1,indx) > threshold*meanv);
    disp(['at given threshold of ' num2str(threshold) ', ' num2str(size(indx,2)) ' particles have a ccc larger than threshold*meanv: ' num2str(threshold*meanv)]);
    itomo_old = 0;
    for indpart = 1:size(motl,2)
        if ((motl(20,indpart) == 1) | (motl(20,indpart) == 2) | (motl(20,indpart) == iclass))
            itomo = motl(5,indpart);
            if itomo_old ~= itomo %wedge stuff - exact weighting according to list
                xx = find(wedgelist(1,:) == itomo);
                minangle= wedgelist(2,xx);
                maxangle= wedgelist(3,xx);
                wedge = av3_wedge(ref,minangle,maxangle);
                if (ibin > 0)
                    bigwedge = av3_wedge(average,minangle,maxangle);
                end;
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
            particle4av = particle;
            if ibin > 0
                particle = tom_bin(particle,ibin);
            end;
            if nargin > 16, % throw away the gold, modified 14/11/05, MB
                [pmean pmax pmin pstd] = tom_dev(particle,'noinfo');
                particle = (particle - pmean)./pstd;
                particle = tom_limit(particle,-limit,limit,'z'); 
            end; % end modify by MB
            % do not shift particle but mask !
            rshift(1) = motl(11,indpart);
            rshift(2) = motl(12,indpart);
            rshift(3) = motl(13,indpart);
            %rshift = [xshift yshift zshift];
            if (ibin == 0)
                rmask = double(tom_rotate(mask,[phi_old,psi_old,the_old]));
            else
                rmask = double(tom_rotate(mask,[phi_old,psi_old,the_old],'linear',rcent));
            end;
            %rshift = tom_pointrotate(rshift,[phi_old,psi_old,the_old]);
            shiftmask = tom_shift(rmask,rshift);
            particle=shiftmask.*particle;
            particle= particle - shiftmask.*(sum(sum(sum(particle)))/npixels);%subtract mean in sphere
            fpart=fftshift(tom_fourier(particle));
            %apply bandpass
            %
            %fprintf(1,'run_av3 size(fpart) %d % d  %d \n',size(fpart));
            fpart= ifftshift(fpart.*mask_x);
            %
            %fpart= ifftshift(tom_spheremask(fpart,lowpass,3) - tom_spheremask(fpart,hipass,2));
            %normalize
            fpart(1,1,1)=0;
            fpart = (size(fpart,1)*size(fpart,2)*size(fpart,3))*fpart/sqrt((sum(sum(sum(fpart.*conj(fpart))))));
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
                        if (ibin == 0)
                            rpart=double(tom_rotate(ref,[phi,psi,the]));
                        else
                            rpart=double(tom_rotate(ref,[phi,psi,the],'linear',rcent));
                        end;
                        %fprintf(1,'run_av3 size(rpart) %d % d  %d \n',size(rpart));
                        %changed back F
                        %rpart=tom_ifourier(ifftshift(tom_spheremask(wedge.*fftshift(tom_fourier(rpart)))));%changed back F

                        if (ibin == 0)
                            rmask=double(tom_rotate(mask,[phi,psi,the]));
                        else
                            rmask=double(tom_rotate(mask,[phi,psi,the],'linear',rcent));
                        end;
                        rpart = rpart.*rmask;%mask with smoothened edges
                        rpart = rpart - rmask.*(sum(sum(sum(rmask.*rpart)))/npixels); %subtract mean in mask
                        fref=fftshift(tom_fourier(rpart));
                        %apply bandpass
                        %fref=ifftshift(tom_spheremask(fref,lowpass,3) - tom_spheremask(fref,hipass,2));
                        % 
                        %fprintf(1,'run_av3 size(fref) %d % d  %d \n',size(fref));
                        fref=ifftshift(wedge.*fref.*mask_z.*mask_x);
                       
                        %
                        fref(1,1,1)=0;
                        %calculate rms - IMPORTANT! - missing wedge!!!
                        % changed back FF
                        % fref = fref/sqrt((sum(sum(sum(fftshift(fref.*conj(fref)).*wedge)))));
                        fref = (size(fref,1)*size(fref,2)*size(fref,3))*fref/sqrt((sum(sum(sum(fref.*conj(fref))))));% to be changed?
                        %ccf = tom_spheremask(real(fftshift(tom_ifourier(fpart.*conj(fref)))),size(ref,1)/5,size(ref,1)/16);
                        %fprintf(1,'run_av3 size(fpart) %d % d  %d \n',size(fpart));
                        ccf = real(fftshift(tom_ifourier(fpart.*conj(fref)))).*mask_y;

                        ccf = ccf/(size(ccf,1).^3);% added for normalization - FF
                        [pos ccctmp] = peak_det_2(real(ccf));
                        if ccctmp > ccc
                            ccc = ccctmp;
                            phi_opt=phi;
                            psi_opt=psi;
                            the_opt=the;
                            tshift = pos-cent;
                            if (strcmp(dispflag,'nodisp')~= 1 )
                                if (size(ccf)>100)
                                    binf = int8(log2(size(ccf)/64));
                                    tom_dspcub(tom_bin(ccf),binf(1));drawnow;
                                else
                                    tom_dspcub(ccf);drawnow;
                                end;
                            end;
                        end;
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
                if (ibin == 0)
                    average = average + double(tom_rotate(tom_shift(particle4av,-tshift),[-psi_opt,-phi_opt,-the_opt]));
                else
                    average = average + double(tom_rotate(tom_shift(particle4av,-tshift*(2^ibin)),[-psi_opt,-phi_opt,-the_opt]));
                end;
                if (strcmp(dispflag,'nodisp')~= 1 )
                    if size(average,1)>100
                        tom_dspcub((tom_bin(average)));drawnow;
                    else
                        tom_dspcub((average));drawnow;
                    end;
                end;
                %weighting - avoid interpolation artefacts
                if (ibin ==0)
                    tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi_opt,-phi_opt,-the_opt])),0.5,1,'z'),0,0.5);
                else
                    tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(bigwedge,[-psi_opt,-phi_opt,-the_opt])),0.5,1,'z'),0,0.5);
                end;
                wei = wei + tmpwei;
                motl(20,indpart)=1;%good particles -> class one
            else
                motl(20,indpart)=2;%bad CCF: kick into class 2
            end;
            disp(['Particle no ' num2str(ifile) ' , Iteration no ' num2str(ind)]);
            if ( rem(indpart,10) == 0 )
                name = [motlfilename '_tmp_' num2str(ind+1) '.em'];
                tom_emwrite(name,motl);
            end;
        end; %endif 
    end;% end particle loop
    name = [motlfilename '_' num2str(ind+1) '.em'];
    tom_emwrite(name,motl);
    ixx = find (motl(20,:)==1);disp([ num2str(size(ixx,2)) ' particles are in class 1 after iteration No. ' num2str(ind)]);
    % do weighting
    lowp = floor(size(average,1)/2)-3;%lowpass
    wei = 1./wei;rind = find(wei > 100000);wei(rind) = 0;% take care for inf
    average = real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(average)).*wei,lowp))));
%     if ibin > 0% fixed FF 03/21/05
%         average = tom_zoom(average,ibin);
%     end;
    name = [refilename '_' num2str(ind+1) '.em'];
    tom_emwrite(name,average);
    disp(['wrote reference ' name]);
end; % end iteration loop

