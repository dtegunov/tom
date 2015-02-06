function average = av3_average_exact(motl, particlefilename,wedgelist,threshold,iclass)
%  AV3_AVERAGE_EXACT averages subtomograms
%   average = av3_average_exact(motl, particlefilename,wedgelist,threshold,iclass)
%
%   AV3_AVERAGE_EXACT is designed for averaging subtomogram 
%   (PARTICLEFILENAME = 'filename'_#no.em) using the parameters of the
%   MOTL. If THRESHOLD is specified only particle with a CCC >
%   thresold*mean(ccc) are included into the average. If ICLASS is
%   specified only particles of this class will be included into the
%   average. 
%
%   motl                motivelist (see AV3_TRANS-ROT-ALIG for format)
%   particlefilename    filename of 3D-particles to be averaged
%                           'particlefilename'_#no.em
%   wedgelist           array containing the tilt range - 1st column: no of
%                           tomogram, 2nd: minimum tilt, 3rd: max tilt
%   threshold           Threshold*mean(ccc) is cutoff for averaging
%   iclass              class of particles - default:0
%   
%
%   SEE ALSO
%   AV3_TRANS-ROT-ALIG, AV3_PHIALIG
%
% last change 03/31/05 FF - doc updated

if nargin<5
    iclass = 0;
end;
if nargin < 4
    threshold = 0;
end;
indx = find (motl(1,:) >= 0); meanv = mean(motl(1,indx));% find ultimate solution!!!!
% indx = find (motl(1,:) > threshold*meanv);
% motl = motl(:,indx);

% define average in advance
%name = [particlefilename '_' num2str(1) '.em'];
%particle = tom_emread(name);particle = particle.Value;
%wei = zeros(size(particle,1),size(particle,2),size(particle,3));
%average = wei;

itomo_old = 0;
icount = 0;
for indpart = 1:size(motl,2) 
    if (motl(1,indpart)>=threshold*meanv & motl(20,indpart) == iclass)
        icount = icount +1;
        itomo = motl(5,indpart);
        xshift = motl(11,indpart);
        yshift = motl(12,indpart);
        zshift = motl(13,indpart);
        tshift = [xshift yshift zshift];
        phi=motl(17,indpart);
        psi=motl(18,indpart);
        the=motl(19,indpart);
        ifile = motl(4,indpart);
        name = [particlefilename '_' num2str(ifile) '.em'];
        particle = tom_emread(name);particle = particle.Value;
        particle = tom_limit(particle,-3,4,'z'); % throw away the gold
        if icount == 1
            wei = zeros(size(particle,1),size(particle,2),size(particle,3));
            average = wei;
        end;
        if itomo_old ~= itomo %wedge stuff - exact weighting according to list
            xx = find(wedgelist(1,:) == itomo);
            minangle= wedgelist(2,xx);
            maxangle= wedgelist(3,xx);
            wedge = av3_wedge(particle,minangle,maxangle);
            itomo_old = itomo;
        end;
        average = average + double(tom_rotate(tom_shift(particle,-tshift),[-psi,-phi,-the]));
        tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi,-phi,-the])),0.5,1,'z'),0,0.5);
        wei = wei + tmpwei;
        disp(['Particle no ' num2str(ifile) ' added to average'  ]);
    end;%if - threshold
end;
lowp = floor(size(average,1)/2)-3;
wei = 1./wei;rind = find(wei > 100000);wei(rind) = 0;% take care for inf
average = real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(average)).*wei,lowp))));
disp(['Averaging finished - ' num2str(icount) ' particles averaged ... '  ]);
