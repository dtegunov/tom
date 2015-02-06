function ccc = av3_constr_ccc(templ, motl, particlefilename, wedgelist,mask,hipass,lowpass,iclass,ibin)
% AV3_CONSTR_CCC computes correlation of particles and template
%
%   ccc = av3_constr_ccc(templ, motl, particlefilename, wedgelist,  ... 
%           mask, hipass, lowpass, iclass, ibin)  
%
%   The cross correlation coefficient is calculated in real space, thus the
%   translation is NOT determined. The CCC that is calculated here is
%   constrained in the following sense: For computing the correlation of
%   any pair of particles, the data is constrained to the part of Fourier
%   space that is sampled by both data. This can be achieved by
%   multiplying a mask in Fourier space. This normalization is only done
%   when FLAG is set to 'norm'. If FLAG is 'unchanged', then the
%   CCC is computed without normalizing the particles to the standard
%   deviation. The default is 'norm' where normalization is performed.
%
% PARAMETERS
%  INPUT
%   templ               template - all particles are correlated
%                       (=projected) with this volume
%   particlefilename    filename of 3D-particles to be aligned and averaged
%                           'particlefilename'_#no.em
%   wedgelist           array containing the tilt range - 1st: no of
%                           tomogram, 2nd: minimum tilt, 3rd: max tilt
%   mask                mask - make sure dims are all right!
%   hipass              hipass - for X-corr
%   lowpass             lowpass - for X-corr
%   iclass              class of particles - fixed to 1 or 2
%   ibin                binning - default 0
%
%  OUTPUT
%   ccc                 constrained cross correlation matrix
%   covar               constrained covariance matrix
%
% EXAMPLE
%
% SEE ALSO
%   TOM_CORR, TOM_ORCD
%
%    Copyright (c) 2004
%    TOM toolbox for Electron Tomography
%    Max-Planck-Institute for Biochemistry
%    Dept. Molecular Structural Biology
%    82152 Martinsried, Germany
%    http://www.biochem.mpg.de/tom
%
%   last change
%   22/06/04 FF

if nargin < 9
    ibin = 0;
end;
if ibin>0
    mask = tom_bin(mask,ibin);
    templ = tom_bin(templ,ibin);
end;
if nargin < 8 
    iclass = 0;
end;
icount = 0;
itomo_old = 0;
npixels = sum(sum(sum(mask)));
for indpart1 = 1:size(motl,2) 
    if ( (motl(20,indpart1) == iclass) | (motl(20,indpart1) == 1) | (motl(20,indpart1) == 2) )
        icount = icount +1;
        itomo = motl(5,indpart1);
        xshift = motl(11,indpart1);
        yshift = motl(12,indpart1);
        zshift = motl(13,indpart1);
        tshift = [xshift yshift zshift]/(2^ibin);
        phi = motl(17,indpart1);
        psi = motl(18,indpart1);
        the = motl(19,indpart1);
        ifile1 = motl(4,indpart1);
        name = [particlefilename '_' num2str(ifile1) '.em'];
        particle = tom_emread(name);
        if ibin > 0
            particle = tom_bin(particle.Value,ibin);
        else
            particle = particle.Value;
        end;
        if itomo_old ~= itomo %wedge stuff - exact weighting according to list
            xx = find(wedgelist(1,:) == itomo);
            minangle= wedgelist(2,xx);
            maxangle= wedgelist(3,xx);
            wedge = av3_wedge(particle,minangle,maxangle);
            itomo_old = itomo;
        end;
        rpart1 = double(tom_rotate(tom_shift(particle,-tshift),[-psi,-phi,-the]));
        frpart1 = fftshift(tom_fourier(rpart1));
        frpart1 = tom_spheremask(frpart1,lowpass,3) - tom_spheremask(frpart1,hipass,2);
        tmpwei1 = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi,-phi,-the])),0.5,1,'z'),0,0.5);
        %wpix = sum(sum(sum(tmpwei1)));
        cpart1 = frpart1.*tmpwei1;
        cpart1 = real(tom_ifourier(ifftshift(cpart1)));
        cpart2 = fftshift(tom_fourier(templ)).*tmpwei1;
        cpart2 = tom_spheremask(cpart2,lowpass,3) - tom_spheremask(cpart2,hipass,2);
        cpart2 = real(tom_ifourier(ifftshift(cpart2)));
        cpart1 = cpart1.*mask;mn1 = (sum(sum(sum(cpart1))))/npixels;
        cpart1 = cpart1 - mask.*mn1;
        stv1 = sqrt(sum(sum(sum(cpart1.*cpart1))));
        cpart1 = cpart1/stv1;
        cpart2 = cpart2.*mask;mn2 = (sum(sum(sum(cpart2))))/npixels;
        cpart2 = cpart2 - mask.*mn2;
        stv2 = sqrt(sum(sum(sum(cpart2.*cpart2))));
        cpart2 = cpart2/stv2;
        ccc(indpart1) = sum(sum(sum(cpart1.*cpart2)));
    end;
    %ccc(indpart1,indpart1) = 1;
    disp(['Correlation computed for particle no ' num2str(ifile1) ' ...']);
end;
%ccc = ccc+eye(size(motl,2));
