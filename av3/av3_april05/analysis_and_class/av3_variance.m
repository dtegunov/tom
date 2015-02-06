function [vari, vari_unw] = av3_variance(motl, particlefilename,wedgelist,threshold,iclass)
% AV3_VARIANCE computes variance of average
%
%   [vari, vari_unw] = av3_variance(motl, particlefilename,wedgelist,threshold,iclass)
%   
%   AV3_VARIANCE calculates spatial variation of density, i.e.
%   var = sum_i((rho(x)-av(x))^2)  - i=index of particle
%   routine takes missing wedge into account.
%
% PARAMETERS
%  INPUT
%   motl                moivelist
%   particlefilename    filename of 3D-particles to be aligned and averaged
%   wedgelist           array containing the tilt range - 1st column: no of
%   threshold           Threshold*mean(ccc) is cutoff for averaging - e.g.
%                           choose 0.5 (mean of ccc from PREVIOUS iteration!)
%   iclass              class of particles - default:0
%
%  OUTPUT
%   vari                variance - weighted according to distribution of
%                           wedges
%   vari_unw            variance - not-weighted
%   
%
% modified 18/02/04
% last change 04/01/05 FF - updated docu

if nargin<5
    iclass = 0;
end;
indx = find (motl(1,:) > 0); meanv = mean(motl(1,indx));
indx = find (motl(1,:) > threshold*meanv);

itomo_old = 0;
icount = 0;
for indpart = 1:size(motl,2) 
    if (motl(1,indpart)>threshold*meanv & motl(20,indpart) == iclass)
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
            vari = wei;
        end;
        if itomo_old ~= itomo %wedge stuff - exact weighting according to list
            xx = find(wedgelist(1,:) == itomo);
            minangle= wedgelist(2,xx);
            maxangle= wedgelist(3,xx);
            wedge = av3_wedge(particle,minangle,maxangle);
            itomo_old = itomo;
        end;
        particle = real(tom_ifourier(ifftshift(fftshift(tom_fourier(particle)).*wedge)));
        rpart = double(tom_rotate(tom_shift(particle,-tshift),[-psi,-phi,-the]));
        average = average + rpart;
        vari = vari + rpart.^2;
        tom_dspcub(vari);
        tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi,-phi,-the])),0.5,1,'z'),0,0.5);
        wei = wei + tmpwei;
        disp(['Particle no ' num2str(ifile) ' added to average'  ]);
    end;%if - threshold
end;
vari_unw = vari-average.^2/icount;
lowp = floor(size(average,1)/2)-3;
wei = 1./wei;rind = find(wei > 100000);wei(rind) = 0;% take care for inf
vari = real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(vari-average.^2/icount)).*wei,lowp))));
disp(['Averaging finished - ' num2str(icount) ' particles averaged ... '  ]);
disp(['start to calculate variance ']);
