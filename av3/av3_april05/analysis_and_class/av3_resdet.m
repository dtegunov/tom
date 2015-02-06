function [average,average2,frc] = av3_resdet(motl, particlefilename,wedgelist,num,nsym,threshold)
% AV3_RESDET determines resolution of average
%
%   [average average2 frc] = av3_resdet(motl, particlefilename,wedgelist,num,nsym,threshold)
%
%   Procedure for resolution determination of particle reconstructed from
%   tomograms. Resolution determination has to be performed according to
%   the resulting Fourier-ring correlation FRC.
%
% PARAMETERS
%  INPUT
%   MOTL                Motive list
%   PARTICLEFILENAME    particlefilenames - individual tomograms are
%                       expected to be called <PARTICLEFILENAME>_no.em 
%   WEDGELIST           list of wedges
%   NUM                 Number of rings for FRC
%   NSYM                SYMMETRY
%   THRESHOLD           threshold for decrimination of particles according
%                       to CCC. 
%
%  OUTPUT
%   AVERAGE             average of 1st set of particles
%   AVERAGE2            average of 2nd set of particles
%   FRC                 Output stack of tom_frc
%
%   11/03/03 FF
%   last change 27/01/05 FF

indx = find (motl(1,:) > 0); meanv = mean(motl(1,indx));
indx = find (motl(1,:) > threshold*meanv);

% define average in advance
name = [particlefilename '_' num2str(motl(4,1)) '.em'];%last change
particle = tom_emread(name);particle = particle.Value;
wei = zeros(size(particle,1),size(particle,2),size(particle,3));
wei2 =wei;
average = wei;
average2 = wei;

itomo_old = 0;
icount = 0;
icount2 = 0;
tag =1;%tag for splitting particles into 2 parts
for indpart = 1:size(motl,2) 
    if (motl(1,indpart)>threshold*meanv ) & (motl(20,indpart)> 0)
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
        if indpart == 1
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
        if tag == 1
            average = average + double(tom_rotate(tom_shift(particle,-tshift),[-psi,-phi,-the]));
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi,-phi,-the])),0.5,1,'z'),0,0.5);
            wei = wei + tmpwei;
            icount = icount +1;
            disp(['Particle no ' num2str(ifile) ' added to average1'  ]);
            tag = 2;
        else
            average2 = average2 + double(tom_rotate(tom_shift(particle,-tshift),[-psi,-phi,-the]));
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi,-phi,-the])),0.5,1,'z'),0,0.5);
            wei2 = wei2 + tmpwei;
            icount2 = icount2 +1;
            disp(['Particle no ' num2str(ifile) ' added to average2'  ]);
            tag = 1;
        end;%if - even/odd
    end;%if - threshold
end;
lowp = floor(size(average,1)/2)-3;
wei = 1./wei;rind = find(wei > 100000);wei(rind) = 0;% take care for inf
average = real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(average)).*wei,lowp))));
wei2 = 1./wei2;rind = find(wei2 > 100000);wei2(rind) = 0;
average2 = real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(average2)).*wei2,lowp))));
disp(['Averaging finished - ' num2str(icount) ' particles averaged to average1 - '  ...
    num2str(icount2) ' particles averaged to average2 ' ]);
%average = tom_symref(average,nsym);
%average2 = tom_symref(average2,nsym);
frc = tom_compare(average, average2, num);
