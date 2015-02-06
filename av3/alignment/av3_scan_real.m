function av3_scan_real(refilename, motlfilename, particlefilename, startindx, iterations, nshifts, ...
    mask, hipass,lowpass, nfold, threshold, wedgelist)
%
%   av3_scan_real(refilename, motlfilename, particlefilename, startindx, iterations, nshifts, ...
%       mask, hipass,lowpass, nfold, threshold, wedgelist)
%
% 10/16/03 FF
% last change 11/05/04 FF

threshold=0;
npixels = sum(sum(sum(mask)));
imask = find(mask>0);
ind = startindx;
name = [refilename '_' num2str(ind) '.em'];
ref = tom_emread(name);
for ind = startindx:startindx+iterations-1
    % prepare reference
    name = [refilename '_' num2str(ind) '.em'];
    ref = tom_emread(name);
    disp(['read in file ' name]);
    [mref dummy1 dummy2 mvar]  = tom_dev(ref.Value,'noinfo');
    ref = (ref.Value - mref)/mvar;
    ref = tom_limit(ref,-3,4,'z'); %throw away the gold
    ref = tom_symref(ref,nfold);

    %initialize average and weighting Fkt
    %ref = tom_bandpass(ref,hipass,lowpass,2);
    average=zeros(size(ref,1),size(ref,2),size(ref,3));
    wei = zeros(size(ref,1),size(ref,2),size(ref,3));
    %MOTL
    name = [motlfilename '_' num2str(ind) '.em'];
    motl = tom_emread(name);
    motl = motl.Value;
    % mean of CCCs - for exclusion of particles
    indx = find (motl(1,:) > 0); meanv = mean(motl(1,indx));
    itomo_old = 0; % tag
    for indpart = 1:size(motl,2)
        if ((motl(20,indpart) == 1) | (motl(20,indpart) == 2))
            itomo = motl(5,indpart);
            if itomo_old ~= itomo %wedge stuff - exact weighting according to list
                xx = find(wedgelist(1,:) == itomo);
                minangle= wedgelist(2,xx);
                maxangle= wedgelist(3,xx);
                wedge = av3_wedge(ref,minangle,maxangle);
                itomo_old = itomo;
            end;%if -wedge
            phi_old=motl(17,indpart);
            psi_old=motl(18,indpart);
            the_old=motl(19,indpart);
            % -- multiply wedge on ref and normalize - only once for old angles! maybe alter --
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi_old,-phi_old,-the_old])),0.5,1,'z'),0,0.5);
            tref = tom_symref(real(tom_ifourier(ifftshift(fftshift(tom_fourier(ref)).*tmpwei))/(size(ref,1)*size(ref,2)*size(ref,3))),nfold);
            tref = mask.*tref;
            mref = sum(sum(sum(tref(imask))))/npixels;
            tref(imask) = tref(imask) - mref;
            mvar = sqrt( sum(sum(sum(tref(imask).*tref(imask)))) );
            tref(imask) = tref(imask)/mvar;
            % shift of particles  - after rotation
            xshift_old = round(motl(11,indpart));
            yshift_old = round(motl(12,indpart));
            zshift_old = round(motl(13,indpart));
            ifile = motl(4,indpart);
            ccc = -1;
            %read particles
            name = [particlefilename '_' num2str(ifile) '.em'];
            particle = tom_emread(name);particle = particle.Value;
            particle = tom_limit(particle,-3,4,'z'); % throw away the gold
            particle4av = particle;
            particle = tom_bandpass(particle,hipass,lowpass,2);
            for ishftx=-nshifts:nshifts
                for ishfty=-nshifts:nshifts
                    for ishftz=-nshifts:nshifts
                        radius = sqrt(ishftx^2+ishfty^2+ishftz^2);
                        if (radius >= nshifts)
                            shift = [xshift_old+ishftx yshift_old+ishfty zshift_old+ishftz];
                            phi = phi_old; psi=psi_old; the=the_old;
                            shift_part = double(tom_rotate(tom_move(particle,-shift),[-psi,-phi,-the]));
                            shift_part = mask.*tom_symref(shift_part,nfold);
                            %tom_dspcub(shift_part);
                            mpart = sum(sum(sum(shift_part(imask))))/npixels;
                            shift_part(imask) = shift_part(imask) - mpart;
                            pvar = sqrt(sum(sum(sum(shift_part(imask).*shift_part(imask)))));
                            shift_part(imask) = shift_part(imask)/pvar;
                            ccctmp = sum(sum(sum(shift_part(imask).*tref(imask))));
                            if ccctmp > ccc
                                ccc = ccctmp;
                                phi_opt=phi;
                                psi_opt=psi;
                                the_opt=the;
                                shift_opt = shift;
                            end;%if - max CCC
                        end; % if - circle around cent
                    end;%for - shifts in Z
                end;%for - shifts in Y
            end;%for - shifts in X
            %  -- enter into MOTL --
            motl(17,indpart)=phi_opt;
            motl(18,indpart)=psi_opt;
            motl(19,indpart)=the_opt;
            motl(11,indpart) = shift_opt(1);
            motl(12,indpart) = shift_opt(2);
            motl(13,indpart) = shift_opt(3);
            motl(1,indpart)=ccc;
            %  -- disciminate particles --
            if (ccc > threshold*meanv) %kick off bad particles
                average = average + double(tom_rotate(tom_move(particle4av,-shift_opt),[-psi_opt,-phi_opt,-the_opt]));
                tom_dspcub((average));
                %weighting - avoid interpolation artefacts
                tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi_opt,-phi_opt,-the_opt])),0.5,1,'z'),0,0.5);
                wei = wei + tmpwei;
                motl(20,indpart)=1;%good particles -> class one
            else
                motl(20,indpart)=2;%bad CCF: kick into class 2
            end;
            disp(['Particle no ' num2str(indpart) ' , Iteration no ' num2str(ind)]);
        end;%if - correct classes
    end;%for - particles
    name = [motlfilename '_' num2str(ind+1) '.em'];
    tom_emwrite(name,motl);
    % do weighting
    lowp = floor(size(average,1)/2)-3;%lowpass
    wei = 1./wei;rind = find(wei > 100000);wei(rind) = 0;% take care for inf
    average = real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(average)).*wei,lowp))));
    name = [refilename '_' num2str(ind+1) '.em'];
    tom_emwrite(name,average);
    disp(['wrote reference ' name]);
end;%for  - iterations
