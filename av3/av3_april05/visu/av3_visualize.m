function av3_visualize(motl, particlefilename, iclass, mode, ifilt, delta)

%AV3_VISUALIZE  visualizes 3D-particles
%
%USAGE
%   av3_visualize(motivelist, filename, iclass, mode, ifilt, delta);
%
%   MOTL                Motive list
%   PARTICLEFILENAME    Filename of particles
%   ICLASS              Class index
%   MODE                Display mode - 'x','y', or 'z'
%   IFILT               Kernel size of filer
%   DELTA               Number of slices - cent-delta:cent+delta are
%                       averaged
%
%   Particles are rotated and shifted according to a motl.
%   last change 11/05/04 FF

for indpart = 1:size(motl,2)
    %if ((motl(20,indpart) == 1) | (motl(20,indpart) == 2) | (motl(20,indpart) == iclass))
    if ((motl(20,indpart) == iclass))
        ifile = motl(4,indpart);
        name = [particlefilename '_' num2str(ifile) '.em'];
        particle = tom_emread(name);
                tshift(1) = motl(11,indpart);
        tshift(2) = motl(12,indpart);
        tshift(3) = motl(13,indpart);
        phi_opt=motl(17,indpart);
        psi_opt=motl(18,indpart);
        the_opt=motl(19,indpart);
        particle  = tom_filter(double(tom_rotate(tom_shift(particle.Value,-tshift),[-psi_opt,-phi_opt,-the_opt])),ifilt);
        cent = floor(size(particle,1)/2)+1;
        switch lower(mode)
            case 'x'
                if (delta > 0) 
                    particle = tom_spheremask(squeeze(sum(particle(cent-delta:cent+delta,:,:),1)),cent-3,1);
                else
                    particle = tom_spheremask(squeeze(particle(cent,:,:)));
                end;
            case 'y'
                if (delta > 0) 
                    particle = tom_spheremask(squeeze(sum(particle(:,cent-delta:cent+delta,:),2)),cent-3,1);
                else
                    particle = tom_spheremask(squeeze(particle(:,cent,:)));
                end;
            case 'z'
                if (delta > 0) 
                    particle = tom_spheremask(squeeze(sum(particle(:,:,cent-delta:cent+delta),3)),cent-3,1);
                else
                    particle = tom_spheremask(squeeze(particle(:,:,cent)));
                end;
        end;%switch
        [mv, dum1, dum2, stdv] = tom_dev(particle,'on');
        tom_imagesc(particle,'fixed');drawnow;
        particle = tom_limit(particle-mv,-3*stdv,3*stdv);
    end;%if
end;%for - motl_indx
