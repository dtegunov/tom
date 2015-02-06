function av3_checkorientations(motl,particlefilename,ifilt, ibin, iclass)
% AV3_CHECKORIENTATIONS displays aligned particles
%
%   av3_checkorientations(motl,particlefilename,ifilt, ibin, iclass)
%
%   Aligned particles are displayed in different TOM_DSPCUB views (Figure 1
%   '0' mode, fig.2 '1' mode and fig.3 '2' mode.
%
% PARAMETERS
%   MOTL        motive list
%   PARTICLE... filename of particle
%   IFILT       FILTER radius - to enhance visibility
%   IBIN        binning
%   ICLASS      class of particles to be displayed 
%
% SEE ALSO 
%   TOM_DSPCUB
%
% last change 03/31/05 FF - updated docu

%tomo = tom_emread(tomofilename);
iclass = 1;
ipart = motl(4,1);
name = [particlefilename '_' num2str(ipart) '.em'];
particle = tom_emread(name);
if ibin > 0
    particle = tom_bin(particle.Value,ibin);
else
    particle = particle.Value;
end;
dims = [size(particle,1), size(particle,2), size(particle,3)];
icount = 0;
for ii=1:size(motl,2)
    icount = icount +1;
    if (motl(20,ii) == iclass)
        xshift = motl(11,ii)/(2^ibin);
        yshift = motl(12,ii)/(2^ibin);
        zshift = motl(13,ii)/(2^ibin);
        tshift = [xshift yshift zshift];
        phi=motl(17,ii);
        psi=motl(18,ii);
        theta=motl(19,ii);
        ifile = motl(4,ii);
        name = [particlefilename '_' num2str(ifile) '.em'];
        particle = tom_emread(name);
        if ibin > 0
            particle = tom_bin(particle.Value,ibin);
        else
            particle = particle.Value;
        end;
        particle = tom_filter(particle,ifilt,'circ');
        [mv, maxv, minv, stdv] = tom_dev(particle,'noinfo');
        particle= ((particle-mv)/stdv);
        particle = double(tom_rotate(tom_shift(particle,-tshift),[-psi,-phi,-theta]));
        figure(1);tom_dspcub(tom_limit(particle,-3,3));title(['particle ' num2str(ifile) ' top']);drawnow;
        figure(2);tom_dspcub(tom_limit(particle,-3,3),1);title(['particle ' num2str(ifile) ' side']);drawnow;
        figure(3);tom_dspcub(tom_limit(particle,-3,3),2);title(['particle ' num2str(ifile) ' side 2']);drawnow;
    end;
end;
