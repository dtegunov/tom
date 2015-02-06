function av3_checkparti(motl,particlefilename,ifilt, ibin,ifile)
%
%   av3_checkparti(motl,particlefilename,ifilt, ibin,ifile)
%
%   MOTL        motive list
%   PARTICLE... filename of particle
%   IFILT       FILTER radius
%   RADIUS        mask applied BEFORE calculating MOI
%   IBIN        binning
%   threshold   thresholding for goldbeads - if min value smaller than
%               threshold
%   RMS         rms - if samller than signal too weak
%
%   idea of av3_mois is to create a list with moments of inertita to
%   perform classification based on the anisotropy of moments of inertia -
%   the procedure is designed for 26S project to distinguish 20S from 26S
%   proteasomes. 
%   Format of moistack:
%   stack of dim 5 x #tomos
%   colums
%   1: index of file
%   2: CCC (from MOTL)
%   3: smallest ev
%   4: middle ev
%   5: highest ev
%   6: class asiigned to 1 for 26S
%
%   last change
%   11/05/04 FF

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
ii = find(motl(4,:) == ifile);
icount = icount +1;
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
figure(1);tom_dspcub(tom_limit(particle,-3,3),1);title(['particle ' num2str(ifile) ' side']);drawnow;
figure(2);tom_dspcub(tom_limit(particle,-3,3));title(['particle ' num2str(ifile) ' top']);drawnow;
figure(3);tom_dspcub(tom_limit(particle,-3,3),2);title(['particle ' num2str(ifile) ' side 2']);drawnow;
