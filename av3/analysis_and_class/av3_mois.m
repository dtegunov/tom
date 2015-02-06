function moistack = av3_mois(motl,particlefilename,ifilt, radius, ibin)
% AV3_MOIS computes principal moments of inertia of particles
%
%   moistack = av3_mois(motl,particlefilename,ifilt, radius, ibin)
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
%
% last change 03/31/05 FF - update docu

%tomo = tom_emread(tomofilename);
ipart = motl(4,1);
name = [particlefilename '_' num2str(ipart) '.em'];
particle = tom_emread(name);
if ibin > 0
    particle = tom_bin(particle.Value,ibin);
else
    particle = particle.Value;
end;
dims = [size(particle,1), size(particle,2), size(particle,3)];
moistack = zeros(5,size(motl,2));
icount = 0;
for ii=1:size(motl,2)
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
    [mv, maxv, minv, stdv] = tom_dev(particle,'noinfo');particle=-((particle-mv)/stdv);
    particle = tom_spheremask(double(tom_rotate(tom_shift(particle,-tshift),[-psi,-phi,-theta])),radius,3);
    indx = find(particle < 0.4); particle(indx) = 0;tom_dspcub(particle,1);drawnow;
    [toi,eigvec,eigval] = tom_moi(particle);
    moistack(1,ii) = motl(4,ii);moistack(2,ii) = motl(1,ii);
    moistack(3,ii) = eigval(1,1);moistack(4,ii) = eigval(2,2);moistack(5,ii) = eigval(3,3);
    if (moistack(5,ii)./(moistack(3,ii)+moistack(4,ii))) > 1.2
        title('26S ? !');drawnow;
    end;
end;
