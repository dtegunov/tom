function projstack = av3_projsfromparticles(motl,particlefilename,direction,nlayers,radius, threshold, ibin)
%
%   projstack = av3_projsfromparticles(motl,particlefilename,direction,nlayers,radius, threshold, ibin)
%

%tomo = tom_emread(tomofilename);
ipart = motl(4,1);
name = [particlefilename '_' num2str(ipart) '.em'];
particle = tom_emread(name);particle = tom_bin(particle.Value,ibin);
dims = [size(particle,1), size(particle,2), size(particle,3)];
if direction ==1
    projstack = zeros(dims(2),dims(3),size(motl,2));
elseif direction ==2
    projstack = zeros(dims(1),dims(3),size(motl,2));
elseif irection == 3
    projstack = zeros(dims(1),dims(2),size(motl,2));
else
    disp('choose direction = 1, 2 or 3 !');
    exit;
end;
indx = find (motl(1,:) > 0); meanv = mean(motl(1,indx));
indx = find (motl(1,:) > threshold*meanv);
icount = 0;
for ii=1:size(motl,2)
    if (motl(1,ii)>threshold*meanv)
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
        particle = tom_emread(name);particle = tom_bin(particle.Value,ibin);
        [mv, maxv, minv, stdv] = tom_dev(particle,'noinfo');particle=(particle-mv)/stdv;
        part = tom_spheremask(double(tom_rotate(tom_shift(particle,-tshift),[-psi,-phi,-theta])),radius,1);
        if direction ==1
            proj = sum(part(floor((dims(1)-nlayers)/2)+1:floor((dims(1)-nlayers)/2)+nlayers,:,:),direction);
        elseif direction ==2
            proj = sum(part(:,floor((dims(2)-nlayers)/2)+1:floor((dims(2)-nlayers)/2)+nlayers,:),direction);
        elseif direction == 3
            proj = sum(part(:,:,floor((dims(3)-nlayers)/2)+1:floor((dims(3)-nlayers)/2)+nlayers),direction);
        else
            disp('choose direction = 1, 2 or 3 !');
            exit;
        end;
        %proj = squeeze(proj);
        %[mv, maxv, minv, stdv] = tom_dev(proj,'noinfo');
        %proj = squeeze(proj-mv/st)
        projstack(:,:,ii) = squeeze(proj);tom_imagesc(squeeze(proj),'fixed');drawnow;
        %part = tom_red(vol.Value,[x-32,y-32,z-32],[64,64,64]);
    end;
end;
