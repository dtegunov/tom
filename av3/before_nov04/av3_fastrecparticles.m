function av3_fastrecparticles(motl, projfile, append, first, last, particlename, voldims, ibin, iclass)
%
%   av3_fastrecparticles(motl, projfile, append, first, last, particlename, voldims, ibin, iclass)
%
%   motl            motive list
%   projfile        name of projection files - (weighted and aligned)
%                       assumed <projfile>no<append>
%   append          appendix of projection file
%   first           no of first projection
%   last            no of last projection
%   particlename    filename of particles
%   voldims         dimensions of particles
%   ibin            binning - voldims corresponfs to binned volume -
%                   default 0 - take care: the output volumes are binned in
%                   dimensions! - see example
%   iclass          class of particles to be reconstructed - if chosen to
%                   '0' all particles are reconstructed - default 0
%
%   example
%       av3_recparticles(motl, 'TEMP_BPP', '',1, 27, 'particle_dim128', [128 128 128],1,1);
%   reconstructs particles of class 1 with binning 1 (respective to
%   projections) from MOTL and projection files 'TEMP_BPP' 1 to 27. The
%   output files will be called 'particle_dim128_#num.em' and have
%   dimensions 128/(2^1)=64! The running number #num is read out from motl
%   (column 4). 
%
%   last change 29.03.04 FF
error(nargchk(7,9,nargin))
if nargin < 9
    iclass = 0;
end;
if nargin < 8
    ibin = 0;
end;
npartis = size(motl,2);
indx = motl(4,:);
if ibin > 0
    voldims = (2^ibin).*voldims;
end;
for iproj=first:last
    projname = [projfile  num2str(iproj) append];
    proj = tom_emread(projname);
    angle_the(iproj) = proj.Header.Tiltangle;
    projstack(:,:,iproj) = proj.Value;
end;
for irun=1:npartis
    ipart = indx(irun);        
    vol = single(zeros(voldims(1), voldims(2), voldims(3)));
    nclass = motl(20,irun);
    offset = [motl(8,irun) motl(9,irun) motl(10,irun)];
    %offset = (offset - [257 257 65])*4;% to be edited...
    if (iclass == 0) | (iclass == nclass) 
        for iproj=first:last
            tom_backproj3d(vol,single(projstack(:,:,iproj)), 0, angle_the(iproj), offset);
        end;
        pname = [particlename '_' num2str(ipart) '.em'];
        if ibin > 0
            vol = tom_bin(vol,ibin);
        end;
        tom_emwrite(pname, vol);
        disp(['particle ' pname ' written']);
    else
        disp(['particle not written']);
    end;
end;