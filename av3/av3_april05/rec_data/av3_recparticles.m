function av3_recparticles(motl, projfile, append, first, last, particlename, voldims, ibin, iclass)
% AV3_RECPARTICLES reconstructs subtomograms according to MOTL
%
%   av3_recparticles(motl, projfile, append, first, last, particlename, voldims, ibin, iclass)
%
%   AV3_RECPARTICLES is used to reconstruct individual subtomograms
%   incorporating particles based on a MOTL and weighted projections
%   specified by PROJFILE, their APPENDix (e.g. '.em'), and their running
%   numbers (FIRST and LAST need to be specified). The particles are
%   reconstructed as PARTICLEFILENAME_#.em, their size is determined by
%   VOLDIMS. The particles can be binned respective to the weighted
%   projections and a class ICLASS can be specified.
%
% PARAMETERS
%   motl            motive list
%   projfile        name of projection files - (weighted and aligned)
%                       assumed <projfile>no
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
%  EXAMPLE
%
%   av3_recparticles(motl, 'TEMP_BPP', '', 1, 27, 'particle_dim128', ...
%               [128 128 128],1,1);
%
%   reconstructs particles of class 1 with binning 1 (respective to
%   projections) from MOTL and projection files 'TEMP_BPP' 1 to 27. The
%   output files will be called 'particle_dim128_#num.em' and have
%   dimensions 128/(2^1)=64! The running number #num is read out from motl
%   (column 4). 
%
% SEE ALSO
%   AV3_FASTRECPARTICLES, AV3_CREATEMOTL, AV3_TRANS-ROT-ALIG
%
%   last change 03/31/05 FF - updated docu
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
for irun=1:npartis
    ipart = indx(irun);        
    vol = single(zeros(voldims(1), voldims(2), voldims(3)));
    nclass = motl(20,irun);
    offset = [motl(8,irun) motl(9,irun) motl(10,irun)];
    %offset = (offset - [257 257 65])*4;% to be edited...
    if (iclass == 0) | (iclass == nclass) 
        for iproj=first:last
            projname = [projfile  num2str(iproj) append ]; %edited by Julio Oct 2009: add append
            proj = tom_emread(projname);
            angle_the = proj.Header.Tiltangle;
            tom_backproj3d(vol,single(proj.Value), 0, angle_the, offset);
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
