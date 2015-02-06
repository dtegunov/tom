function [cc,vol_list] = tom_av3_ccmat(particlefilename, appendix, npart, mask,hipass,lowpass,ibin,alg_param)
% AV3_CCMAT computes correlation- matrix of densities - designed for
% densities not affected by missing wedge etc. The primary use of this
% function is to compare desnities resulting from single particle
% classification.
%
%   cc = av3_ccmat(particlefilename, appendix, npart, mask,hipass,lowpass,ibin)
%
% PARAMETERS
%  INPUT
%   particlefilename    filename of 3D-particles to be correlated
%                           'particlefilename'_#no.em or
%                           'particlefilename'_#no.mrc
%   appendix            appendix for files - em or mrc or hdf
%   npart               number of densities (=last index)
%   mask                ('') mask - make sure dims are all right! ... use '' 2 swich off  
%   hipass              (-1) hipass - for X-corr use -1 to switch off
%   lowpass             (-1) lowpass - for X-corr use -1 to switch off
%   ibin                (0) binning 
%   alg_param           ('') search range  [-2 2 2; 2 2 2]  use '' 2 switch off
% 
%  OUTPUT
%   ccc                 constrained cross correlation matrix
%   vol_list            list of use volumes ...for further processing
%
% EXAMPLE
% mask = tom_spiderread('/fs/sandy03/lv04/pool/pool-titan1/CWT/combined/rec/projm_29_31__34_HiATPyS1mcl_39/output/models/mask_refine.spi');
% mask=mask.Value;
% cc = av3_ccmat('/fs/sandy03/lv04/pool/pool-titan1/CWT/combined/rec/projm_29_31__34_HiATPyS1mcl_39/output/models/align/model_aligned', ...
%        'mrc', 7, mask,0,100,2);
% tree = linkage(squareform(1-cc),'average');
% dendrogram(tree, 'ColorThreshold', .08);
%
% SEE ALSO
%   av3_cccmat, TOM_CORR, TOM_ORCD,tom_tree2chimera
%
%   Copyright (c) 2005-2012
%   Max-Planck-Institute for Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/foerster
%
%
if nargin < 4
    mask = '';
end;

if nargin < 5
    hipass = -1;
end;

if nargin < 6
    lowpass = -1;
end;

if nargin < 7
    ibin = 0;
end;

if ibin>0
    mask = tom_bin(mask,ibin)>0;
end;

if (nargin < 8)
    alg_param='';
end;


cc = zeros(npart,npart);
parfor indpart1 = 1:npart
    cc(indpart1,:)=corr_worker(indpart1,particlefilename, appendix, npart, mask,hipass,lowpass,ibin,alg_param);
    disp(['Correlation computed for particle no ' num2str(indpart1) ' ...']);
end;



%symmetrize matrices
for indpart1 = 1:npart
    cc(indpart1,indpart1) = 1;
    for kk = indpart1:npart-1
        indpart2 = kk + 1;
        cc(indpart2,indpart1) = cc(indpart1,indpart2);
    end;
end;

if (nargout > 1)
    for indpart1 = 1:npart
        vol_list{indpart1}=[particlefilename '' num2str(indpart1) '.' appendix];
    end;
end;

function cc=corr_worker(indpart1,particlefilename, appendix, npart, mask,hipass,lowpass,ibin,alg_param)

if strcmp(appendix,'em')
    name = [particlefilename '' num2str(indpart1) '.em'];
    part1 = tom_emread(name);
elseif strcmp(appendix,'mrc')
    name = [particlefilename '' num2str(indpart1) '.mrc'];
    part1 = tom_mrcread(name);
elseif strcmp(appendix,'hdf')
    name = [particlefilename '' num2str(indpart1) '.hdf'];
    part1 = tom_eman2_read(name);
else
    error('appendix must be em,hdf or mrc');
end;
if ibin > 0
    part1 = tom_bin(part1.Value,ibin);
else
    part1 = part1.Value;
end;
if (lowpass~=-1)
    part1 = tom_bandpass(part1,hipass,lowpass,5);
end;

if (isempty(mask))
    tmp1 = part1;
else
    tmp1 = part1.*mask;
end;

for indpart2 =indpart1:npart
    if strcmp(appendix,'em')
        name = [particlefilename '' num2str(indpart2) '.em'];
        part2 = tom_emread(name);
    elseif (strcmp(appendix,'em'))
        name = [particlefilename '' num2str(indpart2) '.mrc'];
        part2 = tom_mrcread(name);
    elseif (strcmp(appendix,'hdf'))    
        name = [particlefilename '' num2str(indpart2) '.hdf'];
        part2 = tom_eman2_read(name);
    end;
    
    
    if ibin > 0
        part2 = tom_bin(part2.Value,ibin);
    else
        part2 = part2.Value;
    end;
    if (lowpass~=-1)
        part2 = tom_bandpass(part2,hipass,lowpass,5);
    end;
    if (isempty(alg_param)==0)
        part2=tom_av3_align_xmipp(part1,part2,alg_param(1,:),alg_param(1,:),alg_param(1,:),alg_param(2,:),mask);
    end;
    if (isempty(mask))
        tmp2 = part2;
    else
        tmp2 = part2.*mask;
    end;
    if (isempty(mask))
        ccf=tom_corr(tmp1,tmp2,'norm');
    else
        ccf=tom_corr(tmp1,tmp2,'norm',mask);
    end;
    
    [pos val]=tom_peak(ccf,'spline');
    cc(1,indpart2) = val;
end;



