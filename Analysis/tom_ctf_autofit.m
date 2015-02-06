function [e,n] = tom_ctf_autofit(directory, threshold, filterval,startval,split,binning,norings,lowcutoff)
% TOM_CTF_AUTOFIT searches for em images in a given directory and tries to
% do a CTF fit on each image. If successful the new defocus value is
% written to the header field "FocusIncrement", on failure the nominal
% value is copied.
%
% tom_ctf_autofit(directory, threshold, filterval)
% 
%
%  INPUT
% 
%  directory    The full path to the directory containing the images
%  threshold    The maximum difference between fitted and nominal defocus
%               in mu (optional) (default: 2)
%  filterval    The kernel size for a real space quadratic kernel to be
%               applied to each image before CTF fitting (optional) (default 1)
%  startval     start defocus value in nm (optional)
%
%  OUTPUT
%
%EXAMPLE
%   tom_ctf_autofit(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 09/11/06
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

error(nargchk(0, 8, nargin, 'struct'))

if nargin < 8
    lowcutoff = 100;
end

if nargin < 7
    norings = 5;
end

if nargin < 6
    binning = 0;
end

if nargin < 5
    split = 1;
end

if nargin < 4
    startval = [];
end

if nargin < 3
    filterval = 3;
end

if nargin < 2
    threshold = 2;
end

[dircell, flagcell_out] = get_dircontents(directory, {}, {});

e = zeros(length(dircell),5);
n = zeros(length(dircell),4);

for i=1:length(dircell)
    %try
    file = dircell{i};
    header = tom_reademheader(file);
    if ~isempty(startval)
        [Dz, success,envparams,noiseparams] = tom_ctffitter(file, startval, filterval, split, binning, 0, norings, lowcutoff);
    else
        [Dz, success,envparams,noiseparams] = tom_ctffitter(file, header.Header.Defocus./10, filterval, split, binning, 0, norings, lowcutoff);
    end

    e(i,:) = envparams;
    n(i,:) = noiseparams;
    
    %if success == 1 && abs(header.Header.Defocus-Dz) < threshold
        header.Header.FocusIncrement = Dz;
   %else
    %    header.Header.FocusIncrement = header.Header.Defocus;
   %end
   tom_writeemheader(file, header.Header);
   % end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get all the em-files in a directory                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dircell, flagcell_out] = get_dircontents(directory, dircell, flagcell)

%Get list of EM-Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Getting directory list...');

dirlist = dir(directory);
files = size(dirlist,1);
j = size(dircell,2);
if j == 0
    j = 1;
end

%sort out all the em-files
for i = 1:files
    if dirlist(i).isdir == 0
        if isempty(strmatch(dirlist(i).name,dircell))
            if tom_isemfile([directory '/' dirlist(i).name]) == 1
                dircell{j} = dirlist(i).name;
                try
                    flagcell_out{j} = flagcell{j};
                catch
                    flagcell_out{j} = 1;
                end
                j = size(dircell,2) + 1;
            end
        end
    end
    waitbar(i./files,h,[num2str(i), ' of ', num2str(files), ' files scanned.']);
end

close(h);

if size(dircell,2) == 0
    errordlg('No EM files could be found in this directory.');
    return;
end

set(findobj('Tag','fileslider'),'Max',size(dircell,2),'SliderStep',[1./size(dircell,2) 1./size(dircell,2)]);
