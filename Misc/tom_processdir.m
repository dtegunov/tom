function tom_processdir(varargin)
%TOM_PROCESSDIR takes files in input dir, processes and writes them to output dir
%
%   tom_processdir(varargin)
%
%   tom_processdir takes all files in input directory, processes them and
%   writes them to the output dir
%
%PARAMETERS
%
%  INPUT
%   INPUTDIR    input directory
%   OUTPUTDIR   output directory, will be created if it doesn't exist
%   FILETYPE    em, mrc or spi, all files will be saved as em
%   
%   The following parameters can be any combination of the methods:
%   'bin'
%   'norm'
%   'filter'
%   'taper'
%   'mirror'
%   'shift'
%   'xraycorrect'
%   for additional parameters look at TOM_BIN, TOM_NORM, TOM_FILTER, TOM_TAPER, TOM_MIRROR, TOM_SHIFT, TOM_XRAYCORRECT,
%   optional parameters can be left out.
%  
%  OUTPUT
%
%EXAMPLE
%   tom_processdir('/home/korinek/ptest/', '/home/korinek/ptest/processed','filter',5,'bin',1,'norm','3std');
%
%REFERENCES
%
%SEE ALSO
%   tom_processdirgui
%
%   created by AK 05/16/06
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


if nargin < 4
    error('Wrong number of arguments.')
end

indir = varargin{1};
outdir = varargin{2};
type = varargin{3};

if exist(indir,'dir') ~= 7
    error('Input directory does not exist!');
end

if exist(outdir,'dir') ~= 7
    disp('Output directory does not exist, creating new directory...');
    try
        mkdir(outdir);
    catch
        error('Could not create new directory!');
    end
end

h = waitbar(0,'Getting directory list...');

%generate list of files to be processed
filelist = get_dircontents(indir,type);

%read in files
files = size(filelist,2);
for j=1:files
    
    if (strcmp(type,'tif'))
        [aa bb cc]=fileparts([indir '/' filelist{j}]);
        new_name_out=[bb '.em'];
        im = imread([indir '/' filelist{j}]); 
        im = tom_emheader(single(im));
    else
        im = tom_emreadc([indir '/' filelist{j}]);
        new_name_out=filelist{j};
    end
    
    if size(im.Value,1) < 2 || size(im.Value,2) < 2
        disp('Skipping 1D data file.');
        continue;
    end
    im.Value = single(im.Value);
    k = 4;
    
    while k<=nargin
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'bin'
                    waitbar(j./files,h,['binning ... (' num2str(j), ' of ', num2str(files), ' files)']);
                    im.Value = tom_bin(im.Value, varargin{k+1});
                    im.Header.Size = size(im.Value);
                    im.Header.Objectpixelsize = im.Header.Objectpixelsize .* 2^varargin{k+1};
                    k = k + 2;
                    
                case 'norm'
                    waitbar(j./files,h,['norming ... (' num2str(j), ' of ', num2str(files), ' files)']);
                    if k+2 <= nargin && ~ischar(varargin{k+2})
                        im.Value = tom_norm(im.Value, varargin{k+1}, varargin{k+2});                        
                        k = k + 3;
                    else
                        im.Value = tom_norm(im.Value, varargin{k+1});
                        k = k + 2;
                    end
                    
                case 'filter'
                    waitbar(j./files,h,['filtering ... (' num2str(j), ' of ', num2str(files), ' files)']);
                    %anisotropic 
                    if k+1 <= nargin && ischar(varargin{k+1})
                        if k+2 < nargin && ischar(varargin{k+2})
                            im.Value = tom_filter(im.Value, varargin{k+1});
                             k = k + 2;
                        elseif k+3 <= nargin && ischar(varargin{k+3})
                            im.Value = tom_filter(im.Value, varargin{k+1}, varargin{k+2});
                            k = k + 3;
                        elseif k+4 <= nargin && ischar(varargin{k+4})
                            im.Value = tom_filter(im.Value, varargin{k+1}, varargin{k+2}, varargin{k+3});
                            k = k + 4;
                        else
                            im.Value = tom_filter(im.Value, varargin{k+1}, varargin{k+2}, varargin{k+3}, varargin{k+4});
                            k = k + 5;
                        end

                    %fourier/real space
                    else
                        if strcmp(varargin{k+2},'circ') == 0 && strcmp(varargin{k+2},'quadr') == 0
                            im.Value = tom_filter(im.Value, varargin{k+1});
                            k = k + 2;
                        elseif strcmp(varargin{k+3},'real') == 0 && strcmp(varargin{k+3},'fourier') == 0
                            im.Value = tom_filter(im.Value, varargin{k+1}, varargin{k+2});
                            k = k + 3;
                        else
                            im.Value = tom_filter(im.Value, varargin{k+1}, varargin{k+2}, varargin{k+3});
                            k = k + 4;
                        end

                    end

                case 'taper'
                    waitbar(j./files,h,['tapering ... (' num2str(j), ' of ', num2str(files), ' files)']);
                    im.Value = tom_taper(im.Value, varargin{k+1});
                    k = k + 2;
                
                case 'mirror'
                    waitbar(j./files,h,['mirroring ... (' num2str(j), ' of ', num2str(files), ' files)']);
                    im.Value = tom_mirror(im.Value, varargin{k+1});
                    k = k + 2;
                    
                case 'shift'
                    waitbar(j./files,h,['shifting ... (' num2str(i), ' of ', num2str(files), ' files)']);
                    im.Value = tom_shift(im.Value, varargin{k+1});
                    k = k + 2;
                
                case 'xraycorrect'
                    waitbar(j./files,h,['correcting x rays ... (' num2str(j), ' of ', num2str(files), ' files)']);
                    if k+1 > nargin || ischar(varargin{k+1})
                        im.Value = tom_xraycorrect(im.Value);
                        k = k + 1;
                    else
                        im.Value = tom_xraycorrect(im.Value, varargin{k+1});
                        k = k + 2;
                    end
                case 'editheader'
                    waitbar(j./files,h,['editing header ... (' num2str(i), ' of ', num2str(files), ' files)']);
                    im.Header.(varargin{k+1}) = varargin{k+2};
                    k = k + 3;
                otherwise
                    close(h);
                    error(['Method ' varargin{k} ' not implemented.']);
            end
        end
    end

    waitbar(j./files,h,['writing ... (' num2str(j), ' of ', num2str(files), ' files)']);
    tom_emwrite([outdir '/' new_name_out], im);
    
end
close(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get all the em-files in a directory                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dircell] = get_dircontents(directory,type)

dircell = {};
%Get list of EM-Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h = waitbar(0,'Getting directory list...');

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
            if isequal(type,'em')
                cf = tom_isemfile([directory '/' dirlist(i).name]);
            elseif isequal(type,'mrc')
               cf = tom_ismrcfile([directory '/' dirlist(i).name]);
            elseif isequal(type,'spi')
               cf = tom_isspiderfile([directory '/' dirlist(i).name]);
            elseif isequal(type,'tif')
                try
                    t_tmp = imfinfo([directory '/' dirlist(i).name]);
                    
                    if (strcmp(t_tmp.Format,'tif'))
                        cf=1;
                    else
                        cf=0;
                    end;
                catch Me
                    cf=0;
                end;
            else 
                error('file type unknown')
            end

            if  cf == 1
                dircell{j} = dirlist(i).name;
                j = size(dircell,2) + 1;
            end
        end
    end
%    waitbar(i./files,h,[num2str(i), ' of ', num2str(files), ' files scanned.']);
end

%close(h);

if size(dircell,2) == 0
    errordlg('No EM files could be found in this directory.');
    return;
end
