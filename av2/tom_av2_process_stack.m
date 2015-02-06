function stack_new=tom_av2_process_stack(varargin)
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
%   created by fb(eckster) 05/20/08
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


if ischar(varargin{1});
    stack=tom_emread(varargin{1});
else
    stack=varargin{1};
    stack=tom_emheader(stack);
end;


h = waitbar(0,'Getting directory list...');
files=size(stack.Value,3);
im=tom_emheader(stack.Value(:,:,1));

for j=1:size(stack.Value,3)
    
    im.Value=stack.Value(:,:,j);
    k=3;

    while k<=nargin
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'bin'
                    waitbar(j./files,h,['binning ... (' num2str(j), ' of ', num2str(files), ' files)']);
                    im.Value = tom_bin(im.Value, varargin{k+1});
                    im.Header.Size = size(im.Value);
                    im.Header.Objectpixelsize = im.Header.Objectpixelsize ./ 2^varargin{k+1};
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
                        if nargin < 5
                            im.Value = tom_filter(im.Value, varargin{k+1});
                            k = k + 2;
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
                        end;
                        
                        
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
    k=k+1;
    end

    stack_new(:,:,j)=im.Value;
    
end;

if isempty(varargin{2})==0 
    tom_emwrite(varargin{2},stack_new);
end;

close(h);


