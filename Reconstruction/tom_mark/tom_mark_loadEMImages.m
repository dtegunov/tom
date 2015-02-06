function varargout = tom_mark_loadEMImages(varargin)
% This is a helping routine for other functions to load an 
% image stack.
%
% Mode can be the following:
%   'init': Needs a second parameter filenames, an cell array with the 
%           filenames associated to the stack. The image-structure is
%           initalised and 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (nargin < 1)
    error('tom_mark_loadEMImages need parameters.');
end;


options = {'image', 'init', 'header'};
option = '';
if (~isstr(varargin{1}))
    option = options{1};
else
    for (i=1:length(options))
        if (strcmp(varargin{1}, options{i}))
            option = varargin{1};
            varargin = varargin(2:end);
            break;
        end;
    end;    
    if (isempty(option))
        error('Unknown option given');
    end;
end;

if (strcmp(option, 'init'))
    if (length(varargin)<1 || ~iscellstr(varargin{1}))
        error('init expects as second parameter a cell string with filenames');
    end;
    nfilenames = length(varargin{1});
    
    imstr = struct();
    imstr.filenames = varargin{1};
    imstr.usednames = cell(1, nfilenames);
    imstr.Header = cell(1, nfilenames);
    imstr.Value = cell(1,nfilenames);
    imstr.error = false(1, nfilenames);
    imstr.emread_binning = 0;
    imstr.transform = [];
    imstr.cacheorder = [];
    imstr.freesize = 0;
    
else
    if (length(varargin) < 2)
        error(['Not enough input arguments for "' option '". At least the image structure and the index is needed.']);
    end;
    
    imstr = varargin{1};
    idx = varargin{2};
    
    isoption_header = strcmp(option, 'header');
    if (~isoption_header)
        isoption_image = strcmp(option, 'image');
    end;
    
    for (i=idx)
        
        filename = imstr.filenames{i};
        
        if (~strcmp(filename, imstr.usednames{i}))
            imstr.usednames{i} = filename;
            imstr.Header{i} = [];
            imstr.Value{i} = [];
            imstr.cacheorder(imstr.cacheorder==i) = [];
        end;
        
        if (isoption_header)
            if (~imstr.error(i) && isempty(imstr.Header{i}))
                emheader.Header = [];
                for (iii=1:100)
                    try
                        emheader = tom_reademheader(filename);
                        break;
                    catch
                        warning(['Error reading image ' filename '. Retry...']);
                        emheader.Header = [];
                        imstr.error(i) = true;
                        pause(1);
                    end;
                end;
                imstr.Header{i} = emheader.Header;
            end;
            
        elseif (isoption_image)
            if (~imstr.error(i))
                imstr.cacheorder(imstr.cacheorder==i) = [];
                if (isempty(imstr.Value{i}))
                    %disp(['--------------']);
                    %disp(['CACHEORDER: ' num2str(imstr.cacheorder)]);
                    im = [];
                    for (iii=1:100)
                        try
                            im = tom_emreadc(filename,'binning', imstr.emread_binning);
                            %disp(['READ IMAGE #' num2str(i)]);
                            break;
                        catch
                            im = [];
                            warning(['Error reading image ' filename '. Retry...']);
                            pause(1);
                        end;
                    end;
                    if (isempty(im))
                        imstr.error(i) = true; 
                        imstr.Header{i} = [];
                        imstr.Value{i} = [];
                    else
                        imstr.Header{i} = im.Header;

                        im = im.Value';
                        if (isempty(imstr.transform))
                            imstr.Value{i} = im;
                        else
                            imstr.Value{i} = tom_mark_applyTransf('image_pre', im, imstr.transform);
                        end;
                        imstr.cacheorder(end+1) = i;
                        
                        
                        if (imstr.freesize > 0)
                            while (~isempty(imstr.cacheorder) && getfield(whos('imstr'), 'bytes')>imstr.freesize)
                                %disp(['FREE IMAGE' num2str(imstr.cacheorder(1))]);
                                imstr.Value{imstr.cacheorder(1)} = [];
                                imstr.cacheorder(1) = [];
                            end;
                        end;
                        
                    end;
                    %disp(['CACHEORDER: ' num2str(imstr.cacheorder)]);
                    
                else
                    imstr.cacheorder(end+1) = i;
                end;
            end;
        end;
    end;
    
end;


varargout{1} = imstr;



    