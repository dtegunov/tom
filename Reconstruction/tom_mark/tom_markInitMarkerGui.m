function varargout = tom_markInitMarkerGui(varargin)
% TOM_MARKINITMARKERGUI M-file for tom_markInitMarkerGui.fig
%      TOM_MARKINITMARKERGUI, by itself, creates a new TOM_MARKINITMARKERGUI or raises the existing
%      singleton*.
%
%      H = TOM_MARKINITMARKERGUI returns the handle to a new TOM_MARKINITMARKERGUI or the handle to
%      the existing singleton*.
%
%      TOM_MARKINITMARKERGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_MARKINITMARKERGUI.M with the given input arguments.
%
%      TOM_MARKINITMARKERGUI('Property','Value',...) creates a new TOM_MARKINITMARKERGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_markInitMarkerGui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_markInitMarkerGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_markInitMarkerGui

% Last Modified by GUIDE v2.5 25-Apr-2007 15:06:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_markInitMarkerGui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_markInitMarkerGui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before tom_markInitMarkerGui is made visible.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tom_markInitMarkerGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_markInitMarkerGui (see VARARGIN)

if (length(varargin) < 1)
    error('markInitMarkerGui needs as first parameter the filenames of the images or it image size (2-vector).');
else
    if (ischar(varargin{1}))
        handles.data.filenames = {varargin{1}};
    elseif (iscellstr(varargin{1}) && length(varargin{1})>=1)        
        handles.data.filenames = varargin{1};
    else
        error('markInitMarkerGui needs as first parameter the filenames of the images.');
    end;
end;

if (length(varargin) >= 2 && isstruct(varargin{2}) && numel(varargin{2})==1)
    handles.data.config0 = varargin{2};    
else
    if (length(varargin) >=2 && ~isempty(varargin{2}))
        error('markInitMarkerGui needs as second parameter the configuration structure, an empty matrix or nothing.');
    end;
    handles.data.config0 = [];
end;
idx_params = 3;

wait = true;

handles.data.format = 'xy';
handles.data.assignname = '';

while (idx_params <= length(varargin))
    param_recognized = false;
    errors = ['Unrecognized parameter/option'];
    if (ischar(varargin{idx_params}))
        switch (varargin{idx_params})
            case 'nowait'
                wait = false;
                param_recognized = true;
            case 'assignname'
                if (idx_params < length(varargin))
                    idx_params = idx_params + 1;
                    if (isstr(varargin{idx_params}))
                        handles.data.assignname = varargin{idx_params};
                        param_recognized = true;
                    else
                        errors = 'Option ''assignname'' needs a variable name as second parameter.';
                    end;
                else
                    errors = 'Option ''assignname'' needs a second parameter.';
                end;
            case 'wait'
                wait = true;
                param_recognized = true;
            case 'format'
                if (idx_params < length(varargin))
                    idx_params = idx_params + 1;
                    switch (varargin{idx_params})
                        case {'xy','axy','axy12','em','-xy','-axy','-axy12','-em'}
                            handles.data.format = varargin{idx_params};
                            param_recognized = true;
                        otherwise
                            errors = ['Unrecognized format ''' varargin{idx_params} ''''];
                    end;    
                else
                    errors = 'Option ''format'' needs a second parameter.';
                end;
            otherwise
                if (idx_params < length(varargin))
                    set(handles.figure_markInitMarkerGui, varargin{idx_params}, varargin{idx_params+1});
                    idx_params = idx_params + 1;
                    param_recognized = true;
                end;
        end;
    end;
    if (~param_recognized)
        error(errors);
    end;
    idx_params = idx_params + 1;
end;


handles.data.emheader = tom_markGui('callLocalFcn', handles.figure_markInitMarkerGui, 'getFileInfo_emheader', handles.data.filenames);
if (isempty([handles.data.emheader{:}]))
    error('No images files found!');
end;

handles = getFileInfo_info(handles);


handles = loadConfig(handles, handles.data.config0);


handles = loadGui(handles);


guidata(hObject, handles);

if (wait)
    uiwait(handles.figure_markInitMarkerGui);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function waitbar_cancelclick(hObject, eventdata, waitbarh)

%set(waitbarh, 'Userdata', true);
1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs from this function are returned to the command line. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_markInitMarkerGui_OutputFcn(hObject, eventdata, handles) 

if (~isfield(handles.data, 'exitstatus'))
    handles.data.exitstatus = '';
end;

switch (handles.data.exitstatus)
    case {'pushbutton_ok', 'pushbutton_ok_nowait'}
        
        handles = unloadGui(handles);

        delete(findall(0, 'type', 'figure', 'Tag', 'TMWWaitbar_forGetMarkerSet'));
        waitbarh = waitbar(0, 'initializing markerset...', 'CreateCancelBtn', 'set(gcf, ''Userdata'', true);', 'Userdata', false, 'Tag', 'TMWWaitbar_forGetMarkerSet');
        %set(waitbarh, 'WindowStyle', 'modal');
        set(waitbarh, 'WindowStyle', 'normal');
        ms = getMarkerset(handles, waitbarh);
        close(waitbarh);
        delete(waitbarh);
        
        msize = size(ms);
        if (length(msize)==2)
            msize(3) = 1;
        end;
        
        if (~any(strcmp({'xy','axy', 'axy12','em', '-xy','-axy', '-axy12','-em'}, handles.data.format)))
            warning(['formattype not found!: UNEXPECTED ' handles.data.format]);
            handles.data.format = 'xy';
        end;

        if (handles.data.format(1) == '-')
            ms(~(ms > 0)) = -1;
        end;
        
        switch (handles.data.format)
            case {'axy', '-axy'}
                outputms = nan( 3,msize(2),msize(3));
                outputms([2 3], :, :) = ms;
            case {'axy12','em', '-xy','-axy', '-axy12','-em'}
                outputms = nan(12,msize(2),msize(3));
                outputms([2 3], :, :) = ms;
            case {'xy', '-xy'}
                outputms = ms;
            otherwise
                error('hallo');
        end;
        if (size(outputms,1) > 2 && msize(3)>0)
            for (i=1:msize(2))
                if (isempty(handles.data.emheader{i}))
                    outputms(1,i,1) = nan;
                else
                    outputms(1,i,1) = handles.data.emheader{i}.Tiltangle;
                end;
            end;
        end;
        
        if (any(strcmp({'em','-em'}, handles.data.format)))
            outputms = tom_emheader(outputms);
        end
        
        if (strcmp(handles.data.exitstatus, 'pushbutton_ok'))
            if (nargout >= 1)
                varargout{1} = outputms;
            end;
            if (nargout >= 2)
                varargout{2} = handles.data.config;
            end;
        end;
        if (isfield(handles.data, 'assignname') && ~isempty(handles.data.assignname))
            assignin('base', handles.data.assignname, outputms);
        end;
        delete(handles.figure_markInitMarkerGui);
    case {'pushbutton_cancel', 'pushbutton_cancel_nowait'}
        if (strcmp(handles.data.exitstatus, 'pushbutton_cancel'))
            if (nargout >= 1)
                varargout{1} = false;
            end;
            if (nargout >= 2)
                varargout{2} = handles.data.config0;
            end;
        end;
        delete(handles.figure_markInitMarkerGui);
    otherwise
        if (~isfield(handles.data, 'assignname') || isempty(handles.data.assignname))
            warning('gui startet with parameter ''nowait''. Use the parameter ''assignname'' to copy the markerset in the workspace');
        end;
        if (nargout > 0)
            varargout = cell(1,nargout);
        end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% taken and modified from
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=14331
% Draw a circle in a matrix using the integer midpoint circle algorithm
% Does not miss or repeat pixels
% Created by : Peter Bone
% Created : 19th March 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cir = getCircleMask(radius, vv)

if (radius <= 0.5)
    cir = true;
else
    dim = ceil(radius)+1;
    radius2 = radius^2;
    cir4 = true(dim);
    for (i=1:dim)
        for (j=i:dim)
            if ((i-1)^2+(j-1)^2 > radius2)
                cir4(i,j) =false;
                cir4(j,i) =false;
            end;
        end;
    end;
    cir = false(dim*2-1);
    cir(dim:-1:1,dim:-1:1) = cir4;
    cir(dim:-1:1,dim+(1:dim-1)) = cir4(:,2:end);
    cir(dim+(dim-1:-1:1),1:(2*dim-1)) = cir(1:dim-1,:);
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ms = getMarkerset(handles, waitbarh)

if (~exist('waitbarh', 'var'))
    waitbarh = [];
end;



waitbar_abbort = false;

switch (handles.data.config.method)
    case 'random'
        msize = [2, length(handles.data.filenames), 0];
        ms = nan(msize(1),msize(2),length(handles.data.config.sel_images)*handles.data.config.nmarkers);
        
        for (iim=1:length(handles.data.config.sel_images))
            sel_image = handles.data.config.sel_images(iim);
            msize3_old = msize(3);
            waitbar_cnt_tot = msize3_old + (length(handles.data.config.sel_images)-iim+1) * handles.data.config.nmarkers;
            
            if (~isempty(waitbarh))
                waitbar((msize3_old+0)/waitbar_cnt_tot, waitbarh);
            end;
            
            if (~isempty(handles.data.emheader{sel_image}))
                %disp(['Finding ' num2str(handles.data.config.nmarkers) ' random marker points in image ' num2str(sel_image)]);
                imsize = handles.data.emheader{sel_image}.Size([2 1])';
                if (handles.data.config.usemask)
                    mask = logical(handles.data.config.mask{sel_image});
                else
                    mask = [];
                end;
                if (isempty(mask))
                    mask = true(imsize);
                elseif (any(size(mask) ~= imsize))
                    mask = imresize(mask, imsize, 'bilinear');
                end;
                br = 0;
                if (handles.data.config.omitborder)
                    br = handles.data.config.omitborderpixel;
                    if (br < 1)
                        br = br * min(imsize);
                    end;
                    br = round(br);
                    if (~(br >= 0))
                        br = 0;
                    end;
                end;
                if (br > 0)
                    %disp(['   Border: ' num2str(br) ' pixel']);
                    mask = mask((br+1):(end-br), (br+1):(end-br));
                end;
                imsize = size(mask);
                
                if (handles.data.config.mindist_auto)
                    mindist = max(0,round(0.75*sqrt(sum(mask(:))/handles.data.config.nmarkers)));
                else
                    mindist = handles.data.config.mindist;
                end;
                %disp(['   Mindist: ' num2str(mindist) ' pixel']);
                
                clear_mask = getCircleMask(mindist);

                
                [idxy, idxx] = find(mask);
                idxx = uint32(idxx);
                idxy = uint32(idxy);
                idxinv = zeros(size(mask), 'uint32');
                lidx = length(idxx);

                if (lidx > 0)
                    idxa = true(1,lidx);

                    for (iix=1:lidx)
                        idxinv(idxy(iix), idxx(iix)) = iix;
                    end;


                    clear_mask_size = size(clear_mask, 1);
                    clear_mask_centre = fix(size(clear_mask,1)/2)+1;
                    relevant_idx = 1:lidx;
                    for (ima=1:handles.data.config.nmarkers)
                        if (mod(ima,100) && ~isempty(waitbarh))
                            waitbar((msize(3))/waitbar_cnt_tot, waitbarh);
                            waitbar_abbort = get(waitbarh, 'Userdata');
                            if (waitbar_abbort)
                                break;
                            end;
                        end;
                        for (j=1:30)
                            relevant_idxi = ceil(rand()*length(relevant_idx));
                            mpointi = relevant_idx(relevant_idxi);
                            if (idxa(mpointi))
                                break;
                            end;
                        end;
                        if (~idxa(mpointi))
                            relevant_idx = find(idxa);
                            if (length(relevant_idx) < 1)
                                break;
                            end;
                            relevant_idxi = ceil(rand()*length(relevant_idx));
                            mpointi = relevant_idx(relevant_idxi);
                        end;

                        mpoint = double([idxx(mpointi), idxy(mpointi)]);
                        msize(3) = msize(3)+1;
                        ms([1 2], sel_image, msize(3)) = mpoint;                    

                        if (clear_mask_centre == 1)
                            idxa(mpointi) = false;
                            idxinv(mpoint(2), mpoint(1)) = 0;
                        else
                            range_maskx = max(1, mpoint(1)-clear_mask_centre+1) : min(imsize(2), mpoint(1)+clear_mask_centre-1);
                            range_masky = max(1, mpoint(2)-clear_mask_centre+1) : min(imsize(1), mpoint(2)+clear_mask_centre-1);
                            range_clear_maskx = max(1, clear_mask_centre-mpoint(1)+1) : min(clear_mask_size, clear_mask_centre+imsize(2)-mpoint(1));
                            range_clear_masky = max(1, clear_mask_centre-mpoint(2)+1) : min(clear_mask_size, clear_mask_centre+imsize(1)-mpoint(2));

                            for (iix = 1:length(range_maskx))
                                for (iiy = 1:length(range_masky))
                                    if (clear_mask(range_clear_masky(iiy), range_clear_maskx(iix)))
                                        idxii = idxinv(range_masky(iiy), range_maskx(iix));
                                        if (idxii ~= 0)
                                            idxa(idxii) = false;
                                            idxinv(range_masky(iiy), range_maskx(iix)) = 0;
                                        end;
                                    end;
                                end;
                            end;

                            %figure(2);hold off;imagesc(idxinv3);axis off equal;colormap gray;hold on;plot(squeeze(ms(1,sel_image,:)),squeeze(ms(2,sel_image,:)),'g+');drawnow;
                        end;
                    end;
                    if (br > 0) 
                        ms([1 2],sel_image,msize3_old+1:msize(3)) = ms([1 2],sel_image,msize3_old+1:msize(3)) + br;
                    end;
                    if (waitbar_abbort)
                        break;
                    end;
                end;
            end;
        end;
        ms = ms([1 2],1:msize(2),1:msize(3));
        ms(~(ms>0)) = nan;

    case 'grid'
        
        ms = {};
        for (iim=1:length(handles.data.config.sel_images))
            sel_image = handles.data.config.sel_images(iim);
            
            if (~isempty(waitbarh))
                if (get(waitbarh, 'Userdata'))
                    break;
                end;
                waitbar((iim-1)/length(handles.data.config.sel_images), waitbarh);
            end;
            
            if (~isempty(handles.data.emheader{sel_image}))
                %disp(['Finding ' num2str(handles.data.config.nmarkers) ' random marker points in image ' num2str(sel_image)]);
                imsize = handles.data.emheader{sel_image}.Size([2 1])';
                if (handles.data.config.usemask)
                    mask = logical(handles.data.config.mask{sel_image});
                else
                    mask = [];
                end;
                if (isempty(mask))
                    mask = true(imsize);
                elseif (any(size(mask) ~= imsize))
                    mask = imresize(mask, imsize, 'bilinear');
                end;
                br = 0;
                if (handles.data.config.omitborder)
                    br = handles.data.config.omitborderpixel;
                    if (br < 1)
                        br = br * min(imsize);
                    end;
                    br = round(br);
                    if (~(br >= 0))
                        br = 0;
                    end;
                end;
                if (br > 0)
                    %disp(['   Border: ' num2str(br) ' pixel']);
                    mask = mask((br+1):(end-br), (br+1):(end-br));
                end;
                imsize = size(mask);
                
                if (handles.data.config.mindist_auto)
                    mindist = max(0,round(sqrt(sum(mask(:))/handles.data.config.nmarkers)));
                else
                    mindist = handles.data.config.mindist;
                end;
                %disp(['   Mindist: ' num2str(mindist) ' pixel']);
                
                left_upper_corner = fix(mod(imsize, mindist)/2) + 1;
               
                msize = [2, length(handles.data.filenames), 0];
                msi = nan(2, msize(2), prod(ceil(imsize / (mindist-1))));

                for (ix = left_upper_corner(1):mindist:imsize(1))
                    for (iy = left_upper_corner(2):mindist:imsize(2))
                        if (mask(iy, ix))
                            msize(3) = msize(3)+1;
                            msi(1:2, sel_image, msize(3)) = [ix, iy];
                        end;
                    end;
                end;
                
                if (br > 0) 
                    ms{iim} = msi([1 2],:,1:msize(3)) + br;
                else
                    ms{iim} = msi([1 2],:,1:msize(3));
                end;
            end;
            
        end;
        
        i = 0;
        for (iim = 1:length(ms))
            i = i + size(ms{iim}, 3);
        end;        
        mms = nan(2, length(handles.data.filenames), i);
        i = 0;
        for (iim = 1:length(ms))
            mms(1:2, :, (i+1):(i+size(ms{iim}, 3))) = ms{iim};
            i = i + size(ms{iim}, 3);
        end;        
        ms = mms;

    case 'contrast';
        msize = [2, length(handles.data.filenames), 0];
        ms = nan(msize(1),msize(2),length(handles.data.config.sel_images)*handles.data.config.nmarkers);

        cma = handles.data.config.nmarkers;
        for (iim=1:length(handles.data.config.sel_images))
            msize3_old = msize(3);

            sel_image = handles.data.config.sel_images(iim);

            waitbar_cnt_tot = msize3_old + (length(handles.data.config.sel_images)-iim+1) * handles.data.config.nmarkers;
            
            if (~isempty(waitbarh))
                waitbar((msize3_old+0)/waitbar_cnt_tot, waitbarh);
            end;
            
            
            % READING THE IMAGE
            try
                im = tom_emreadc(handles.data.filenames{sel_image}, 'binning', handles.data.config.method_contrast_binning);
            catch
                im = [];
            end;
            if (isempty(im))
                continue;
            end;
            im = tom_mark_applyTransf('image_pre', im.Value', handles.data.config.method_contrast_filter);
            imsize = size(im);
            

            % READING THE MASK and resize.
            if (handles.data.config.usemask)
                mask = logical(handles.data.config.mask{sel_image});
            else
                mask = [];
            end;
            if (isempty(mask))
                mask = true(imsize);
            elseif (any(size(mask) ~= imsize))
                mask = imresize(mask, imsize, 'bilinear');
            end;
            
            factor = min(handles.data.emheader{sel_image}.Size([1 2])' ./ imsize);
            br = 0;
            if (handles.data.config.omitborder)
                br = handles.data.config.omitborderpixel;
                if (br < 1)
                    br = br * min(handles.data.emheader{sel_image}.Size([1 2]));
                end;
                br = round(br / factor);
                if (~(br >= 0))
                    br = 0;
                end;
            end;
            if (br > 0)
                %disp(['   Border: ' num2str(br) ' pixel']);
                mask = mask((br+1):(end-br), (br+1):(end-br));
                im = im((br+1):(end-br), (br+1):(end-br));
            end;
            if (~isempty(im))
                imsize = size(mask);


                if (handles.data.config.mindist_auto)
                    mindist = max(0,round(0.45*sqrt(sum(mask(:))/cma)));
                else
                    mindist = handles.data.config.mindist / factor;
                    if (handles.data.config.mindist > 0 && handles.data.config.mindist < 1)
                        mindist = mindist * min(handles.data.emheader{sel_image}.Size([1 2]));
                    end;    
                end;
                %disp(['   Mindist: ' num2str(mindist) ' pixel']);

                clear_mask = getCircleMask(mindist);
                clear_mask_size = size(clear_mask, 1);
                clear_mask_centre = fix(size(clear_mask,1)/2)+1;



                %newms = tom_mark_initMarker(im, handles.data.config.nmarkers, use_minDistance, useborder);

                [eimage, threshold, bx, by] = edge(im, 'sobel');

    %            cim = abs(bx) .* abs(by);
                cim = bx.^2 .* by.^2;
                cim_min = min(cim(:));
                cim = (cim-cim_min)/(max(cim(:))-cim_min);
                cim(~mask) = 0;


                [cim_sort, cim_sort_idx] = sort(cim(:), 1, 'descend');
                cim_sort_idx = uint32(cim_sort_idx);
                [cim_sort_idxy, cim_sort_idxx] = ind2sub(imsize, cim_sort_idx);
                cim_sort_idxx = uint32(cim_sort_idxx);
                cim_sort_idxy = uint32(cim_sort_idxy);

                cim_sort_idxmax = find(cim_sort==0, 1);
                cim_inv = zeros(imsize, 'uint32');

                for (i=1:cim_sort_idxmax)
                    cim_inv(cim_sort_idxy(i), cim_sort_idxx(i)) = i;
                end;

                ima = 0;
                ii = 1;
                while (ima < cma && ii < cim_sort_idxmax)
                    if (mod(ima,100) && ~isempty(waitbarh))
                        waitbar((msize(3))/waitbar_cnt_tot, waitbarh);
                        waitbar_abbort = get(waitbarh, 'Userdata');
                        if (waitbar_abbort)
                            break;
                        end;
                    end;


                    while (ii < cim_sort_idxmax)
                        if (cim_sort(ii) > 0)
                            break;
                        end;
                        ii = ii+1;
                    end;
                    if (ii < cim_sort_idxmax)
                        mpoint = [cim_sort_idxx(ii) cim_sort_idxy(ii)];
                        msize(3) = msize(3)+1;
                        ms([1 2], sel_image, msize(3)) = mpoint;

                        if (clear_mask_centre == 1)
                            cim_sort(ii) = 0;
                            cim_inv(mpoint(2), mpoint(1)) = 0;
                        else
                            range_maskx = max(1, mpoint(1)-clear_mask_centre+1) : min(imsize(2), mpoint(1)+clear_mask_centre-1);
                            range_masky = max(1, mpoint(2)-clear_mask_centre+1) : min(imsize(1), mpoint(2)+clear_mask_centre-1);
                            range_clear_maskx = max(1, clear_mask_centre-mpoint(1)+1) : min(clear_mask_size, clear_mask_centre+imsize(2)-mpoint(1));
                            range_clear_masky = max(1, clear_mask_centre-mpoint(2)+1) : min(clear_mask_size, clear_mask_centre+imsize(1)-mpoint(2));

                            for (iix = 1:length(range_maskx))
                                for (iiy = 1:length(range_masky))
                                    if (clear_mask(range_clear_masky(iiy), range_clear_maskx(iix)))
                                        idxii = cim_inv(range_masky(iiy), range_maskx(iix));
                                        if (idxii ~= 0)
                                            cim_sort(idxii) = false;
                                            cim_inv(range_masky(iiy), range_maskx(iix)) = 0;
                                        end;
                                    end;
                                end;
                            end;

                            %imtt = im;imtt(cim_inv==0) = min(im(:));figure(1);hold off;imagesc(imtt);axis off equal;colormap gray;hold on;plot(squeeze(ms(1,sel_image,1:msize(3))),squeeze(ms(2,sel_image,1:msize(3))),'g+');drawnow;
                        end;


                        ima = ima+1;
                    end;
                end;

                if (br > 0) 
                    ms([1 2],sel_image,msize3_old+1:msize(3)) = ms([1 2],sel_image,msize3_old+1:msize(3)) + br;
                end;
                ms([1 2],sel_image,msize3_old+1:msize(3)) = tom_mark_applyTransf('markerset_post', ms([1 2],sel_image,msize3_old+1:msize(3)), handles.data.config.method_contrast_filter, handles.data.config.method_contrast_binning);

                if (waitbar_abbort)
                    break;
                end;
            end;            
        end;
        ms = ms([1 2],1:msize(2),1:msize(3));
        ms(~(ms>0)) = nan;
        
    otherwise
        error(['Unknown method ' handles.data.config.method]);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_cancel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_finished(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.exitstatus = get(hObject, 'Tag');
guidata(hObject, handles);

waiting = strcmp(get(handles.figure_markInitMarkerGui, 'waitstatus'), 'waiting');

uiresume(handles.figure_markInitMarkerGui);

if (~waiting)
    handles.data.exitstatus = [handles.data.exitstatus '_nowait'];
    tom_markInitMarkerGui_OutputFcn(hObject, eventdata, handles);
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads the configuration from the structure config0
% and saves it into handles.data.config.
% If fields are wrong or unspecified, a defaultvalue is
% set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = loadConfig(handles, config0)

lfilenames = length(handles.data.filenames);
minn = inf();
minidx = [];
for (i=1:lfilenames)
    if (~isempty(handles.data.emheader{i}) && abs(handles.data.emheader{i}.Tiltangle) <= minn)
        minn = abs(handles.data.emheader{i}.Tiltangle);
        minidx = i;
    end;
end;
configdefault.sel_images = minidx;
configdefault.nmarkers = 50;
configdefault.omitborder = true;
configdefault.omitborderpixel = 0.05;
configdefault.usemask = true;
configdefault.mask = cell(1,lfilenames);
configdefault.maskwidth = 500;
configdefault.mindist_auto = true;
configdefault.mindist = 0;
configdefault.method = 'random';
configdefault.method_contrast_binning = 0;
configdefault.method_contrast_filter = [];


if (~isstruct(config0) || isempty(config0))
    config0 = struct();
end;

if (isfield(config0, 'sel_images'))
    sel_images = unique(config0.sel_images((config0.sel_images>0) & (config0.sel_images<=lfilenames)));
    handles.data.config.sel_images = [];
    for (i=1:length(sel_images))
        if (~isempty(handles.data.emheader{sel_images(i)}))
            handles.data.config.sel_images = [handles.data.config.sel_images sel_images(i)];
        end;
    end;
    if (isempty(handles.data.config.sel_images))
        handles.data.config.sel_images = configdefault.sel_images;
    end;
else
    handles.data.config.sel_images = configdefault.sel_images;
end;

if (isfield(config0, 'nmarkers'))
    handles.data.config.nmarkers = config0.nmarkers;
else
    handles.data.config.nmarkers = configdefault.nmarkers;
end;

if (isfield(config0, 'omitborder'))
    handles.data.config.omitborder = config0.omitborder;
else
    handles.data.config.omitborder = configdefault.omitborder;
end;

if (isfield(config0, 'omitborderpixel'))
    handles.data.config.omitborderpixel = config0.omitborderpixel;
else
    handles.data.config.omitborderpixel = configdefault.omitborderpixel;
end;

if (isfield(config0, 'usemask'))
    handles.data.config.usemask = config0.usemask;
else
    handles.data.config.usemask = configdefault.usemask;
end;
if (isfield(config0, 'mask') && iscell(config0.mask) && length(config0.mask)==lfilenames)
    handles.data.config.mask = config0.mask;
else
    handles.data.config.mask = configdefault.mask;
end;
if (isfield(config0, 'maskwidth'))
    handles.data.config.maskwidth = config0.maskwidth;
else
    handles.data.config.maskwidth = configdefault.maskwidth;
end;
if (isfield(config0, 'mindist_auto'))
    handles.data.config.mindist_auto = config0.mindist_auto;
else
    handles.data.config.mindist_auto = configdefault.mindist_auto;
end;
if (isfield(config0, 'mindist'))
    handles.data.config.mindist = config0.mindist;
else
    handles.data.config.mindist = configdefault.mindist;
end;
if (~isfield(config0, 'method'))
    handles.data.config.method = configdefault.method;
elseif (any(strcmp(config0.method, {'random', 'contrast', 'grid'})));
    handles.data.config.method = config0.method;
else
    handles.data.config.method = configdefault.method;
end;
if (isfield(config0, 'method_contrast_binning'))
    handles.data.config.method_contrast_binning = config0.method_contrast_binning;
else
    handles.data.config.method_contrast_binning = configdefault.method_contrast_binning;
end;
 if (isfield(config0, 'method_contrast_filter'))
    handles.data.config.method_contrast_filter = config0.method_contrast_filter;
else
    handles.data.config.method_contrast_filter = configdefault.method_contrast_filter;
end;   

handles.data.configdefault = configdefault;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = loadGui(handles)


set(handles.edit_nmarkers, 'String', num2str(handles.data.config.nmarkers));
set(handles.checkbox_omitborder, 'Value', handles.data.config.omitborder);
set(handles.edit_omitborderpixel, 'String', num2str(handles.data.config.omitborderpixel));
set(handles.checkbox_usemask, 'Value', handles.data.config.usemask);
set(handles.radiobutton_mindist_automatic, 'Value', handles.data.config.mindist_auto);
set(handles.radiobutton_mindist_manual, 'Value', ~handles.data.config.mindist_auto);
set(handles.edit_mindist, 'String', num2str(handles.data.config.mindist));

set(handles.radiobutton_method_random, 'Value', strcmp(handles.data.config.method,'random'));
set(handles.radiobutton_method_contrast, 'Value', strcmp(handles.data.config.method,'contrast'));
set(handles.radiobutton_method_grid, 'Value', strcmp(handles.data.config.method,'grid'));
set(handles.edit_method_contrast_binning, 'String', num2str(handles.data.config.method_contrast_binning));



handles = Callback_listbox_images(handles.checkbox_usemask, [], handles); 


value = [];
cnt = 0;
existing_filename = false(1,length(handles.data.filenames));
for i=1:length(existing_filename)
    existing_filename(i) = ~isempty(handles.data.emheader{i});
    if (existing_filename(i))
        cnt = cnt+1;
        if (any(i==handles.data.config.sel_images))
            value(end+1) = cnt;
        end;
    end;
end;
existing_filename_idx = find(existing_filename);

set(handles.listbox_images, 'Value', value, 'Userdata', existing_filename_idx(value));    





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = unloadGui(handles)


[handles.data.config.nmarkers, valError] = getNumberFromEditField(handles.edit_nmarkers, handles.data.configdefault.nmarkers, true, @checkIntervalInteger, {0,50000});
handles.data.config.omitborder = get(handles.checkbox_omitborder, 'Value');
[handles.data.config.omitborderpixel, valError] = getNumberFromEditField(handles.edit_omitborderpixel, handles.data.configdefault.omitborderpixel, true, @checkOmitBoarder, {1000});
handles.data.config.usemask = get(handles.checkbox_usemask, 'Value');
handles.data.config.mindist_auto = get(handles.radiobutton_mindist_automatic, 'Value');
[handles.data.config.mindist, valError] = getNumberFromEditField(handles.edit_mindist, handles.data.configdefault.mindist, true, @checkIntervalInteger, {0,1000});


if (get(handles.radiobutton_method_contrast,'Value'))
    handles.data.config.method = 'contrast';
elseif (get(handles.radiobutton_method_random,'Value'))
    handles.data.config.method = 'random';
elseif (get(handles.radiobutton_method_grid,'Value'))
    handles.data.config.method = 'grid';
end;
[handles.data.config.method_contrast_binning, valError] = getNumberFromEditField(handles.edit_method_contrast_binning, handles.data.configdefault.method_contrast_binning, true, @checkIntervalInteger, {0,8});

handles.data.config.sel_images = get(handles.listbox_images, 'Userdata');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = checkIntervalInteger(min, max, value)
res = value>=min && value <=max && round(value)==value;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = checkOmitBoarder(max, value)
res = value>0 && (value<1 || (value <=max && round(value)==value));        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val, valError] = getNumberFromEditField(hObject, defaultval, setDefaultIfError, hFunction, varsFunction)

s = get(hObject, 'String');
n = str2double(s);

if (~isnan(n))
    varsFunction = {varsFunction{:}, n};
    if (~isempty(hFunction) && ~hFunction(varsFunction{:}))
        n = nan;
    end;
end;
if (isnan(n))
    val = defaultval;
    valError = true;    
    if (setDefaultIfError)
        set(hObject, 'String', num2str(defaultval));
    end;
else
    val = n;
    valError = false;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_uicontextmenu_listbox_images(hObject, eventdata, handles) 

currentpoint = get(handles.figure_markInitMarkerGui, 'CurrentPoint');
listbox_position = get(handles.listbox_images, 'Position');
listbox_extent = get(handles.listbox_images, 'Extent');

value = fix((listbox_position(2)+listbox_position(4)-currentpoint(2))/listbox_extent(4)*length(get(handles.listbox_images, 'String')) + get(handles.listbox_images,'ListboxTop'));

set(handles.listbox_images, 'Value', value);
handles = Callback_listbox_images(handles.listbox_images, eventdata, handles);
sel_image = get(handles.listbox_images, 'Userdata');

if (get(handles.checkbox_usemask, 'Value'))
    if (isempty(handles.data.config.mask{sel_image}))
        set([handles.uimenu_image_showmask, handles.uimenu_image_clearmask, handles.uimenu_image_invertmask], 'Visible', 'off');
    else
        set([handles.uimenu_image_showmask, handles.uimenu_image_clearmask, handles.uimenu_image_invertmask], 'Visible', 'on');
    end;        
    set([handles.uimenu_image_addmask], 'Visible', 'on');
else
    set([handles.uimenu_image_showmask, handles.uimenu_image_invertmask, handles.uimenu_image_clearmask, handles.uimenu_image_addmask], 'Visible', 'off');
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_uimenu_image_openmask(hObject, eventdata, handles) 
sel_image = get(handles.listbox_images, 'UserData');
if (numel(sel_image)~=1 || sel_image<1 || sel_image>length(handles.data.filenames) || isempty(handles.data.emheader{sel_image}))
    return;
end;







try
    im = tom_emreadc(handles.data.filenames{sel_image}, 'binning', max(0,fix(log2(min(handles.data.emheader{sel_image}.Size([1 2])) / handles.data.config.maskwidth)))); 
catch
    return;
end;
imsize = size(im.Value);
im = imresize(im.Value', [handles.data.config.maskwidth * [imsize(2) imsize(1)] / max(imsize)], 'bilinear');


immax = max(im(:));
immin = min(im(:));

im = (im-immin) / (immax-immin);


oldmask = handles.data.config.mask{sel_image};
mask = imresize(logical(oldmask), size(im), 'bilinear');

h = findobj('Tag','figure_markInitMarkerGuiMASK');
if (isempty(h))
    h = figure();
end;
set(h, 'NumberTitle','off', 'WindowStyle', 'normal','Tag','figure_markInitMarkerGuiMASK', 'Name',['Show mask for image ' num2str(sel_image)], 'Menubar', 'figure');
figure(h);

im(~mask) = 0;

imagesc(im);
colormap('gray');
axis('equal', 'off');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_uimenu_image_invertmask(hObject, eventdata, handles) 


sel_image = get(handles.listbox_images, 'UserData');
if (numel(sel_image)~=1 || sel_image<1 || sel_image>length(handles.data.filenames) || isempty(handles.data.emheader{sel_image}))
    return;
end;

handles.data.config.mask{sel_image} = ~handles.data.config.mask{sel_image};
%handles = Callback_listbox_images(handles.uimenu_image_addmask, eventdata, handles); 
guidata(handles.figure_markInitMarkerGui, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_uimenu_image_addmask(hObject, eventdata, handles) 

sel_image = get(handles.listbox_images, 'UserData');
if (numel(sel_image)~=1 || sel_image<1 || sel_image>length(handles.data.filenames) || isempty(handles.data.emheader{sel_image}))
    return;
end;

try
    im = tom_emreadc(handles.data.filenames{sel_image}, 'binning', max(0,fix(log2(min(handles.data.emheader{sel_image}.Size([1 2])) / handles.data.config.maskwidth)))); 
catch
    return;
end;
imsize = size(im.Value);
im = imresize(im.Value', [handles.data.config.maskwidth * [imsize(2) imsize(1)] / max(imsize)], 'bilinear');


immax = max(im(:));
immin = min(im(:));

im = (im-immin) / (immax-immin);


oldmask = handles.data.config.mask{sel_image};
if (isempty(oldmask))
    oldmask = false(size(im));
else
    oldmask = imresize(logical(oldmask), size(im), 'bilinear');
end;



changed = false;
while (true)
    
    axes(findobj('Tag', 'axes_tom_markInitMarkerGui_dummy_roi'))
    
    im2 = im;
    im2(oldmask) = 1;
    

    h = findobj('Tag','figure_markInitMarkerGuiMASK');
    if (isempty(h))
        h = figure();
    else
        figure(h);
    end;
    set(h, 'NumberTitle','off', 'WindowStyle', 'modal','Tag','figure_markInitMarkerGuiMASK', 'Name',['Show mask for image ' num2str(sel_image)], 'Menubar', 'figure');

    newmask = [];
    try
        newmask = roipoly(im2);
    catch
    end;

    if (isempty(newmask))
        if (changed)
            newmask = oldmask;
        else
            break;
        end;
    else
        changed = true;
        newmask = oldmask | newmask;
    end;
    
    im2 = im;
    im2(~newmask) = 0;

    imagesc(im2);
    colormap('gray');
    axis('equal', 'off');


    answ = questdlg('Take the selected mask?', 'Add mask', 'Yes', 'Continue', 'Cancel', 'Yes');

    if (strcmp(answ, 'Continue'))
        oldmask = newmask;
    else
        if (strcmp(answ, 'Yes'))
            handles.data.config.mask{sel_image} = newmask;
            handles = Callback_listbox_images(handles.uimenu_image_addmask, eventdata, handles); 
            guidata(handles.figure_markInitMarkerGui, handles);
        end;
        break;
    end;
end;

if (ishandle(h))
    delete(h);
end;
return;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_uimenu_image_clearmask(hObject, eventdata, handles) 

sel_image = get(handles.listbox_images, 'UserData');
if (numel(sel_image)~=1 || sel_image<1 || sel_image>length(handles.data.filenames) || isempty(handles.data.emheader{sel_image}))
    return;
end;

handles.data.config.mask{sel_image} = [];
handles = Callback_listbox_images(handles.uimenu_image_addmask, eventdata, handles); 
guidata(handles.figure_markInitMarkerGui, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = Callback_listbox_images(hObject, eventdata, handles, value) 

value = get(handles.listbox_images, 'Value');


existing_filename = false(1,length(handles.data.filenames));
for i=1:length(existing_filename)
    existing_filename(i) = ~isempty(handles.data.emheader{i});
end;
existing_filename_idx = find(existing_filename);


if (hObject == handles.listbox_images)
    set(handles.listbox_images, 'Userdata', existing_filename_idx(value));    
else
    if (get(handles.checkbox_usemask,'Value'))
        for (i=1:length(handles.data.filenames))
            if (isempty(handles.data.config.mask{i}))
                s{i} = [handles.data.fileinfo(i).dispname ' (no mask defined)'];
            else
                s{i} = [handles.data.fileinfo(i).dispname ' (mask defined)'];
            end;
        end;
    else
        s = { handles.data.fileinfo(:).dispname };
    end;
    set(handles.listbox_images, 'String', s(existing_filename_idx), 'Userdata', existing_filename_idx(value));
end;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = getFileInfo_info(handles, idx)

if (~exist('idx','var'))
    idx = 1:length(handles.data.filenames);
end;

handles.data.fileinfo = struct('filename', {}, 'dirname', {}, 'dispname', {});

currn = 1;
for (i=1:length(idx))
    [pathstr, name, ext] = fileparts(handles.data.filenames{idx(i)});
    %[pathstr, name, ext, versn] = fileparts(handles.data.filenames{idx(i)});
    versn='';
    
    handles.data.fileinfo(i).filename = [name, ext, versn];
    handles.data.fileinfo(i).dirname = pathstr;
    if (isempty(handles.data.emheader{i}))
        handles.data.fileinfo(i).dispname = [sprintf('%3d :        : ', i) name ext versn];
    else
        handles.data.fileinfo(i).dispname = [sprintf('%3d : %8.4f° : ', i, handles.data.emheader{i}.Tiltangle) name ext versn];
    end;    

end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_checkbox_usemask(hObject, eventdata, handles)


Callback_listbox_images(handles.checkbox_usemask, [], handles);

















