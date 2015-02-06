%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_markGui(varargin)
% TOM_MARKGUI M-file for tom_mark_autoMarker.fig
%      TOM_MARK_AUTOMARKER, by itself, creates a new TOM_MARK_AUTOMARKER or raises the existing
%      singleton*.
%
%      H = TOM_MARK_AUTOMARKER returns the handle to a new TOM_MARK_AUTOMARKER or the handle to
%      the existing singleton*.
%
%      TOM_MARK_AUTOMARKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_MARK_AUTOMARKER.M with the given input arguments.
%
%      TOM_MARK_AUTOMARKER('Property','Value',...) creates a new TOM_MARK_AUTOMARKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_mark_autoMarker_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_markGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_mark_autoMarker

% Last Modified by GUIDE v2.5 07-Jun-2007 18:56:34

% Begin initialization code - DO NOT EDIT

%TIME0 = cputime();
%tic;

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_markGui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_markGui_OutputFcn, ...
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

%disp(['-- tom_markGui took ' num2str(cputime()-TIME0) ' seconds!']);
%toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before tom_mark_autoMarker is made visible.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tom_markGui_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_mark_autoMarker (see VARARGIN)


%set(handles.mainmenu_project_exporthandles, 'Visible', 'off');


accessImageInfoDispOptions(handles.pushbutton_image_info_dispoptions, 4);
set(handles.listbox_images, 'Units', get(handles.uipanel_listbox_images_position, 'Units'));

% tom_imagesc creates a uimenu if it does not exist.
% Create invisible dummy here...
if (isempty(findobj('label', 'Process')))
    uimenu (gcf,'label', 'Process', 'Visible', 'off');    	    
end;
% Dummy to use tom_load_tiltseries
if (isempty(findobj('Tag', 'menu'))) 
    uicontrol(gcf, 'TooltipString', 'Dummy for tom_fig_menu', 'Tag', 'menu', 'Visible', 'off');
end;

handles.fcn_local.getTiltangles = @getTiltangles;
handles.fcn_local.getListboxValue = @getListboxValue;

handles = setEmptyProject(handles);
sel_images = [];
sel_markers = [];


set(handles.figure_markGui, 'Menubar', 'none')
viewmenufcn(handles.figure_markGui, 'FigureToolbar');

hTB = findall(handles.figure_markGui, 'Tag', 'FigureToolBar', 'Type', 'uitoolbar');
hChilds = findall(handles.figure_markGui, 'Parent', hTB);
hActive = findall(hChilds, 'flat', 'Tag', 'Exploration.ZoomOut', '-or', 'Tag', 'Exploration.ZoomIn', '-or', 'Tag', 'Exploration.Pan', '-or', 'Tag', 'Exploration.DataCursor');
hDeActive = setdiff(hChilds, hActive);

set(hDeActive,'Visible', 'off', 'Enable', 'off');
set([hTB; hActive],'Visible', 'on');





updateListboxDisplayLockImages(handles.listbox_images, handles.uipanel_listbox_images_position, length(handles.data.markersets));
updateListboxDisplay(handles);
handles = updateSelection(handles, sel_images, sel_markers);


set(handles.mainmenu_markerset, 'Enable', getOnOff(~isempty(handles.data.filenames)));

set(handles.listbox_markersets, 'KeyPressFcn', @KeyPressFcn_listbox_markersets);
set(handles.figure_markGui, 'KeyPressFcn', @KeyPressFcn_figure_markGui);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sstatus = getOnOff(status)
if (status)
    sstatus = 'on';
else
    sstatus = 'off';
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = callLocalFcn(hObjectDummy, fname, varargin)
[varargout{1:nargout}] = feval(fname, varargin{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets the field "data" from the handles structure
% as there is an new and empty project. deletes nearly all 
% data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = setEmptyProject(handles)




set(handles.axes_main, 'Userdata', struct('uri', '', 'wrect', [0 0; 0 0], 'imreadbinning', 0, 'force', false));

handles.data.filenames = {};
handles.data.transf = [];
handles.data.sel_markerset = [];
handles.data.findmatchesparam = struct();
handles.data.emheader = getFileInfo_emheader(handles.data.filenames);
handles.data.fileinfo = getFileInfo_info(handles, 1:length(handles.data.filenames));
handles.data.markersets = struct('name', {},'filename', {}, 'markerset', {}, 'findmatchesparam', {});
handles.data.clipboard.markerset = nan(2, 0, 0);


if (~isfield(handles, 'workingdir') || ~exist(handles.workingdir, 'dir'))
    handles.workingdir = pwd();
end;


handles.data.initMarkerGuiData = [];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs from this function are returned to the command line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_markGui_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL,INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ms = getDefaultMarkerset(names, nImages)

lnames = length(names);

i = 0;
s = '';
while (i <= lnames && isempty(s))
    i = i+1;
    s = ['< ' num2str(i) ' >'];
    for (j=1:lnames)
        if (strcmp(s, names{j}))
            s = '';
            break
        end;
    end;
end;

ms = struct;
ms.name = s;
ms.filename = 'markerfile.em';
ms.markerset = nan(2,nImages,0);
ms.findmatchesparam = struct('comparelist', {zeros(0, 3)}, 'matchchain', {nan(1, nImages, 0)});







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Searches for the best matching existing filename.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filen = getExistingFilenameFromString(s, workingdir)

if (~exist('workingdir', 'var') || ~exist(workingdir,'dir'))
    workingdir = pwd();
end;

if (isempty(s))
    filen = fullfile(workingdir, 'a');
    return;
end;

dirprefix = '';
if (s(1) ~= filesep())
    if (ispc() && any(ismember(['a':'z','A':'Z'], s(1))) && length(s)>1 &&  s(2)==':')
        % is it an M$-absolute filename????
    else
        dirprefix = workingdir;
    end;
end;


[s, pname, pext] = fileparts(s);

pversn='';

while (~isempty(s))
    if (exist(fullfile(dirprefix, s), 'dir'))
        break;
    end;
    [s, pname, pext] = fileparts(s);
end;

filen = fullfile(dirprefix, s, [pname, pext, pversn]);








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ButtonDownFcn_info_marker_HighlightSelectedMarker(hObject, eventdata, handles) 


sel_marker = getListboxValue(handles.listbox_markers);
sel_image = getListboxValue(handles.listbox_images);

if (isempty(handles.data.markersets) || isempty(sel_image) || isempty(sel_marker))
    return;
end;

sel_marker = sel_marker(1);
sel_image = sel_image(1);

marker = handles.data.markersets(handles.data.sel_markerset).markerset(:,sel_image,sel_marker);




if (all(marker>0))
    t = timerfind('Name', 'HighlightMarkerTimer');
    if (length(t)~=1)
        delete(findobj('Type', 'line', 'Tag', 'HighlightMarker'));
        delete(t);
        
        h = [];
        t = timer('Name', 'HighlightMarkerTimer', 'Period', 0.01, 'ExecutionMode', 'fixedRate', 'Busymode', 'drop');
    else
        stop(t);
        h = get(t, 'Userdata');
    end;
    
    
    limits = [get(handles.axes_main, 'XLim'); get(handles.axes_main, 'YLim')];
    if (marker(1)>=limits(1,1) && marker(1)<=limits(1,2) && ...
        marker(2)>=limits(2,1) && marker(2)<=limits(2,2))
        if (ishandle(h))
            set(h, 'X', marker(1), 'Y', marker(2), 'Color', getColorFromRange(sel_marker, size(handles.data.markersets(handles.data.sel_markerset).markerset, 3)), ...
                'MarkerSize', 500, 'LineWidth', 6); 
        else
            hold(handles.axes_main, 'on'); 
            h = plot(handles.axes_main, marker(1), marker(2), 'o', 'Tag', 'HighlightMarker', 'Color', getColorFromRange(sel_marker, size(handles.data.markersets(handles.data.sel_markerset).markerset, 3)), ...
                     'MarkerSize', 500, 'LineWidth', 6); 
            hold(handles.axes_main, 'off'); 
        
            set(t, 'TimerFcn', { @ButtonDownFcn_text_infomarker_number_TimerFcn, h }, 'Userdata', h);
        end;

        start(t);
    end;
end;

if (get(handles.checkbox_infomarker_set_marker, 'Value'))
    userdata = get(handles.checkbox_infomarker_set_marker, 'Userdata');
    if (isempty(userdata))
        userdata.direction = 1;
    else
        userdata.direction = userdata.direction * -1;
    end;
    set(handles.checkbox_infomarker_set_marker, 'Userdata', userdata, 'BackgroundColor', getDirectionLabelColor(userdata.direction));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bgc = getDirectionLabelColor(direction)

if (direction == 1)
    bgc = 'green';
else
    bgc = [17, 196, 255]/255;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ButtonDownFcn_text_infomarker_number_TimerFcn(obj, event, handle)

finished = ~ishandle(handle);
if (~finished)
    msize = get(handle, 'MarkerSize');
    if (msize > 15)
        set(handle, 'MarkerSize', max(1, fix(msize*0.80)));
%        set(handle, 'Color', rand(1,3));
    else
        delete(handle);
        finished = true;
    end;
end;

if (finished)
    stop(obj);
end;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_markers_paste(hObject, eventdata, handles) 

if (isempty(handles.data.clipboard.markerset) || isempty(handles.data.markersets))
    return;
end;

sel_markers = getListboxValue(handles.listbox_markers);
sel_markerset = handles.data.sel_markerset;

msize = size3(handles.data.markersets(sel_markerset).markerset);
msize2 = size3(handles.data.clipboard.markerset);

if (msize2(2) ~= msize(2) || msize2(3) < 1 || msize2(1) ~= msize(1))
    return;
end;

if (msize(3) > 0)
    markerset = cat(3, handles.data.markersets(sel_markerset).markerset, handles.data.clipboard.markerset);
else
    markerset = handles.data.clipboard.markerset;
end;


% Paste/modify the matchchain
newmatchchain = nan([1, msize2(2), msize2(3)]);
newmatchchain(all(isfinite(handles.data.clipboard.markerset), 1)) = 0;
handles.data.markersets(sel_markerset).findmatchesparam.matchchain = cat(3, handles.data.markersets(sel_markerset).findmatchesparam.matchchain, newmatchchain);


handles.data.markersets(sel_markerset).markerset = markerset;

sel_markers = [(msize(3)+1):(msize(3)+msize2(3)), sel_markers];

handles.data.fileinfo = getFileInfo_info(handles, 1:length(handles.data.filenames));
sel_images = getListboxValue(handles.listbox_images);

updateListboxDisplay(handles);

updateSelection(handles, sel_images, sel_markers);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns a vector of length 3 with the size of the markerset.
% (used to substitute the buildin-function siye, which returns an vector of
% length 2 in case of ms = 2xNx1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msize = size3(ms)
msize = [0, 0, 0];
[msize(1), msize(2), msize(3)] = size(ms);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copies the selected markers into the "clipboard"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_markers_copy(hObject, eventdata, handles) 

sel_markers = getListboxValue(handles.listbox_markers);
if (isempty(sel_markers))
    return;
end

handles.data.clipboard.markerset = handles.data.markersets(handles.data.sel_markerset).markerset(:, :, sort(sel_markers));
guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Puts the selected markers into the clipboard and deletes
% them from the current markerset.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_markers_cut(hObject, eventdata, handles) 

idxcut = sort(getListboxValue(handles.listbox_markers));
if (isempty(idxcut))
    return;
end
sel_markerset = handles.data.sel_markerset;
msize = size3(handles.data.markersets(sel_markerset).markerset);

idxrem = sort(setdiff(1:msize(3), idxcut));

handles.data.clipboard.markerset = handles.data.markersets(sel_markerset).markerset(:,:,idxcut);

% Paste/modify the matchchain
handles.data.markersets(sel_markerset).markerset(:,:,idxcut) = [];
handles.data.markersets(sel_markerset).findmatchesparam.matchchain(:,:,idxcut) = [];
handles.data.markersets(sel_markerset).findmatchesparam.comparelist ...
    (ismember(handles.data.markersets(sel_markerset).findmatchesparam.comparelist(:,1), idxcut), :) = [];


handles.data.fileinfo = getFileInfo_info(handles, 1:length(handles.data.filenames));

sel_images = getListboxValue(handles.listbox_images);
updateListboxDisplay(handles);
updateSelection(handles, sel_images, []);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in listbox_images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = Callback_listbox_images(hObject, eventdata, handles) %#ok<DEFNU,INUSL>


selectionType = get(handles.figure_markGui, 'SelectionType');
sel_images = getListboxValue(handles.listbox_images);

switch (selectionType) 
    case {'open'}
        if (length(sel_images)==1)
            newname = inputdlg(['Enter the new filename for "' handles.data.filenames{sel_images(1)} '".'], 'Change filename', 1, handles.data.filenames(sel_images(1)), 'on'); 
            if (length(newname)==1 && ~strcmp(newname, handles.data.filenames{sel_images(1)}))
                newname = newname{1};
                change = true;
                if (exist(newname, 'file') ~= 2)
                    change = strcmp(questdlg('There exists not file under the given name. Change anyway?', 'Change filename'), 'Yes');
                end;
                if (change)
                    handles.data.filenames{sel_images} = newname;
                    handles.data.emheader(sel_images) = getFileInfo_emheader({newname});
                    handles.data.fileinfo(sel_images) = getFileInfo_info(handles, sel_images);
                    updateListboxDisplay(handles);
                end;
            end;
        end;
    otherwise
end;

updateSelection(handles, sel_images, getListboxValue(handles.listbox_markers));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in listbox_images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = Callback_pushbutton_info_image_open_embrowse(hObject, eventdata, handles) %#ok<DEFNU,INUSL>

sel_images = getListboxValue(handles.listbox_images);
if (isempty(sel_images))
    warndlg('No image selected', 'open tom_embrowse', 'modal');
    return;
end;

newpath = fileparts(handles.data.filenames{sel_images(1)});
if (~exist(newpath, 'dir'))
    warndlg(['Imagepath does not exist (' newpath ')'], 'open tom_embrowse', 'modal');
    return;
end;

oldcd = cd();
cd(newpath);
tom_embrowse;
cd(oldcd);











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_images_sort(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
warning('TODO: check');
errordlg('NOT YET IMPELEMTED', '', 'modal');
return;


lfilenames = length(handles.data.filenames);

tiltangles = nan(1,lfilenames);
for (i=1:lfilenames)
    if (~isempty(handles.data.cache(i).emheader))
        tiltangles(i) = handles.data.cache(i).emheader.Tiltangle;
    end;
end;

[tiltangles_sortet, tiltangles_index] = sort(tiltangles);

handles = change_imageset_to_subtupel_of_images(handles, tiltangles_index);

sel_images = getListboxValue(handles.listbox_images);
for (i=1:length(sel_images))
    sel_images(i) = find(tiltangles_index==sel_images(i));
end;

sel_markers = getListboxValue(handles.listbox_markers);
updateListboxDisplay(handles);

updateSelection(handles, sel_images, sel_markers);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns a cell array containing the em-headers read from
% filenames. It uses tom_reademheader.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function emheaders = getFileInfo_emheader(filenames)

lfilenames = length(filenames);
emheaders = cell(1, lfilenames);
for (i=1:lfilenames)
    tmpheader.Header = [];
    try
        tmpheader = tom_reademheader(filenames{i});
    catch
        tmpheader.Header = [];
    end;
    if (~isempty(tmpheader))
        emheaders{i} = tmpheader.Header;
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fileinfo = getFileInfo_info(handles, idx)

fileinfo = struct('filename', {}, 'dirname', {}, 'dispname', {});
for (i=1:length(idx))
    [pathstr, name, ext] = fileparts(handles.data.filenames{idx(i)});
    fileinfo(i).filename = [name, ext];
    fileinfo(i).dirname = pathstr;
    fileinfo(i).dispname = getImageName(handles, idx(i));
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = updateSelection(handles, sel_images, sel_markers)

if (exist('sel_markers','var'))
    sel_images = setListboxValue(handles.listbox_images, sel_images);
    sel_markers = setListboxValue(handles.listbox_markers, sel_markers);
else
    if (exist('sel_images', 'var'))
        sel_images = setListboxValue(handles.listbox_images, sel_images);
    else
        sel_images = getListboxValue(handles.listbox_images);

    end;
    sel_markers = getListboxValue(handles.listbox_markers);
end;

sel_images_sort = sort(sel_images);
sel_markers_sort = sort(sel_markers);


% Set the content of the listbox with the IMAGES WHERE THE SELECTED MARKERS
% OCCUR!
if (~isempty(sel_markers))
    existing_markerset_idx = all(handles.data.markersets(handles.data.sel_markerset).markerset(:,:,sel_markers_sort) > 0, 1);

    if (get(handles.checkbox_images_marked_invert, 'Value'))
        index = find(~any(existing_markerset_idx, 3));
    else
        index = find(any(existing_markerset_idx, 3));
    end;
    lindex = length(index);

    sel_images_subset = find(ismember(index, intersect(sel_images_sort, index)));
    s = { handles.data.fileinfo(index).dispname };
    listboxTop = max(1,min(get(handles.listbox_images_marked, 'ListBoxTop'), lindex));
else
    sel_images_subset = [];
    s = {};
    listboxTop = 1;
    index = [];
end;
set(handles.listbox_images_marked, 'String', s, 'Value', sel_images_subset, ...
    'ListBoxTop', listboxTop, 'Userdata', index);




% Set the content of the listbox with the MARKERS WHICH ARE IN THE SELECTED
% IMAGES!
if (isempty(handles.data.markersets) || isempty(handles.data.markersets(handles.data.sel_markerset).markerset) || isempty(sel_images_sort))
    s = {};
    sel_markers_subset = [];
    listboxTop = 1;
    index = [];
else
    existing_markerset_idx = all(handles.data.markersets(handles.data.sel_markerset).markerset(:,sel_images_sort,:) > 0, 1);
    if (get(handles.checkbox_images_markers_invert, 'Value'))
        index = find(~any(existing_markerset_idx, 2));
    else
        index = find(any(existing_markerset_idx, 2));
    end;
    lindex = length(index);
    s = get(handles.listbox_markers, 'String');
    s = s(index);
    sel_markers_subset = find(ismember(index, intersect(sel_markers_sort, index)));
    listboxTop = max(1,min(get(handles.listbox_images_markers, 'ListBoxTop'),lindex));
end;
set(handles.listbox_images_markers, 'String', s, 'Value', sel_markers_subset, ...
    'ListBoxTop', listboxTop, 'Userdata', index');
    
handles = updateImageDisplay(handles, sel_images, sel_markers);


guidata(handles.figure_markGui, handles);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_images_down.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_images_down(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
warning('TODO: check');
errordlg('NOT YET IMPELEMTED', '', 'modal');
return;

sel_images = getListboxValue(handles.listbox_images);
if (isempty(sel_images))
    return;
end;
sel_images_sort = sort(sel_images, 'descend');
lsel_images_sort = length(sel_images_sort);

lfilenames = length(handles.data.filenames);
index = 1:lfilenames;


for (i=1:lsel_images_sort)
    if (sel_images_sort(i) < lfilenames-i+1)
        index([sel_images_sort(i)+1,sel_images_sort(i)]) = index([sel_images_sort(i),sel_images_sort(i)+1]);
    end;
end;

handles = change_imageset_to_subtupel_of_images(handles, index);

for (i=1:length(sel_images))
    sel_images(i) = find(index==sel_images(i));
end;

sel_markers = getListboxValue(handles.listbox_markers);

updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_images_up.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_images_up(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
warning('TODO: check');
errordlg('NOT YET IMPELEMTED', '', 'modal');
return;

sel_images = getListboxValue(handles.listbox_images);
if (isempty(sel_images))
    return;
end;
sel_images_sort = sort(sel_images, 'ascend');
lsel_images_sort = length(sel_images_sort);

lfilenames = length(handles.data.filenames);
index = 1:lfilenames;


for (i=1:lsel_images_sort)
    if (i<sel_images_sort(i))
        index([sel_images_sort(i)-1,sel_images_sort(i)]) = index([sel_images_sort(i),sel_images_sort(i)-1]);
    end;
end;


handles = change_imageset_to_subtupel_of_images(handles, index);

for (i=1:length(sel_images))
    sel_images(i) = find(index==sel_images(i));
end;

sel_markers = getListboxValue(handles.listbox_markers);

updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quits tom_markGui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_project_quit(hObject, eventdata, handles) 
answ = questdlg('Do you really want to quit?', 'Quit tom_markGui','Yes','No','Yes');
if (~strcmp(answ, 'Yes'))
    return;
end 
close(handles.figure_markGui);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_images_delete.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_project_importhandles(hObject, eventdata, handles) 

answ = questdlg('Do you really want to replace the application data with the content from a variable from the workspace? The old data will be lost and in case of inconsistency tom_markGui will probably crash :)', 'Import handles', 'Continue', 'No', 'No');
if (~strcmp(answ, 'Continue'))
    return;
end;

vars = evalin('base', 'whos');

s = {};
for (i=1:length(vars))
    if (strcmp(vars(i).class, 'struct'))
        fields = evalin('base', ['fieldnames(' vars(i).name ');']);
        if (sum(strcmp(fields, 'data')) == 1)
            s{end+1} = vars(i).name;
        end;
    end;
end;

if (isempty(s))
    errordlg('In the base workspace there is no usable struct-variable with a datafield "data".', 'Import handles', 'modal');
    return;
end;

[selection, answ] = listdlg('PromptString', 'Select the workspace variable', 'SelectionMode', 'single', 'ListString' , s, 'InitialValue', 1, 'Name', 'Import handles', 'ListSize', [300, 200]);
if (~answ)
    return;
end;

handles.data = evalin('base', [s{selection} '.data;']);

handles.data.emheader = getFileInfo_emheader(handles.data.filenames);
if (isempty(handles.data.markersets))
    handles.data.sel_markerset = [];
elseif (handles.data.sel_markerset < 1 || handles.data.sel_markerset>length(handles.data.markersets))
    handles.data.sel_markerset = 1;
end;

tiltangles = Inf(1, length(handles.data.filenames));
for (i=1:length(tiltangles))
    if (~isempty(handles.data.emheader{i}))
        tiltangles(i) = handles.data.emheader{i}.Tiltangle;
    end;
end;
[minangle, sel_images] = min(abs(tiltangles));


if (isempty(handles.data.markersets))
    sel_markers = [];
else
    sel_markers = 1:size(handles.data.markersets(handles.data.sel_markerset).markerset, 3);
end;


cla(handles.axes_main);
userdata = get(handles.axes_main, 'Userdata');
userdata.force = true;
set(handles.axes_main, 'Userdata', userdata);

set(handles.mainmenu_markerset, 'Enable', getOnOff(~isempty(handles.data.filenames)));

updateListboxDisplayLockImages(handles.listbox_images, handles.uipanel_listbox_images_position, length(handles.data.markersets));
updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_images_add.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_images_add(hObject, eventdata, handles) %#ok<DEFNU,INUSL>
warning('TODO: check');
errordlg('NOT YET IMPELEMTED', '', 'modal');
return;
% There seams to be a bug in uigetfile with MultiSelect 'on'.
% If the user sellects a file with the mouse, nothing is inserted to the
% editfield for the filename. When the user then presses OK, uigetfile
% complains, that no file is selected.
% As current workaround use tom_mark_uiGetFiles!
%[filename, pathname] = uigetfile({'*.em', 'EM-Image'; '*.*',  'All Files (*.*)'}, 'Add files', 'MultiSelect', 'on');
[filename, pathname] = tom_mark_uiGetFiles();


if (isequal(filename,0) || isequal(pathname,0))
    return;
end;
if (ischar(filename))
    filename = {filename};
end;

lfilename = length(filename);

lfilenames = length(handles.data.filenames);
index = (lfilenames+1) : (lfilenames+lfilename);

for (i=1:lfilename)
    handles.data.filenames{end+1} = fullfile(pathname, filename{i});
end;
handles.data.markerset(:,index,:) = nan;

sel_markers = getListboxValue(handles.listbox_markers);

set(handles.mainmenu_markerset, 'Enable', getOnOff(~isempty(handles.data.filenames)));

updateListboxDisplay(handles);
updateSelection(handles, index, sel_markers);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_images_delete.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_images_delete(hObject, eventdata, handles) %#ok<DEFNU,INUSL>
warning('TODO: check');
errordlg('NOT YET IMPELEMTED', '', 'modal');
return;
lfilenames = length(handles.data.filenames);
value = getListboxValue(handles.listbox_images);
lvalue = length(value);

if (lvalue >= 1)
    if (strcmp(questdlg(['Remove the ' num2str(lvalue) ' filenames from the list?'], 'Delete'), 'Yes')) 
        index = setdiff(1:lfilenames, value);
        handles = change_imageset_to_subtupel_of_images(handles, index);
    else
        return;
    end;
elseif (lfilenames > 0)
    if (strcmp(questdlg(['Remove all the ' num2str(length(handles.data.filenames)) ' filenames from the list?'], 'Delete'), 'Yes'))
        handles = deldata_all_images(handles);
    else
        return;
    end;
else
    return;
end;

sel_markers = getListboxValue(handles.listbox_markers);

updateListboxDisplay(handles);
updateSelection(handles, [], sel_markers);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_images_load.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_images_load(hObject, eventdata, handles) %#ok<DEFNU,INUSL>

uiwait(tom_load_tiltseries);

m = findobj('Tag','menu');

if (length(m) ~= 1 || ~strcmp(get(m, 'TooltipString'), 'Dummy for tom_fig_menu'))
    disp('TODO: REPLACE tom_load_tiltseries by something better :) ...');
    error('Error in Callback_pushbutton_files_load. Maybe conflicting guiobject ''menu'' around. ABORT');
    return;
end;
param = get(m,'Userdata');
set(m, 'Userdata', 0);

if (strcmp(param.newproj_cancel, 'no'))
    
    param.myfirstnb = str2double(param.myfirstnb);
    param.mylastnb = str2double(param.mylastnb);

    filenames = {};
    
    if (isnan(param.myfirstnb) || param.myfirstnb < 1)
        warning('Callback_pushbutton_files_load: myfirstnb has a wrong value (must be a positive number)'); %#ok<WNTAG>
        return;
    elseif (isnan(param.mylastnb) || param.mylastnb < param.myfirstnb)
        warning('Callback_pushbutton_files_load: mylastnb has a wrong value (must be a positive number)'); %#ok<WNTAG>
        return;
    else
        for (i=param.myfirstnb:param.mylastnb)
            f = fullfile(param.mypathname, [param.myfilename num2str(i) param.myext]);
            if (exist(f, 'file') ~= 2)
                warning(['Selected file "' f '" were not found']); %#ok<WNTAG>
            end;
            filenames{end+1} = f; %#ok<AGROW>
        end;
%        warning('TODO: the markerfilename is not regarded!');
%        if (isempty(get(handles.edit_markerfile_filename, 'String')) || ...
%            strcmp(questdlg('Replace markerfile-entry in gui-element with new one?', 'Replace'), 'Yes'))
%            set(handles.edit_markerfile_filename, 'String', param.myfilemarker_default);
%        else
%        end;
    end;
    
    if (~isempty(handles.data.filenames)) 
        answ = questdlg(['There are already ' num2str(length(handles.data.filenames)) ' images loaded. Replace them, or append new one?'], 'Loading files', 'Replace', 'Append', 'Cancel', 'Append');
        if (strcmp(answ, 'Cancel'))
            return;
        elseif (strcmp(answ, 'Replace'))
            handles.data.filenames = {};
            handles.data.emheader = {};
            handles.data.fileinfo = struct('filename', {}, 'dirname', {}, 'dispname', {});
        end;
    end;
    new_index = length(handles.data.filenames)+(1:length(filenames));

    handles.data.filenames(new_index) = filenames;
    
    handles.data.emheader(new_index) = getFileInfo_emheader(filenames);
    handles.data.fileinfo(new_index) = getFileInfo_info(handles, new_index);
    
    cla(handles.axes_main);

    set(handles.mainmenu_markerset, 'Enable', getOnOff(~isempty(handles.data.filenames)));

    updateListboxDisplay(handles);
    updateSelection(handles, new_index, []);

end;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_project_load(hObject, eventdata, handles) %#ok<INUSL,DEFNU>


filename = get(handles.mainmenu_project_load, 'Userdata');
if (isempty(filename))
    filename = fullfile(handles.workingdir, 'markGui_project.mat');
end;
if (~exist(filename, 'file'))
    filename = getExistingFilenameFromString(filename);
end;


[filename, pathname, filterindex] = uigetfile({'*.mat','MAT-files (*.mat)'; ...
                                  '*.*',  'All Files (*.*)'}, ...
                                 'Select projectfile', filename);

if (isequal(filename, 0) || isequal(pathname, 0))
    return;
end;
fullfilename = fullfile(pathname, filename);


try
    data = load(fullfilename);
catch
    errordlg(['Error reading the project data from the mat-file "' filename '": ' err.message], 'Loading project', 'modal');
    return;
end;

% %, 'markGui_SelectMarker'
% datafields = {'filenames', 'transf', 'markersets', 'initMarkerGuiData', 'findmatchesparam','markGui_SelectMarker'};
% if (~all(isfield(data, datafields)) || length(fieldnames(data)) ~= length(datafields))
%     warning('Loading datafield');
%     errordlg(['Error loading the projectconfiguration from "' filename '". Wrong datafields!'], 'Loading project', 'modal');
%     return;
% end;

if (length(handles.data.filenames) == length(data.filenames) && ~isempty(data.filenames) && ~isempty(data.markersets) && any(~strcmp(handles.data.filenames, data.filenames)))
    answ = questdlg('The project you are about to load contains as many images as the currently loaded. Do you want to use the current images instead using the filenames saved in the new project?', 'Loading Project', 'Yes', 'No', 'Cancel', 'No');
    if (strcmp(answ, 'Cancel'))
        return;
    elseif (strcmp(answ, 'Yes'))
        data.filenames = handles.data.filenames;
    end;
end;

handles.data = data;
handles.data.emheader = getFileInfo_emheader(handles.data.filenames);
if (isempty(handles.data.markersets))
    handles.data.sel_markerset = [];
else
    handles.data.sel_markerset = 1;
end;
handles.data.fileinfo = getFileInfo_info(handles, 1:length(handles.data.filenames));
handles.data.clipboard.markerset = nan(2, 0, 0);





handles.workingdir = pathname;
set(handles.mainmenu_project_load, 'Userdata', fullfilename);

if (isempty(handles.data.markersets))
    sel_markers = [];
else
    sel_markers = 1:size(handles.data.markersets(handles.data.sel_markerset).markerset, 3);
end;

tiltangles = getTiltangles(handles.data.emheader);
[minangle, sel_images] = min(abs(tiltangles));


cla(handles.axes_main);
userdata = get(handles.axes_main, 'Userdata');
userdata.force = true;
set(handles.axes_main, 'Userdata', userdata);

set(handles.mainmenu_markerset, 'Enable', getOnOff(~isempty(handles.data.filenames)));

updateListboxDisplayLockImages(handles.listbox_images, handles.uipanel_listbox_images_position, length(handles.data.markersets));
updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [code] = accessImageInfoDispOptions(hObject, offset)

if (~exist('offset', 'var') || isempty(offset))
    code = get(hObject, 'Userdata');
else
    states = { 'binning'; 'sample'; 'binning (scaled)'; 'sample (scaled)' };
    userdata = {['b' '-']; ['s' '-']; ['b' 's']; ['s' 's']};

    lstates = length(states);
    s = get(hObject, 'String');

    switch (offset)
        case -1
            istate = mod(find(strcmp(states, s)), lstates) + 1;
        case -2
            istate = mod(find(strcmp(states, s))-2, lstates) + 1;
        otherwise 
            if (numel(offset)==1 && uint8(offset)==offset && offset >= 1 && offset <= lstates)
                istate = offset;
            else
                istate = 1;
            end;
    end;
    code = userdata{istate};
    set(hObject, 'String', states{istate}, 'userdata', code);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ButtonDownFcn_pushbutton_image_info_dispoptions(hObject, eventdata, handles) %#ok<INUSL,INUSL,DEFNU>

if (strcmp(get(handles.figure_markGui, 'SelectionType'), 'alt'))
    accessImageInfoDispOptions(hObject, -2);
else
    accessImageInfoDispOptions(hObject, -1);
end;

userdata = get(handles.axes_main, 'Userdata');
userdata.force = true;
set(handles.axes_main, 'Userdata', userdata);

updateImageDisplay(handles, ...
    getListboxValue(handles.listbox_images), ...
    getListboxValue(handles.listbox_markers));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_project_save(hObject, eventdata, handles) %#ok<INUSL,INUSL,DEFNU>


filename = get(handles.mainmenu_project_load, 'Userdata');
if (isempty(filename))
    filename = fullfile(handles.workingdir, 'markGui_project.mat');
end;
if (~exist(filename, 'file'))
    filename = getExistingFilenameFromString(filename);
end;

[filename, pathname, filterindex] = uiputfile({'*.mat','MAT-files (*.mat)'; ...
                                  '*.*',  'All Files (*.*)'}, ...
                                 'Select projectfile', filename);
if (isequal(filename,0) || isequal(pathname,0)) 
    return;
end;
fullfilename = fullfile(pathname, filename);

data = handles.data;
data = rmfield(data, {'clipboard','fileinfo','emheader', 'sel_markerset'});


try
    save(fullfilename, '-struct', 'data');
catch
    err = lasterror();
    errordlg(['Error saving file "' filename '": ' err.message], 'Saving project', 'modal');
    return;
end;

handles.workingdir = pathname;
set(handles.mainmenu_project_load, 'Userdata', fullfilename);

guidata(hObject, handles);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = updateImageDisplay(handles, sel_images, sel_markers, currentaxes)


if (isempty(handles.data.markersets))
    msize = [2, length(handles.data.filenames), 0];
    markerset = nan(msize);
else
    markerset = handles.data.markersets(handles.data.sel_markerset).markerset;
    msize = size(markerset);
    if (length(msize) == 2)
        msize(3) = 1;
    end;
end;

if (isempty(sel_images))
    sel_image = [];
else
    sel_image = sel_images(1);
end;
if (isempty(sel_markers))
    sel_marker = [];
else
    sel_marker = sel_markers(1);
end;

if (~exist('currentaxes', 'var'))
    currentaxes = handles.axes_main;
end;
axes(currentaxes);

children = findobj(get(currentaxes, 'Children'), 'flat', 'Type', 'line', '-or', 'Type', 'text');
delete(children);
himage = findobj(get(currentaxes, 'Children'), 'flat', 'Type', 'image');
if (isempty(sel_image))
    set(currentaxes, 'Visible', 'off');
    set(himage, 'Visible', 'off');
else
    if (isempty(handles.data.emheader{sel_image}))
        set(currentaxes, 'Visible', 'off');
        set(himage, 'Visible', 'off');
    else
        userdata = get(currentaxes, 'Userdata');
        wrect = ([get(currentaxes,'XLim');get(currentaxes,'YLim')]);
        
        if (userdata.force || ~strcmp(userdata.uri, handles.data.filenames{sel_image}) || any(any(wrect ~= userdata.wrect)))

            imsize = handles.data.emheader{sel_image}.Size([1 2]);

            gca_width = get(currentaxes, 'Position');
            gca_width = gca_width([3 4]);

            if (isempty(himage))
                rect = [[1, 1]', imsize];
                wrect = rect;
            else
                if (get(handles.checkbox_image_info_clamp, 'Value') && ~isempty(sel_marker))
                    if (currentaxes == handles.axes_main)
                        holdclamp = handles.checkbox_image_info_clamp;
                    else
                        holdclamp = findall(get(currentaxes, 'Parent'),'Tag', 'clamp_dummy');
                    end;
                    markerold = get(holdclamp, 'Userdata');
                    marker = handles.data.markersets(handles.data.sel_markerset).markerset(:,sel_image,sel_marker);
                    if (~all(marker>0))
                        set(holdclamp, 'Userdata', []);
                    else
                        set(holdclamp, 'Userdata', struct('sel_marker', sel_marker, 'marker' , marker));
                        if (~isempty(markerold) && markerold.sel_marker==sel_marker)
                            mdiff = marker - markerold.marker;
                            wrect(1,:) = wrect(1,:) + mdiff(1);
                            wrect(2,:) = wrect(2,:) + mdiff(2);
                        end;
                    end;
                    rect = [ max(1, floor(wrect(:,1))), min(imsize, ceil(wrect(:,2))) ];
                else
                    rect = [ max(1, floor(wrect(:,1))), min(imsize, ceil(wrect(:,2))) ];
                end;
            end;

            dispoptions_userdata = accessImageInfoDispOptions(handles.pushbutton_image_info_dispoptions);
            
%             imreadbinning = max(0, fix(log2(min((wrect(:,2) - wrect(:,1)) ./ gca_width'))));
%             if (dispoptions_userdata(1) == 's')
%                 rect = [[1, 1]', imsize];
%                 %im = tom_emreadc(handles.data.filenames{sel_image}, 'subregion', [rect(:,1)', 1], [rect(:,2)'-rect(:,1)', 0], 'resample', [1 1 1]*2^imreadbinning);
%                 im = tom_emreadc(handles.data.filenames{sel_image}, 'resample', [1 1 1]*2^imreadbinning);
%             else
%                 im = tom_emreadc(handles.data.filenames{sel_image}, 'subregion', [rect(:,1)', 1], [rect(:,2)'-rect(:,1)', 0], 'binning', imreadbinning);    
%             end;

            imreadbinning = [max([1;1], fix((wrect(:,2) - wrect(:,1)) ./ gca_width')); 1];

            if (dispoptions_userdata(1) == 's')
                im = tom_emreadc3(handles.data.filenames{sel_image}, [([rect(:,1)', 1]-1), ([rect(:,2)'-rect(:,1)', 0]+1)], imreadbinning, []);
            else
                im = tom_emreadc3(handles.data.filenames{sel_image}, [([rect(:,1)', 1]-1), ([rect(:,2)'-rect(:,1)', 0]+1)], [], imreadbinning);
            end;

            imValue = single(im.Value');
            
            if (dispoptions_userdata(2) == 's')
                v_min = min(imValue(:));
                v_max = max(imValue(:));
                v_meanv = mean(imValue(:));
                v_std = std2(imValue);
                range = [v_meanv-(5*v_std), v_meanv+(5*v_std)];
                himage = imagesc(rect(1,:), rect(2,:), imValue, range);
            else
                himage = imagesc(rect(1,:), rect(2,:), imValue);
            end;
            
            axis equal;
            colormap gray;

            xlim([wrect(1,1), wrect(1,2)]);
            ylim([wrect(2,1), wrect(2,2)]);

            userdata.uri = handles.data.filenames{sel_image};
            userdata.wrect = wrect;
            userdata.imreadbinning = imreadbinning;
            userdata.force = false;
            set(currentaxes, 'Userdata', userdata);
        else
            set(himage, 'Visible', 'on');            
        end;
        
        % Plot the markers.
        if (~isempty(sel_markers))
            hold('on');
            markers = markerset(:,sel_image,sel_markers);
            active_markers = find((markers(1,1,:)>=wrect(1,1)) & (markers(2,1,:)>=wrect(2,1)) & (markers(1,1,:)<=wrect(1,2)) & (markers(2,1,:)<=wrect(2,2)))';
            lactive_markers = length(active_markers);
            if (lactive_markers <= 100)
                for (i=active_markers)
                    marker = markers(1:2,1,i);
                    if (all(marker>=0,1))
                        color = getColorFromRange(sel_markers(i), msize(3));
                        if (true)
                            if (i==1)
                                plot(marker(1), marker(2), 'd', 'Color', color, 'MarkerSize', 15, 'LineWidth', 3);
                            else
                                plot(marker(1), marker(2), 'x', 'Color', color, 'MarkerSize', 12, 'LineWidth', 3);
                            end;
                            text(marker(1), marker(2), ['\leftarrow ' num2str(sel_markers(i))], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontName', 'FixedWidth');
                        end;
                    end;
                end;
            else
                color = getColorFromRange(sel_markers(1), msize(3));
                marker_highlight = active_markers(1) == 1;
                if (marker_highlight)
                    active_markers(1) = [];
                end;
                plot(squeeze(markers(1,1,active_markers)), squeeze(markers(2,1,active_markers)), 'x', 'Color', [1 1 1]-color, 'MarkerSize', 12, 'LineWidth', 3);
                if (marker_highlight)
                    marker = markers(1:2,1,1);
                    plot(marker(1), marker(2), 'd', 'Color', color, 'MarkerSize', 15, 'LineWidth', 3);
                    text(marker(1), marker(2), ['\leftarrow ' num2str(sel_markers(1))], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontName', 'FixedWidth');
                end;
            end;
            hold('off');
        end;
        
        set(currentaxes, 'Visible', 'on');
    end;    
end;




hHIDE = [ ...
        %handles.text_info_image_number, ...
        %handles.text_info_image_number_label, ...
        %handles.text_info_image_name, ...
        %handles.text_info_image_name_label, ...
        handles.text_info_image_size, ...
        handles.text_info_image_size_label, ...
        handles.text_info_image_angle, ...
        handles.text_info_image_angle_label, ...
        ];
    
if (isempty(sel_image))
    set(hHIDE , 'Visible', 'off');
    set(handles.text_info_image_number, 'String', 'none');
    set([handles.text_info_image_name_label, handles.text_info_image_name], 'Visible', 'off');
    uri = '';
else
    filename = handles.data.filenames{sel_image};
    uri = [sprintf('%3i', sel_image) '_' filename];
    olduri = get(handles.uipanel_image_info, 'Userdata');
    if (~strcmp(uri, olduri))
        emheader = handles.data.emheader{sel_image};

        set(handles.text_info_image_number, 'String', num2str(sel_image));
        set([handles.text_info_image_name, handles.text_info_image_name_label], 'Visible', 'on');
        set(handles.text_info_image_name, 'String', handles.data.fileinfo(sel_image).filename);
        
        if (isempty(emheader)) 
            set(hHIDE , 'Visible', 'off');
            set(handles.text_info_image_name, 'Background', 'red');
        else
            set(hHIDE , 'Visible', 'on');
            set(handles.text_info_image_name, 'Background', get(handles.uipanel_image_info, 'Background'));

            DataType = cellstr(['byte    ';'short   ';'        ';'long int';'float   ';'        ';'        ';'complex '; 'double  ']);
            set(handles.text_info_image_size, 'String', [num2str(emheader.Size(1)) ' x ' num2str(emheader.Size(2)) ' x ' num2str(emheader.Size(3))]);
            set(handles.text_info_image_angle, 'String', emheader.Parameter(19)./1000.0);
            
            %set(handles.text_info_image_type, 'String', DataType{emheader.Magic(4)});
            %set(handles.text_info_image_defocus, 'String', emheader.Defocus);
            %set(handles.text_info_image_tiltaxis, 'String', emheader.Tiltaxis);
            %set(handles.text_info_image_objpixsize, 'String', emheader.Objectpixelsize);
            %set(handles.text_info_image_comment, 'String', reshape(emheader.Comment, [1,numel(emheader.Comment)]));

            %fileinfo = dir(filename);
            %set(handles.text_info_image_lastmod, 'String', fileinfo.date); 
        end;
    end;
end;
set(handles.uipanel_image_info, 'Userdata', uri);

if (isempty(sel_marker))
    set(handles.text_infomarker_position, 'Visible', 'off');
    set(handles.text_infomarker_position_label, 'Visible', 'off');
    set(handles.text_infomarker_color, 'Background', get(handles.uipanel_infomarker_panel, 'Background'));
    set(handles.text_infomarker_number, 'String', 'none');
else
    marker = markerset(:, sel_image, sel_marker);
    if (any(marker>0))
        sposition = ['(' num2str(marker(1)) ',' num2str(marker(2)) ')'];
    else
        sposition = 'no point';
    end;
    set(handles.text_infomarker_position, 'Visible', 'on', 'String', sposition);
    set(handles.text_infomarker_position_label, 'Visible', 'on');
    set(handles.text_infomarker_color, 'Background', getColorFromRange(sel_marker, msize(3)));
    set(handles.text_infomarker_number, 'String', getMarkerName(handles, handles.data.sel_markerset, sel_marker));
end;

if (currentaxes == handles.axes_main)
    set(get(currentaxes, 'Children'), 'ButtonDownFcn', 'tom_markGui(''Callback_axes_main'',gcbo,[],guidata(gcbo))'); 
end;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For multiselection-listboxes save the order in which the
% items are selected in userdata.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newval = setListboxValue(hObject, value)

if (isempty(value))
    newval = [];
else
    oldval = get(hObject, 'Userdata');
    idx = ismember(value, oldval);
    newval = [value(~idx), value(idx)];
end;

set(hObject, 'Userdata', newval, 'Value', value);
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For multiselection-listboxes return the indexes of the
% selected items in the order as they were selected.
% Uses userdata!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value = getListboxValue(hObject)

value = get(hObject, 'Value');
oldval = get(hObject, 'Userdata');


oldval = oldval(ismember(oldval, value));

value = [setdiff(value, oldval), oldval];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function color = getColorFromRange(index, clength)

if (clength < 1 || index>clength)
    color = [1 0 0];
else
    %h = (0:clength-1)'/clength;
    %color = hsv2rgb([h ones(clength,2)]);
    color = hsv2rgb([(index-1)/clength 1 1]);
end;    
%cmap = hsv(clength);
%color = cmap(index, :);











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_project_new(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

answ = questdlg('Start new project? Unsaved data will be lost!', 'New project', 'Yes', 'No', 'Cancel', 'No');
if (~strcmp(answ, 'Yes'))
    return;
end;

handles = setEmptyProject(handles);

set(handles.mainmenu_markerset, 'Enable', getOnOff(~isempty(handles.data.filenames)));

updateListboxDisplay(handles);
updateSelection(handles, [], []);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_exportworkspace(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

if (isempty(handles.data.markersets))
    warndlg('No markerset available.', 'Save em-file', 'modal');
    return;
end;

ms = handles.data.markersets(handles.data.sel_markerset);

msize = size(ms.markerset);
if (length(msize)==2)
    msize(3) = 1;
end;
if (msize(3) == 0)
    warndlg('The selected markerset contains no markers.', 'Save em-file', 'modal');
    return;
end;

s = {['as 2x' num2str(msize(2)) 'x' num2str(msize(3)) ' markerset with unknown coordinates NaN'], ...
     ['as 3x' num2str(msize(2)) 'x' num2str(msize(3)) ' markerset with unknown coordinates -1 and saved Tiltangles'], ...
     ['in EM-Format']};
 
[selection, ok] = listdlg('ListString', s', 'SelectionMode', 'single', 'ListSize', [600,150], 'InitialValue', 1, 'Name', 'tom_markGui: Export Markerset', 'PromptString','choose Format:');

if (~ok)
    return;
end;

basename = regexprep(ms.name, '[^\w]', '_');
if (length(basename) < 1 || basename(1)=='_')
    basename = ['ms' basename];
end;


idx = 1;
name = basename;
while (true)
    if (~evalin('base', ['exist(''' name ''',''var'')']))
        break;
    end;    
    name = [basename '_' num2str(idx)];
    idx = idx+1;
end;


switch (selection)    
    case {1,2,3}
        ms_ = ms.markerset;
        if (selection>1)
            [Tiltangle, imsize] = getTiltangles(handles.data.emheader);

            invalididx = ~all(round(ms_) >= 1 & round(ms_) <= repmat(imsize, [1,1,msize(3)]), 1);
            j = sum(all(isfinite(ms_(1:2,invalididx)), 1));
            if (j > 0)
                warndlg(['There were ' num2str(j) ' markerpositions outside the image dimension. They were removed from the result (i.e. to -1). Use export as NaN to keep them.'], 'Export markerset', 'modal');
            end;
            ms_(1:2, invalididx) = -1;
            ms_2 = nan(3, msize(2), msize(3));
            ms_2(1,1:msize(2),1) = Tiltangle;
            ms_2(2:3,:,:) = ms_;

            if (selection == 3)
                ms_ = tom_emheader(ms_2);
            else
                ms_ = ms_2;
            end;
        end;
        
        assignin('base', name, ms_);
        
    otherwise
        warning('unexpected!!')
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_export_tomrec3d(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

if (isempty(handles.data.markersets))
    warndlg('No markerset available.', 'Save em-file', 'modal');
    return;
end;

%if (isempty(handles.data.markersets(handles.data.sel_markerset).filename) || ~exist(handles.data.markersets(handles.data.sel_markerset).filename, 'file'))
%    errordlg('You have to save the markerset into an em-file before exporting to a tom_Rec3d-project. Use export->to EM-file.');
%    return;
%end;


load TESTNAME
Rec3dProject_orig = Rec3dProject;
clear('Rec3dProject');

ProjectionsName = cell(1, length(handles.data.filenames));
ProjectionsPath = cell(1, length(handles.data.filenames));

 
i = 1;
j = [];
base = [];
for (i=1:length(handles.data.filenames))
    n = regexp(handles.data.filenames{i}, ['^(?<pathname>.*)' filesep() '(?<filename>[^' filesep() ']*)$'], 'names');
    if (isempty(n))
        errordlg(['Invalid filename for image #' num2str(i) ': ' handles.data.filenames{i}]);
        return;
    end;
    
    if (~exist(fullfile(n.pathname, n.filename)))
        errordlg(['tom_Rec3d expects all images to exist! (image #' num2str(i) ': ' handles.data.filenames{1} ').'], 'Export to tom_Rec3d', 'modal');
        return;
    end;

    if (i == 1)
        if (~exist(n.pathname))
            errordlg(['The path for image #1 does not exist (' handles.data.filenames{1} ').'], 'Export to tom_Rec3d', 'modal');
            return;
        end;
        ProjectionsPath1 = n.pathname;
    else
        if (~strcmp(ProjectionsPath1, n.pathname))
            errordlg(['The path for image #' num2str(i) ' is different from the first. tom_Rec3d expects all em-image in the same direcotry'], 'Export to tom_Rec3d', 'modal');
            return;
        end;
    end;
    
    n.filename_div = regexp(n.filename, '^(?<name>.*[^1-9])(?<num>[1-9][0-9]*)(?<ext>.*)$', 'names');
    
    if (isempty(n.filename_div))
        errordlg(['The filename of image #' num2str(i) ' (' n.filename ') has not the expected format for tom_Rec3d'], 'Export to tom_Rec3d', 'modal');
        return;
    end;
    
    if (str2double(n.filename_div.num) ~= i)
        errordlg(['The filename of image #' num2str(i) ' (' n.filename ') has not the expected format for tom_Rec3d: wrong numbering!'], 'Export to tom_Rec3d', 'modal');
        return;
    end;
    
    if (i == 1)
        ProjectionsName1 = n.filename_div.name;
        ProjectionsExt1 = n.filename_div.ext;
    else
        if (~strcmp(n.filename_div.name, ProjectionsName1))
            errordlg(['The filename of image #' num2str(i) ' (' n.filename ') has not the expected format for tom_Rec3d: different base name!'], 'Export to tom_Rec3d', 'modal');
            return;
        end;
        if (~strcmp(n.filename_div.ext, ProjectionsExt1))
            errordlg(['The filename of image #' num2str(i) ' (' n.filename ') has not the expected format for tom_Rec3d: different extension!'], 'Export to tom_Rec3d', 'modal');
            return;
        end;
    end;
    
    ProjectionsName{i} = n.filename;
    ProjectionsPath{i} = n.pathname;
end;

[Tiltangle, imsize] = getTiltangles(handles.data.emheader, false); 

i = find(~isfinite(Tiltangle));
if (~isempty(i))
    errordlg(['Not all em-images exist. This is mandatory for tom_Rec3d (image #' num2str(i(1)) ': ' handles.data.filenames{i(1)} ').'], 'Export to tom_Rec3d', 'modal');
    return;
end;
i = find(any(imsize(1,1) ~= imsize, 1));
if (~isempty(i))
    errordlg(['Not em-images must have the same dimension! (#' num2str(i(1)) ': ' handles.data.filenames{i(1)} ' is not ' num2str(imsize(1,1)) ').'], 'Export to tom_Rec3d', 'modal');
    return;
end;

ms = handles.data.markersets(handles.data.sel_markerset);
fversn='';
[fpathstr, fname, fext] = fileparts(ms.filename); 


ms_markerset = ms.markerset;

i = ms_markerset < 0.5 | ms_markerset > (imsize(1,1) + 0.5);
if (any(i(:)))
    warndlg(['There are some markerpositions outside the image-dimensions. They are not exported!'], 'Export to tom_Rec3d', 'modal');
end;
ms_markerset(:, ~all(ms_markerset >= 0.5 & ms_markerset <= (imsize(1,1) + 0.5), 1)) = -1;
ms_markerset_perm = permute(ms_markerset, [3 2 1]);
ms_markerset_orig = cat(1, zeros(1, size(ms_markerset,2), size(ms_markerset,3)), ms_markerset);
ms_markerset_orig(1,:,1) = Tiltangle;
ms_markerset_orig = struct('Value', ms_markerset_orig, 'Header', struct);

Rec3dProject = struct( ...
    'ProjectName', handles.data.markersets(handles.data.sel_markerset).name, ...
    'ProjectStatus', 'loaded', ...
    'ProjectVersion', '1.0', ...
    'PROJECT', struct(...
        'TiltingGeometry', 'singleaxis', ...
        'ProjectionsName1', ProjectionsName1, ...
        'ProjectionsPath1', [ProjectionsPath1 filesep()], ...
        'ProjectionsExt1', ProjectionsExt1, ...
        'MarkerfileName1', [fname, fext, fversn], ...
        'MarkerfilePath1', fpathstr, ...
        'MarkerfileExt1', '', ...
        'NumOfProj1', length(handles.data.filenames), ...
        'NumOfProjOriginal1', length(handles.data.filenames), ...
        'RefProj1', find(abs(Tiltangle) == min(abs(Tiltangle)), 1), ...
        'Tiltangles1', Tiltangle, ...
        'MarkersOnProj1', ms_markerset_perm, ...
        'MarkersOnProjOriginal1', ms_markerset_perm, ...
        'Markerfile1', ms_markerset_orig, ...
        'MarkerfileOriginal1', ms_markerset_orig, ...
        'ProjectionsName2', '', ...
        'ProjectionsPath2', '', ...
        'ProjectionsExt2', '', ...
        'MarkerfileName2', '', ...
        'MarkerfilePath2', '', ...
        'MarkerfileExt2', '', ...
        'NumOfProj2', 0, ...
        'NumOfProjOriginal2', 0, ...
        'RefProj2', 0, ...
        'Tiltangles2', 0, ...
        'MarkersOnProj2', 0, ...
        'MarkersOnProjOriginal2', 0, ...
        'Markerfile2', 0, ...
        'MarkerfileOriginal2', 0, ...
        'NumOfProj', length(handles.data.filenames), ...
        'NumOfMarkers', size(ms.markerset, 3), ...
        'Imdim', imsize(1,1) ...
        ), ...
    'ALIGNMENT', struct( ...
        'AlignmentMethod', 'rigidbody', ...
        'ReferenceMarker', 1, ...
        'Origin1', [0 0 0], ...
        'm3d1', 0, ...
        'Tiltaxis1', 0, ...
        'tx1', 0, ...
        'ty1', 0, ...
        'isoscale1', 1, ...
        'WarpDone1', 'no', ...
        'WarpAlignment1', 0, ...
        'Origin2', [0 0 0], ...
        'm3d2', 0, ...
        'Tiltaxis2', 0, ...
        'tx2', 0, ...
        'ty2', 0, ...
        'isoscale2', 1, ...
        'WarpDone2', 'no', ...
        'WarpAlignment2', 0, ...
        'RotMatrix', 0, ...
        'Psi', 0, ...
        'Theta', 0, ...
        'Phi', 0 ...
        ), ...
    'ALGRESIDUALS', struct( ...
        'ResidualMatrix1', 0, ...
        'Sigma1', 0, ...
        'AveragePerProjection1', 0, ...
        'AveragePerMarker1', 0, ...
        'ResidualMatrix2', 0, ...
        'Sigma2', 0, ...
        'AveragePerProjection2', 0, ...
        'AveragePerMarker2', 0, ...
        'EulerAnglesResidual', 0, ...
        'MaximumResidual', 0, ...
        'ResidualSpheres', 0, ...
        'AverageResidualSphere', 0 ...
        ), ...
    'NAMEPATHEXT', struct(...
        'ReconstructionName', '', ...
        'ReconstructionPath', '', ...
        'ReconstructionExt', '', ...
        'TempFilesName', '', ...
        'TempFilesPath', '', ...
        'TempFilesExt', '' ...
        ), ...
    'PARAMETER', struct( ...
        'ProjectionsName', {ProjectionsName}, ...
        'ProjectionsNameOriginal', {ProjectionsName}, ...
        'ProjectionsPath', {ProjectionsPath}, ...
        'ProjectionsPathOriginal', {ProjectionsPath}, ...
        'Tiltangles', 0, ...
        'Tiltaxis', 0, ...
        'ProjDir', 0, ...
        'tx', 0, ...
        'ty', 0, ...
        'isoscale', 1 ...
        ), ...
    'VOLUME', struct( ...
        'SizeX', imsize(1,1), ...
        'SizeY', imsize(1,1), ...
        'SizeZ', imsize(1,1), ...
        'PreBinning', 0, ...
        'PostBinning', 0 ...
        ), ...
    'METHOD', struct(...
        'ReconstructionMethod', 'WBP', ...
        'Normalization', 'phase', ...
        'SmoothBorders', 0, ...
        'AlignTiltaxis', 'ProjDir', ...
        'Handedness', 0, ...
        'ApplyWeighting', 'on', ...
        'WeightingMethod', 'exact', ...
        'ObjectThickness', imsize(1,1), ...
        'Taper', 'off', ...
        'Iterations', 12, ...
        'Relaxation', 0.0500, ...
        'Pathlength', 'ones' ...
        ), ...
    'FILTER', struct( ...
        'Apply', 1, ...
        'Value', struct( ...
            'times', 1, ...
            'low', 0, ...
            'high', 1024, ...
            'smooth', 0, ...
            'space', 'real', ...
            'method', 'quadr', ...
            'radius', 0 ...
            ), ...
        'Type', 'bandpass' ...
        ), ...
    'DETAIL', struct(...
        'DetailMode', 'off', ...
        'OverviewSizeZ', 0, ...
        'OverviewPreBinning', 0, ...
        'OverviewPostBinning', 0, ...
        'NumberOfDetails', 0, ...
        'DetailCoordinates', [], ...
        'av3_alignstruct', [] ...
        ), ...
    'PARALLEL', struct( ...
        'jobmanager', 'default_jobmanager', ...
        'packageloss', 0, ...
        'number_of_tasks', 1, ...
        'workers', struct( ...
            'min', 1, ...
            'max', 16 ...
            ), ...
        'timeout', 3600, ...
        'restart_workers', 0 ... 
        ), ...
    'SIRTRESIDUALS', struct(    ...
        'EuclidianDistancePerProjection', 0, ...
        'EuclidianDistancePerIteration', 0, ...
        'dProjVal', 0, ...
        'DifVolVal', 0, ...
        'RecVolVal', 0 ...
        ));

filename = fullfile(fpathstr, [fname '_Rec3d.mat']);                             
[filename, pathname, filterindex] = uiputfile({'*.mat','MATLAB-files (*.mat)'; ...
                                  '*',  'All Files (*)'}, ...
                                 'Export tom_Rec3d-Project', filename);                         
if (isequal(filename, 0) || isequal(pathname, 0))
    return;
end;
 
uri = fullfile(pathname, filename);
save(uri, 'Rec3dProject');

                             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_saveemfile(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

if (isempty(handles.data.markersets))
    warndlg('No markerset available.', 'Save em-file', 'modal');
    return;
end;

ms = handles.data.markersets(handles.data.sel_markerset);

msize = size(ms.markerset);
if (length(msize)==2)
    msize(3) = 1;
end;
if (msize(3) == 0)
    warndlg('The selected markerset contains no markers.', 'Save em-file', 'modal');
    return;
end;


filename = ms.filename;
if (isempty(filename))
    filename = fullfile(handles.workingdir, 'markerfile.em');
elseif (~exist(filename, 'file'))
    %[fpathstr, fname, fext, fversn] = fileparts(filename);
    fversn='';
    [fpathstr, fname, fext] = fileparts(filename);
    if (~exist(fpathstr,'dir'))
        filename = fullfile(handles.workingdir, [fname fext, fversn]);
    end;
end;
[filename, pathname, filterindex] = uiputfile({'*.em','EM-files (*.em)'; ...
                                  '*',  'All Files (*)'}, ...
                                 'Select markerfile', filename);

if (isequal(filename, 0) || isequal(pathname, 0))
    return;
end;
 
uri = fullfile(pathname, filename);

idx = ~all(ms.markerset > 0, 1);
ms2 = zeros(3, msize(2), msize(3));

for (i=1:msize(2))
    if (isempty(handles.data.emheader{i}))
        ms2(1, i, 1) = 0;
        %idx(1, i, :) = true;
    else
        ms2(1, i, 1) = handles.data.emheader{i}.Tiltangle;
    end;    
end;
ms2([2 3], 1:msize(2), 1:msize(3)) = double(ms.markerset);
ms2([2 3], idx) = -1;


ms2 = tom_emheader(ms2);
tom_emwrite(uri, ms2);

handles.workingdir = pathname;

handles.data.markersets(handles.data.sel_markerset).filename = uri;
guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_loademfile(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

if (isempty(handles.data.filenames))
    errordlg('Insert first the images to the project!', 'Loading Markerfile', 'modal');
    return;
end;

filename = get(handles.mainmenu_markersets_load_emfile, 'Userdata');
if (isempty(filename))
    filename = 'markerfile.em';
end;
filename = getExistingFilenameFromString(filename, handles.workingdir);


[filename, pathname] = uigetfile({'*.em','EM-files (*.em)'; ...
                                  '*.*',  'All Files (*.*)'}, ...
                                 'Select markerfile', filename);

if (isequal(filename, 0) || isequal(pathname, 0))
    return;
end;
uri = fullfile(pathname, filename);

try
    mfile = tom_emread(uri);
catch
    lerror = lasterror();
    errordlg(['Error reading markerfile: ' lerror.message], 'Loading Markerfile', 'modal');
    return;
end;


sel_markerset = length(handles.data.markersets) + 1;

handles.data.sel_markerset = sel_markerset;
if (length(filename)>=4 && strcmp(filename(end-2:end),'.em'))
    handles.data.markersets(sel_markerset).name = filename(1:end-3);
else
    handles.data.markersets(sel_markerset).name = filename;
end;
handles.data.markersets(sel_markerset).filename = uri;
handles.data.markersets(sel_markerset).markerset = [];
handles.data.markersets(sel_markerset).findmatchesparam = struct();
[error, handles] = setMarkerSet(handles, mfile);
if (error)
    return;
end;
handles.data.markersets(sel_markerset).findmatchesparam.comparelist = zeros(0,3);
msize = size3(handles.data.markersets(sel_markerset).markerset);
handles.data.markersets(sel_markerset).findmatchesparam.matchchain = nan(1, msize(2), msize(3));
handles.data.markersets(sel_markerset).findmatchesparam.matchchain(all(isfinite(handles.data.markersets(sel_markerset).markerset), 1)) = 0;

handles.workingdir = pathname;
set(handles.mainmenu_markersets_load_emfile, 'Userdata', filename);
handles.workingdir = pathname;

handles.data.fileinfo = getFileInfo_info(handles, 1:length(handles.data.filenames));

sel_images = getListboxValue(handles.listbox_images);
sel_markers = 1:size(handles.data.markersets(sel_markerset).markerset, 3);

updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_info_image_resetdisplay(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

cla(handles.axes_main);

userdata = get(handles.axes_main, 'Userdata');
userdata.force = true;
set(handles.axes_main, 'Userdata', userdata);

updateSelection(handles, ...
    getListboxValue(handles.listbox_images), ...
    getListboxValue(handles.listbox_markers));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads a markerset from the base workspace.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_loadworkspace(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

lfilenames = length(handles.data.filenames);

if (lfilenames < 1)
    errordlg('Insert first the images to the project!', 'Loading Markerfile', 'modal');
    return;
end;


vars = evalin('base', 'whos');
idx = false(1, length(vars));

for (i=1:length(idx))
    switch (vars(i).class)
        case {'struct'}
            if (all(vars(i).size == [1 1]))
                tmpstruct = evalin('base', vars(i).name);
                if (isfield(tmpstruct, 'Header') && isfield(tmpstruct,'Value') && isnumeric(tmpstruct.Value))
                    msize = size(tmpstruct.Value);
                    if (msize(1)>=2 && msize(1)<=12 && msize(2)==lfilenames)
                        idx(i) = true;
                    end;
                end;
                clear('tmpstruct');
            end;
        case {'double', 'single', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'}
            msize = vars(i).size;
            if (msize(1)>=2 && msize(1)<=12 && msize(2)==lfilenames && length(msize)>=2 && length(msize)<=3)
                idx(i) = true;
            end;
    end;       
end;

if (~any(idx))
    warndlg(['In the base workspace there were no variable found, which can be used as markerset (either an em-structure or a (2-12)x' num2str(lfilenames) 'xM-array).'], 'Loading markerset from workspace', 'modal')
    return;
end;
vars = vars(idx);


[selection,ok] = listdlg('PromptString', 'Select the workspace variable',...
                         'SelectionMode','single',...
                         'Name', 'Load markerset', ...
                         'ListString', { vars(:).name }, ...
                         'Listsize', [300 200]);
if (~ok)
    return;
end;

try
    ms = evalin('base', vars(selection).name);
catch
    errordlg(['An error occured loading the variable ' vars(selection).name ' from the base workspace']);
    return;
end;

if (~isstruct(ms))
    ms = double(ms);
end;

sel_markerset = length(handles.data.markersets) + 1;

handles.data.sel_markerset = sel_markerset;
handles.data.markersets(sel_markerset).name = vars(selection).name;
handles.data.markersets(sel_markerset).filename = fullfile(pwd(), [vars(selection).name '.em']);
handles.data.markersets(sel_markerset).markerset = [];
handles.data.markersets(sel_markerset).findmatchesparam = struct();
[error, handles] = setMarkerSet(handles, ms);
if (error)
    return;
end;
handles.data.markersets(sel_markerset).findmatchesparam.comparelist = zeros(0,3);
msize = size3(handles.data.markersets(sel_markerset).markerset);
handles.data.markersets(sel_markerset).findmatchesparam.matchchain = nan(1, msize(2), msize(3));
handles.data.markersets(sel_markerset).findmatchesparam.matchchain(all(isfinite(handles.data.markersets(sel_markerset).markerset), 1)) = 0;



handles.data.fileinfo = getFileInfo_info(handles, 1:length(handles.data.filenames));

sel_images = getListboxValue(handles.listbox_images);
sel_markers = 1:size(handles.data.markersets(sel_markerset).markerset, 3);

updateListboxDisplay(handles);

updateSelection(handles, sel_images, sel_markers);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [error, handles] = setMarkerSet(handles, ms)

error = true;

if (isstruct(ms))
    if (~isfield(ms, 'Value'))
        return;
    end;
    ms = ms.Value;
end;

msize = size3(ms);

if (msize(2) ~= length(handles.data.filenames))
    errordlg(['The markerfile contains data for ' num2str(msize(2)) ' images. Load first all images in the appropriate order!'], 'Loading Markerfile', 'modal');
    return;   
end;
 
if (msize(1) ~= 2)
    tiltangles = getTiltangles(handles.data.emheader);
    errtxt = {};
    for (i=1:msize(2))
        if (isnan(tiltangles(i)))
            errtxt{end+1} = ['The ' num2str(i) 'th image could not be loaded successfully. You should load first the images!']; 
        elseif (abs(ms(1,i,1) - tiltangles(i)) > 1e-5)
            errtxt{end+1} = ['The ' num2str(i) 'th image has an other tiltangle then in the markerfile (' num2str(handles.data.emheader{i}.Tiltangle) ' vs. ' num2str(ms(1,i,1)) ')!']; %#ok<AGROW>
        end;
    end;

    if (~isempty(errtxt))
        disp(strvcat(errtxt)); %#ok<VCAT>
        answ = questdlg('There were warnings during analysing the markerfile. You can ignore them or abort (see the matlab command window for warnings).', 'Loading Markerfile', 'Ignore', 'Cancel', 'Ignore');
        if (~strcmp(answ, 'Ignore'))
            return;
        end;
    end;
    
    ms = ms([2 3], 1:msize(2), 1:msize(3));
end;


ms(ms<0) = nan;

handles.data.markersets(handles.data.sel_markerset).markerset = ms;


error = false;









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = getMarkerName(handles, sel_markerset, sel_marker)

s = [num2str(sel_marker) ' (in ' num2str(sum(all(handles.data.markersets(sel_markerset).markerset(:,:,sel_marker) > 0, 1))) ' images)'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = getImageName(handles, sel_image)

%[fp_path, fp_name, fp_ext, fp_versn] = fileparts(handles.data.filenames{sel_image});

[fp_path, fp_name, fp_ext] = fileparts(handles.data.filenames{sel_image});

fp_versn='';

if (isempty(handles.data.emheader{sel_image}))
    s = [sprintf('%3d: ', sel_image) fp_name fp_ext fp_versn '*'];
else
    s = [sprintf('%3d: ', sel_image) fp_name fp_ext fp_versn ' (' num2str(handles.data.emheader{sel_image}.Tiltangle) ')'];
end;

lmarkersets = length(handles.data.markersets);
if (lmarkersets > 0)
    n = 0;
    nsel = 0;
    for (i=1:lmarkersets)
        if (~isempty(handles.data.markersets(i).markerset))
            if (i == handles.data.sel_markerset)
                nsel = sum(all(handles.data.markersets(i).markerset(:,sel_image,:)>0, 1));
                n = n + nsel;
            else
                n = n + sum(all(handles.data.markersets(i).markerset(:,sel_image,:)>0, 1));
            end;
        end;
    end;
    s = [s ' (' num2str(nsel) '/' num2str(n) 'm)'];
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializes the guifield listbox_markers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateListboxDisplayLockImages(hLB, hPN, lmarkersets)
posLS = get(hLB, 'Position');
posPN = get(hPN, 'Position');
if ((lmarkersets==0 && posLS(1)<posPN(1)) || ...
    (lmarkersets>=1 && posLS(1)>posPN(1)))            
    set(hLB, 'Position', posPN);
    set(hPN, 'Position', posLS);
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializes the guifield listbox_markers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateListboxDisplay(handles)


lmarkersets = length(handles.data.markersets);
if (lmarkersets < 2)
    updateListboxDisplayLockImages(handles.listbox_images, handles.uipanel_listbox_images_position, lmarkersets);
    if (lmarkersets == 0)
        handles.data.clipboard.markerset = nan(2, 0, 0);    
    end;
end;

% First the image parameters...
if (isempty(handles.data.filenames))
    s = {};
else
    s = { handles.data.fileinfo(:).dispname };
end;
setListboxValue(handles.listbox_images, []);
set(handles.listbox_images, ...
    'String', s, ...
    'ListBoxTop', max(1, min(get(handles.listbox_images, 'ListBoxTop'), length(handles.data.filenames))));


% Then the markersets...
s = {};
for (i=1:lmarkersets)
    nel = numel(handles.data.markersets(i).markerset);
    if (nel == 0)
        s{i} = [num2str(i,2) ': ' handles.data.markersets(i).name];
    else
        s{i} = [num2str(i,2) ': ' handles.data.markersets(i).name ' (' num2str(size(handles.data.markersets(i).markerset, 3)) ', ' ...
                sprintf('%.1f%%', 200*sum(sum(all(handles.data.markersets(i).markerset>0, 1)))/nel) ')'];
    end;
end;


set(handles.listbox_markersets, ...
    'String', s, ...
    'Value', handles.data.sel_markerset, ...
    'ListBoxTop', max(1, min(get(handles.listbox_markersets, 'ListBoxTop'), lmarkersets)));


s_markerlist = {};
if (lmarkersets >= 1)
    markerset = handles.data.markersets(handles.data.sel_markerset);
    msize = size(markerset.markerset);
    if (length(msize)==2)
        msize(3) = 1;
    end;    

    for (i=1:msize(3))
        s_markerlist{i} = getMarkerName(handles, handles.data.sel_markerset, i);
    end;
else
    msize = [2 0 0];
end;

setListboxValue(handles.listbox_markers, []);
set(handles.listbox_markers, 'String', s_markerlist, 'ListBoxTop', max(1, min(get(handles.listbox_markers, 'ListBoxTop'), msize(3))));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function KeyPressFcn_figure_markGui(hObject, eventdata) 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function KeyPressFcn_listbox_markersets(hObject, eventdata) 


handles = guidata(hObject);

if (isempty(eventdata.Modifier) && strcmp(eventdata.Key, 'delete') && ~isempty(handles.data.markersets))    
    Callback_mainmenu_markersets_delete(hObject, [], handles);    
elseif (isempty(eventdata.Modifier) && strcmp(eventdata.Key, 'insert') && ~isempty(handles.data.filenames))
    Callback_mainmenu_markersets_new(hObject, [], handles);    
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_listbox_markersets(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

if (isempty(handles.data.markersets))
    return;
end;

guidata_updated = false;

sel_markerset = get(handles.listbox_markersets, 'Value');

selectionType = get(handles.figure_markGui, 'SelectionType');
switch (selectionType) 
    case {'open'}
        newname = inputdlg(['Enter the new name for the ' num2str(sel_markerset) '. markerset.'], 'Change markerset name', 1, { handles.data.markersets(sel_markerset).name }, 'on'); 
        if (length(newname)==1 && ~strcmp(newname{1}, handles.data.markersets(sel_markerset).name))
            handles.data.markersets(sel_markerset).name = newname{1};
            guidata_updated = true;
        end;
    otherwise
end;


if (sel_markerset ~= handles.data.sel_markerset)
    
    msize3 = size(handles.data.markersets(sel_markerset).markerset, 3);
    
    forceidx = getIndexRangeFromString(get(handles.edit_infomarker_selnumber, 'String'), msize3);


    sel_markersold = getListboxValue(handles.listbox_markers);
    if (get(handles.checkbox_infomarker_selnumber, 'Value'))
        sel_markers = forceidx;
    else
        if (size(handles.data.markersets(handles.data.sel_markerset).markerset, 3) == size(handles.data.markersets(sel_markerset).markerset, 3))
            sel_markers = [forceidx, setdiff(sel_markersold, forceidx)];            
        else
            sel_markers = [forceidx, setdiff(1:msize3, forceidx)];            
        end;
    end;

    if (~isempty(sel_markersold) && sel_markersold(1)<=msize3 && sel_markersold(1)>0)
        i = find(sel_markersold(1) == sel_markers);
        if (~isempty(i))
            sel_markers = sel_markers([i, setdiff(1:length(sel_markers), i)]);
        end;
    end;
    
    handles.data.sel_markerset = sel_markerset;
    
    
    
    
%new_sel_markers = getIndexRangeFromString(get(handles.edit_infomarker_selnumber, 'String'), msize(3));    
 %   getIndexRangeFromString(get(handles.edit_infomarker_selnumber, 'String'), size(handles.data.markersets(sel_markerset).markerset, 3))
    
    guidata_updated = true;
else
    sel_markers = getListboxValue(handles.listbox_markers);
end;

if (guidata_updated)
    handles.data.fileinfo = getFileInfo_info(handles, 1:length(handles.data.filenames));
    sel_images = getListboxValue(handles.listbox_images);
    updateListboxDisplay(handles);
    
    updateSelection(handles, sel_images, sel_markers);
    drawnow;
    %set(handles.listbox_markersets, 'BackgroundColor', 'White');    
    
end;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_listbox_images_marked(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

sel_markerimages = get(handles.listbox_images_marked, 'Value');
if (isempty(sel_markerimages))
    return;
end;

sel_markers = getListboxValue(handles.listbox_markers);
sel_markers_sort = sort(sel_markers);

%index = find(any(all(handles.data.markersets(handles.data.sel_markerset).markerset(:,:,sel_markers_sort) > 0, 1), 3));
index = get(handles.listbox_images_marked, 'Userdata');

sel_images = getListboxValue(handles.listbox_images);
new_sel_images = index(sel_markerimages);

updateSelection(handles, [setdiff(new_sel_images, sel_images), sel_images(ismember(sel_images, new_sel_images))], sel_markers);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_listbox_images_markers(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

sel_imagemarkers = get(handles.listbox_images_markers, 'Value');
if (isempty(sel_imagemarkers))
    return;
end;
sel_images = getListboxValue(handles.listbox_images);
sel_images_sort = sort(sel_images);

%index = find(all(all(handles.data.markersets(handles.data.sel_markerset).markerset(:,sel_images_sort,:) > 0, 1), 2))';
index = get(handles.listbox_images_markers, 'Userdata');

sel_marker = getListboxValue(handles.listbox_markers);
new_sel_marker = index(sel_imagemarkers);


updateSelection(handles, sel_images, [setdiff(new_sel_marker, sel_marker), sel_marker(ismember(sel_marker, new_sel_marker))]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_initms(hObject, eventdata, handles) %#ok<INUSL,DEFNU>

if (isempty(handles.data.filenames))
    errordlg('Insert first the images to the project!', 'Init markerfile', 'modal');
    return;
end;

if (~isfield(handles.data, 'initMarkerGuiData') || ~isfield(handles.data.initMarkerGuiData, 'sel_images') || isempty(handles.data.initMarkerGuiData.sel_images))
    handles.data.initMarkerGuiData.sel_images = sort(getListboxValue(handles.listbox_images));
end;

[markerset, initMarkerGuiData] = tom_markInitMarkerGui(handles.data.filenames, handles.data.initMarkerGuiData, 'wait', 'WindowStyle', 'modal');
%[markerset, initMarkerGuiData] = tom_markInitMarkerGui(handles.data.filenames, handles.data.initMarkerGuiData, 'wait', 'WindowStyle', 'normal');


if (isnumeric(markerset))
    handles.data.initMarkerGuiData = initMarkerGuiData;


    sel_markerset = length(handles.data.markersets) + 1;

    handles.data.sel_markerset = sel_markerset;
    msize = size3(markerset);
    
    i = 1;
    if (sel_markerset == 1)
        s = ['< init ' num2str(i) ' >'];
    else
        while (i<=sel_markerset)
            s = ['< init ' num2str(i) ' >'];
            if (any(strcmp({ handles.data.markersets(:).name }, s)))
                i = i+1;
            else
                break;
            end;
        end;
    end;
    handles.data.markersets(sel_markerset).name = s;
    handles.data.markersets(sel_markerset).filename = ['init' num2str(i) '.em'];
    handles.data.markersets(sel_markerset).markerset = markerset;
    handles.data.markersets(sel_markerset).findmatchesparam = struct('comparelist', {zeros(0, 3)}, 'matchchain', {nan(1, msize(2), msize(3))});
    handles.data.markersets(sel_markerset).findmatchesparam.matchchain(all(isfinite(markerset), 1)) = 0;
    

    handles.data.fileinfo = getFileInfo_info(handles, 1:length(handles.data.filenames));

    sel_images = getListboxValue(handles.listbox_images);
    sel_markers = 1:size(markerset, 3);

    updateListboxDisplay(handles);
    updateSelection(handles, sel_images, sel_markers);

end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_axes_main(hObject, eventdata, handles)

sel_images = getListboxValue(handles.listbox_images);
sel_markers = getListboxValue(handles.listbox_markers);





point1 = get(handles.axes_main, 'CurrentPoint');
rbbox;
point2 = get(handles.axes_main,'CurrentPoint'); 

point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);   
p2 = p1 + offset;

markerset = handles.data.markersets(handles.data.sel_markerset).markerset;

msize = size(markerset);
if (length(msize) == 2)
    msize(3) = 1;
end;

limits = [get(handles.axes_main, 'XLim'); get(handles.axes_main, 'YLim')];


if (all(offset ./ (limits(:,2)-limits(:,1))' < 0.025))
    if (get(handles.checkbox_infomarker_set_marker, 'Value'))
        selectiontype = get(handles.figure_markGui, 'SelectionType');

        userdata = get(handles.checkbox_infomarker_set_marker, 'Userdata');
        
        switch (selectiontype)
            case {'alt', 'extend', 'open'}
                if (strcmp(selectiontype, 'extend'))
                    userdata.direction = userdata.direction * -1;
                    set(handles.checkbox_infomarker_set_marker, 'Userdata', userdata, 'BackgroundColor', getDirectionLabelColor(userdata.direction));
                end;
                sel_images = max(1,min(msize(2), sel_images(1) + userdata.direction));
            case {'normal'}
                disp(['reset marker #' num2str(sel_markers(1)) ' in image #' num2str(sel_images(1)) ' from [' num2str(markerset(1,sel_images(1), sel_markers(1))) ', ' num2str(markerset(2,sel_images(1), sel_markers(1))) '] to [' num2str(p1(1)) ', ' num2str(p1(2)) ']']);
                handles.data.markersets(handles.data.sel_markerset).markerset(1:2,sel_images(1), sel_markers(1)) = p1;
            otherwise
        end;
        %setMarkerPreview(handles, sel_images, sel_markers);
    else
        if (strcmp(get(handles.figure_markGui, 'SelectionType'), 'alt'))
            consider_idx = 1:msize(3);
        else
            consider_idx = sel_markers;
        end;
        markers_idx = find(all(markerset(:,sel_images(1),consider_idx) > 0, 1));

        distance = markerset(:,sel_images(1),consider_idx(markers_idx));
        distance(1,:,:) = distance(1,:,:) - p1(1);
        distance(2,:,:) = distance(2,:,:) - p1(2);
        distance = sum(distance .^ 2, 1);
        %distance = sqrt(distance);

        [mindist, minindex] = min(distance);

        sel_marker = consider_idx(markers_idx(minindex));
        if (isempty(sel_marker))
            sel_markers = [];
        else
            sel_markers = [sel_marker sel_markers(sel_markers~=sel_marker)];
        end;
    end;
else
    sel_markers_new = find(markerset(1,sel_images(1),:)>=p1(1) & markerset(1,sel_images(1),:)<=p2(1) & markerset(2,sel_images(1),:)>=p1(2) & markerset(2,sel_images(1),:)<=p2(2))';
    if (strcmp(get(handles.figure_markGui, 'SelectionType'), 'alt'))
        sel_markers = [sel_markers_new(~ismember(sel_markers_new, sel_markers)) sel_markers];
    else
        sel_markers = sel_markers_new;
    end;
end;

updateSelection(handles, sel_images, sel_markers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setMarkerPreview(handles, varargin)

markerset = handles.data.markersets(handles.data.sel_markerset).markerset;
msize = size(markerset);
if (length(msize)==2)
    msize(3) = 1;
end;

hide = false;
if (length(varargin) == 1 && strcmp(varargin{1},'hide'))
    hide = true;
else    
    sel_images = varargin{1};
    sel_markers = varargin{2};
end;
    

tag{1,1} = 'prev'; tag{1,2} = -1;
tag{2,1} = 'next'; tag{2,2} =  1;
currentfigure = get(0, 'CurrentFigure');

xlim = get(handles.axes_main, 'XLim');
ylim = get(handles.axes_main, 'YLim');
for (i=1:length(tag))
    h = findobj('Type', 'figure', 'Tag', ['figure_tom_markGui_setmarker_' tag{i,1}]);
    if (length(h)~=1 || ~ishandle(h))
        delete(h);
        h = figure('Tag', ['figure_tom_markGui_setmarker_' tag{i,1}]);
        uicontrol('Parent', h, 'Style','Pushbutton', 'Visible', 'off', 'Tag', ['clamp_dummy']);
    end;
    if (hide)
        set(h, 'Visible', 'off');
    else
        sel_image_preview = sel_images(1) + tag{i,2};
        if (sel_image_preview >= 1 && sel_image_preview <= msize(2))
            set(h, 'Name', ['tom_markGui: set marker #' num2str(sel_markers(1)) ', image #' num2str(sel_image_preview)], 'Visible', 'on');
            set(0, 'CurrentFigure', h);
            currentaxes = gca;
            set(currentaxes, 'Units', 'pixels')
            
            hposition = get(h, 'Position');
            borderdistance = 25;
            set(currentaxes, 'Position', [borderdistance,borderdistance,hposition(3:4) - borderdistance]);

            userdata = get(currentaxes, 'UserData');
            if (isempty(userdata))
                userdata = struct('image', '', 'markers', []);
                userdata.uri = '';
                userdata.wrect = ([0 1; 0 1]);
                userdata.force = true;
                set(currentaxes, 'UserData', userdata);
            end;
            
            currxlim = get(currentaxes,'XLim');
            currylim = get(currentaxes,'YLim');
            if ((currxlim(2)-currxlim(1)) ~= (xlim(2)-xlim(1)) || ...
                (currylim(2)-currylim(1)) ~= (ylim(2)-ylim(1)))
                set(currentaxes, 'XLim', xlim);
                set(currentaxes, 'YLim', ylim);
            end;
            handles = updateImageDisplay(handles, sel_image_preview, sel_markers, currentaxes);
            
        else
            set(h, 'Visible', 'off');
        end;
    end;
end;
set(0, 'CurrentFigure', currentfigure);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_checkbox_infomarker_set_marker(hObject, eventdata, handles)


if (get(handles.checkbox_infomarker_set_marker, 'Value'))
    userdata = get(handles.checkbox_infomarker_set_marker, 'Userdata');
    if (isempty(userdata))
        userdata.direction = 1;
    end;
    set(handles.checkbox_infomarker_set_marker, 'Userdata', userdata, 'BackgroundColor', getDirectionLabelColor(userdata.direction));
    %sel_images = getListboxValue(handles.listbox_images);
    %sel_markers = getListboxValue(handles.listbox_markers);
    %setMarkerPreview(handles, sel_images, sel_markers)
else
    set(handles.checkbox_infomarker_set_marker, 'BackgroundColor', get(handles.uipanel_infomarker_panel, 'Background'));
    %setMarkerPreview(handles, 'hide')
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_markers_remove(hObject, eventdata, handles)

sel_markers = getListboxValue(handles.listbox_markers);
sel_images = getListboxValue(handles.listbox_images);

msize = size3(handles.data.markersets(handles.data.sel_markerset).markerset);

if (isempty(sel_markers) || isempty(sel_images) || ~any(any(all(handles.data.markersets(handles.data.sel_markerset).markerset(1:2, sel_images, sel_markers) > 0, 1))))
    return;
end;

matchchain_cascade = false;
matchchain = handles.data.markersets(handles.data.sel_markerset).findmatchesparam.matchchain(1,setdiff(1:msize(2),sel_images),sel_markers);
for (i=1:length(sel_markers))
    if (any(ismember(matchchain(1,:,i), sel_images)))
        matchchain_cascade = true;            
        break;
    end;
end;

if (matchchain_cascade)
    answ = questdlg(['Some of the selected ' num2str(length(sel_markers)) ' markers have known positions in images other then the seleted ' num2str(length(sel_images)) ', which derive from a matching based on those you are about to delete. Only delete the selected positions or cascade?'], 'Remove markers', 'Cascade', 'Only selected', 'Cancel', 'Cascade');
    if (~any(strcmp(answ, {'Cascade', 'Only selected'})))
        return;
    end;
    matchchain_cascade = strcmp(answ, 'Cascade');
else
    answ = questdlg(['Are you sure, that you want to remove the selected ' num2str(length(sel_markers)) ' markers in the selected ' num2str(length(sel_images)) ' images?'], 'Remove markers', 'Yes', 'No', 'Yes');
    if (~strcmp(answ, 'Yes'))
        return;
    end;
end;

setzeromask = false(1, msize(2), msize(3));
setzeromask(1,sel_images, sel_markers) = true;
if (matchchain_cascade)
    warning('TODO: check correctness of cascaded deleting!');
    matchchain = handles.data.markersets(handles.data.sel_markerset).findmatchesparam.matchchain;
    for (i=sel_markers)
        remidx = sel_images;
        matchchaini = matchchain(1,:,i);
        while (true)
            matchchaini(remidx) = nan;
            remidx = find(ismember(matchchaini, remidx));
            if (isempty(remidx))
                break;
            end;
            setzeromask(1, remidx, i) = true;
        end;
    end;
end;
handles.data.markersets(handles.data.sel_markerset).markerset(1:2, setzeromask) = nan;
handles.data.markersets(handles.data.sel_markerset).findmatchesparam.matchchain(setzeromask) = nan;


handles.data.fileinfo(sel_images) = getFileInfo_info(handles, sel_images);
updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_markers_new(hObject, eventdata, handles)


if (isempty(handles.data.markersets))
    return;
end;


answ = inputdlg('How many new markers do you want to create?', 'New marker', 1, {num2str(1)});

if (length(answ)~=1)
    return;
end;


n = str2double(answ{1});
if (~isfinite(n) || round(n)~=n || n<1)
    return;
end;





sel_markerset = handles.data.sel_markerset;


msize = size3(handles.data.markersets(sel_markerset).markerset);

handles.data.markersets(sel_markerset).markerset(1:2,1:msize(2),(msize(3)+1):(msize(3)+n)) = nan;

handles.data.fileinfo = getFileInfo_info(handles, 1:length(handles.data.filenames));

sel_images = getListboxValue(handles.listbox_images);
sel_markers = [(msize(3)+1):(msize(3)+n), sort(getListboxValue(handles.listbox_markers))];

updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_markers_delete(hObject, eventdata, handles)

if (isempty(handles.data.markersets))
    return;
end;

sel_markerset = handles.data.sel_markerset;

sel_markers = sort(getListboxValue(handles.listbox_markers));

msize = size3(handles.data.markersets(sel_markerset).markerset);

if (isempty(sel_markers))
    answ = questdlg(['Are you sure, that you want to delete all ' num2str(msize(3)) ' markers from the markerset?'], 'Delete markers', 'Yes', 'No', 'Cancel', 'Yes');
    sel_markers = 1:msize(3);
else
    answ = questdlg(['Are you sure, that you want to delete the selected ' num2str(length(sel_markers)) ' markers?'], 'Delete markers', 'Yes', 'No', 'Cancel', 'Yes');
end;


if (~strcmp(answ, 'Yes'))
    return;
end;
index = sort(setdiff(1:msize(3), sel_markers));

handles.data.markersets(sel_markerset).markerset = handles.data.markersets(sel_markerset).markerset(1:2, :, index);
handles.data.markersets(sel_markerset).findmatchesparam.matchchain = handles.data.markersets(sel_markerset).findmatchesparam.matchchain(1, :, index);
handles.data.markersets(sel_markerset).findmatchesparam.comparelist(ismember(handles.data.markersets(sel_markerset).findmatchesparam.comparelist(:,1), sel_markers), :) = [];


handles.data.fileinfo = getFileInfo_info(handles, 1:length(handles.data.filenames));

sel_images = getListboxValue(handles.listbox_images);

updateListboxDisplay(handles);
updateSelection(handles, sel_images, []);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deletes the selected markerset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_delete(hObject, eventdata, handles)

if (isempty(handles.data.markersets))
    warndlg('No markerset available', 'Delete markerset', 'modal');    
    return;
end;

sel_markerset = handles.data.sel_markerset;

msize = size(handles.data.markersets(sel_markerset).markerset);
if (length(msize)==2)
    msize(3) = 1;
end;

if (msize(3) > 0) 
    answ = questdlg(['Are you sure, that you want to delete the ' num2str(sel_markerset) '. markerset ("' handles.data.markersets(sel_markerset).name '") ' ...
                     '(' num2str(msize(3)) ' markers and ' num2str(sum(sum(all(handles.data.markersets(sel_markerset).markerset > 0, 1)))) ' known positions)?'], 'Delete markerset', 'Yes', 'No', 'Cancel', 'Yes');

    if (~strcmp(answ, 'Yes'))
        return;
    end;
end;

handles.data.markersets(sel_markerset) = [];
handles.data.sel_markerset = max(1, handles.data.sel_markerset - 1);
if (isempty(handles.data.markersets))
    sel_markers = [];
else
    sel_markers = 1:size(handles.data.markersets(handles.data.sel_markerset).markerset, 3);
end;

handles.data.fileinfo = getFileInfo_info(handles, 1:length(handles.data.filenames));

set(handles.mainmenu_markerset, 'Enable', getOnOff(~isempty(handles.data.filenames)));


sel_images = getListboxValue(handles.listbox_images);
updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = toggleState(state)

if (strcmp(state,'on'))
    state = 'off';
else
    state = 'on';
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a new and empty markerset and selectes it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_new(hObject, eventdata, handles)


if (isempty(handles.data.filenames))
    errordlg('Insert first the images to the project!', 'Loading Markerfile', 'modal');
    return;
end;

handles.data.sel_markerset = length(handles.data.markersets) + 1;
newms = getDefaultMarkerset({ handles.data.markersets(:).name }, length(handles.data.filenames));
newmsfn = fieldnames(newms);
for (i=1:length(newmsfn))
    handles.data.markersets(handles.data.sel_markerset).(newmsfn{i}) = newms.(newmsfn{i});
end;


sel_markers = 1:size(handles.data.markersets(handles.data.sel_markerset).markerset, 3);

handles.data.fileinfo = getFileInfo_info(handles, 1:length(handles.data.filenames));

sel_images = getListboxValue(handles.listbox_images);
updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_infomarker_selectmarker(hObject, eventdata, handles, step) %#ok<INUSL,DEFNU>



sel_markers = getListboxValue(handles.listbox_markers);
lsel_markers = length(sel_markers);
[sel_markers_sort, sel_markers_sortidx] = sort(sel_markers);

switch (step)
    case {-1, +1} % Step prev
        if (lsel_markers <= 0)
            return;
        end;
        sel_marker = sel_markers(1);        
        idx = mod(find(sel_marker==sel_markers_sort) + step -1, lsel_markers) + 1;
        new_sel_marker = sel_markers_sort(idx);
        
        new_sel_markers = [sel_markers_sort(idx) sel_markers];
        new_sel_markers(sel_markers_sortidx(idx)+1) = [];
        
    case 0
        s = get(handles.edit_infomarker_selnumber, 'String');
        if (isempty(s) || isempty(handles.data.markersets))
            return;
        end;
        msize = size(handles.data.markersets(handles.data.sel_markerset).markerset);
        if (length(msize)==2)
            msize(3) = 1;
        end;
        
        new_sel_markers = getIndexRangeFromString(s, msize(3));
        
        if (~get(handles.checkbox_infomarker_selnumber, 'Value'))
            new_sel_markers = [new_sel_markers sel_markers];
        end;
        [idx, new_sel_markersidx] = unique(new_sel_markers,'first');
        new_sel_markers = new_sel_markers(sort(new_sel_markersidx));        
    otherwise
        error('UNEXPECTED!');
end;

updateSelection(handles, getListboxValue(handles.listbox_images), new_sel_markers);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residx = getIndexRangeFromString(s, n)


smask_delimiter = s < '0' | s > '9';
smask_range = [s == '-', false];
smask_delimiteidx = [0, find(smask_delimiter), length(s)+1];
idx = [];
range = [];
addrange = false;
for (i=1:length(smask_delimiteidx)-1)
    if (smask_delimiteidx(i)+1 < smask_delimiteidx(i+1))
        idx(end+1) = str2double(s((smask_delimiteidx(i)+1):smask_delimiteidx(i+1)-1));
        range(end+1) = addrange;
        addrange = false;
    end;
    if (smask_range(smask_delimiteidx(i+1)))
        addrange = true;
    end;
end;
range(end+1) = addrange;
residx = [];
for (i=1:length(idx))
    if (range(i) || range(i+1))
        if (range(i) && i==1)
            residx = [residx, 1:idx(1)];
        end;
        if (range(i+1))
            if (i==length(idx))
                residx = [residx, idx(i):n];
            else
                residx = [residx, idx(i):idx(i+1)];
            end;                    
        end;
    else
        residx(end+1) = idx(i);
    end;
end;

residx(residx>n | residx<1) = [];

[residxsort, i] = unique(residx, 'first');
residx = residx(sort(i));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_markers_select(hObject, eventdata, handles)

if (isempty(handles.data.markersets) || isempty(handles.data.markersets(handles.data.sel_markerset).markerset))
    return;
end;

[cfg, sel_markers] = tom_markGui_SelectMarker();

if (isempty(cfg))
    return;
end;

handles.data.markGui_SelectMarker = cfg;

sel_images = getListboxValue(handles.listbox_images);

updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_markers_selinvert(hObject, eventdata, handles)

if (isempty(handles.data.markersets) || isempty(handles.data.markersets(handles.data.sel_markerset).markerset))
    return;
end;

sel_markers = getListboxValue(handles.listbox_markers);
sel_images = getListboxValue(handles.listbox_images);

sel_markers_new = setdiff(1:size(handles.data.markersets(handles.data.sel_markerset).markerset, 3), sel_markers);

updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers_new);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_merge(hObject, eventdata, handles)

if (isempty(handles.data.markersets))
    warndlg('No markerset available', 'Clone/Merge markerset', 'modal');    
    return;
end;


if (length(handles.data.markersets) == 1)
    selected = 1;
else
    s = get(handles.listbox_markersets, 'String');
    
    [selected, answ] = listdlg('PromptString', 'Select the markerset(s) to clone/merge', 'SelectionMode', 'multiple', 'ListString' , s, 'InitialValue', handles.data.sel_markerset, 'Name', 'Clone/Merge markerset', 'ListSize', [300, 200]);
    
    if (~answ || isempty(selected))
        return;
    end;    
end;

sel_markerset = length(handles.data.markersets) + 1;
if (length(selected) == 1)
    handles.data.markersets(sel_markerset) = handles.data.markersets(selected);
else
    
    handles.data.markersets(sel_markerset).name = ['merge <' strrep(num2str(selected), '  ', ',') '>'];
    handles.data.markersets(sel_markerset).filename = 'markerfile.em';
    handles.data.markersets(sel_markerset).markerset = cat(3, handles.data.markersets(selected).markerset);
    handles.data.markersets(sel_markerset).findmatchesparam.matchchain = nan(1, size(handles.data.markersets(sel_markerset).markerset, 2), size(handles.data.markersets(sel_markerset).markerset, 3));

    comparelist = zeros(0, 3);
    baseidx = 0;
    for (i=selected)
        msize = size3(handles.data.markersets(i).markerset);
        
        handles.data.markersets(sel_markerset).findmatchesparam.matchchain(1,1:msize(2), baseidx + (1:msize(3))) = handles.data.markersets(i).findmatchesparam.matchchain;

        if (~isempty(handles.data.markersets(i).findmatchesparam.comparelist))
            comparelist = cat(1, comparelist, handles.data.markersets(i).findmatchesparam.comparelist + baseidx);
        end;
        baseidx = baseidx + msize(3);
    end;
    handles.data.markersets(sel_markerset).findmatchesparam.comparelist = comparelist;
end;

sel_images = getListboxValue(handles.listbox_images);
sel_markers = getListboxValue(handles.listbox_markers);

updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_findparameters(hObject, eventdata, handles) 

findmatchesparam = handles.data.findmatchesparam;

if (isempty(handles.data.markersets))
    msset = cell(0, 2);
else
    findmatchesparam.comparelist = handles.data.markersets(handles.data.sel_markerset).findmatchesparam.comparelist;
    findmatchesparam.matchchain = handles.data.markersets(handles.data.sel_markerset).findmatchesparam.matchchain;
    msset = cat(2, get(handles.listbox_markersets, 'String'), { handles.data.markersets(:).markerset }');
end;

[findmatchesparam, answ] = tom_markFindMatchesParam(length(handles.data.filenames), findmatchesparam, 'markersets', msset);
 
if (~answ)
    return;
end;

handles.data.findmatchesparam = rmfield(findmatchesparam, {'comparelist','matchchain'});
guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_findstart(hObject, eventdata, handles) 


if (isempty(handles.data.markersets))
    warndlg('No markerset available.', 'Fill markerset by correlation of neighbours', 'modal');
    return;
end;

sel_markerset = handles.data.sel_markerset;

markerset = handles.data.markersets(sel_markerset);

findmatchesparam = handles.data.findmatchesparam;
findmatchesparam.comparelist = markerset.findmatchesparam.comparelist;
findmatchesparam.matchchain = markerset.findmatchesparam.matchchain;

findmatchesparam_default = tom_markFindMatchesParamHelperFcn('getDefaultConfig', length(handles.data.filenames));

findmatchesparam = tom_markFindMatchesParamHelperFcn('parseConfig', length(handles.data.filenames), findmatchesparam, findmatchesparam_default);


delete(findall(0, 'type', 'figure', 'Tag', 'TMWWaitbar_tom_markGui_findmatches'));
hwaitbar = waitbar(0, 'please wait... filling markerset', 'CreateCancelBtn', 'set(gcf, ''Userdata'', true);', 'Userdata', false, 'Name', 'markGui: finding matches', 'Tag', 'TMWWaitbar_tom_markGui_findmatches');
%set(waitbarh, 'WindowStyle', 'modal');
set(hwaitbar, 'WindowStyle', 'normal');

if (findmatchesparam.comparelist_use)
    comparelist = findmatchesparam.comparelist;
else
    comparelist = [];
end;

tic;
[fmarkerset, fmatchchain, fcomparelist] = ...
    tom_mark_findmatches(   handles.data.filenames, markerset.markerset, findmatchesparam.ncorrelate, comparelist, ...
                            findmatchesparam, findmatchesparam.imreadbinning, findmatchesparam.filter, findmatchesparam.verbose, hwaitbar);
toc;                        

close(hwaitbar);
delete(hwaitbar);

handles.data.sel_markerset = length(handles.data.markersets) + 1;
handles.data.markersets(handles.data.sel_markerset).name = handles.data.markersets(sel_markerset).name;
handles.data.markersets(handles.data.sel_markerset).filename = handles.data.markersets(sel_markerset).filename;
handles.data.markersets(handles.data.sel_markerset).markerset = fmarkerset;
handles.data.markersets(handles.data.sel_markerset).findmatchesparam.comparelist = fcomparelist;

idx = ~(fmatchchain>0);
idx(~(findmatchesparam.matchchain>0)) = false;

try
fmatchchain(idx) = findmatchesparam.matchchain(idx);
handles.data.markersets(handles.data.sel_markerset).findmatchesparam.matchchain = fmatchchain;
catch
    warning('---');
end
sel_markers = 1:size(handles.data.markersets(handles.data.sel_markerset).markerset, 3);


sel_images = getListboxValue(handles.listbox_images);
updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_findsplit(hObject, eventdata, handles) 

if (isempty(handles.data.markersets))
    warndlg('No markerset available.', 'Split markerset', 'modal');
    return;
end;

sel_markerset = handles.data.sel_markerset;
markerset = handles.data.markersets(sel_markerset).markerset;

[tiltangle, imsize] = getTiltangles(handles.data.emheader);


while (true)
    threshold = ceil(0.09 * mean(imsize(:)));

    threshold = inputdlg('Enter threshold in pixels. This is the maximal value that a reprojected point may differ before splitting the marker', ...
        'Split markerset using RANSAC', 1, {num2str(threshold)});

    if (numel(threshold)~=1 || isempty(threshold{1})) 
        return; 
    end

    threshold = str2double(threshold{1});
    if (isfinite(threshold) && threshold > 0)
        break;
    end;
    uiwait(errordlg('Enter a positive distance in pixels', 'Split markerset', 'modal'));
end;


delete(findall(0, 'type', 'figure', 'Tag', 'TMWWaitbar_tom_markGui_splitset'));
hwaitbar = waitbar(0, 'please wait... splitting markerset', 'CreateCancelBtn', 'set(gcf, ''Userdata'', true);', 'Userdata', false, 'Name', 'markGui: splitting markerset', 'Tag', 'TMWWaitbar_tom_markGui_splitset');
%set(waitbarh, 'WindowStyle', 'modal');
set(hwaitbar, 'WindowStyle', 'normal');

markerset = tom_mark_cvaf_split_msset(markerset, threshold, 1000, hwaitbar);

if (isempty(markerset) || (ishandle(hwaitbar) && get(hwaitbar, 'UserData')))
    markerset = [];
end;

close(hwaitbar);
delete(hwaitbar);


if (isempty(markerset))
    return;
end;


handles.data.sel_markerset = length(handles.data.markersets) + 1;
handles.data.markersets(handles.data.sel_markerset).name = [handles.data.markersets(sel_markerset).name '/' num2str(threshold)];
handles.data.markersets(handles.data.sel_markerset).filename = 'markerset.em';
handles.data.markersets(handles.data.sel_markerset).markerset = markerset;
handles.data.markersets(handles.data.sel_markerset).findmatchesparam = struct('comparelist', {zeros(0, 3)}, 'matchchain', {nan(1, size(markerset,2), size(markerset,3))});


sel_markers = 1:size(handles.data.markersets(handles.data.sel_markerset).markerset, 3);


sel_images = getListboxValue(handles.listbox_images);
updateListboxDisplay(handles);
updateSelection(handles, sel_images, sel_markers);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function answ = view_correspondence(filenames, markerset, iim, idxinliers, distance, P)

idxex = squeeze(all(all(markerset(1:2,:,:) > 0, 1), 2));
markerset = markerset(1:2,:,idxex);
idxinliers = find(ismember(find(idxex), idxinliers));
idxoutliers = setdiff(1:size(markerset,3), idxinliers);

h = findobj('Type', 'figure', 'Tag', 'tom_markGui_findsplit');
if (length(h) ~= 1)
    delete(h);
    h = figure('Tag', 'tom_markGui_findsplit');
else
    figure(h);                
end;
set(h, 'UserData', 'abbort', 'Visible', 'on');

hp = findobj('Style', 'pushbutton', 'Parent', h, 'Tag', 'pushbutton_tom_markGui_findsplit_prev');
if (length(hp) ~= 1)
    delete(hp);
    hp = uicontrol('Parent', h, 'Style', 'pushbutton', ...
                    'Callback', 'set(gcf,''userdata'',''prev'');uiresume();', ...
                    'String', '<<', 'Tag', 'pushbutton_tom_markGui_findsplit_prev'); 
end;

oldhp = hp;
hp = findobj('Style', 'pushbutton', 'Parent', h, 'Tag', 'pushbutton_tom_markGui_findsplit_next');
if (length(hp) ~= 1)
    delete(hp);
    hp = uicontrol('Parent', h, 'Style', 'pushbutton', ...
                    'Callback', 'set(gcf,''userdata'',''next'');uiresume();', 'Position', get(oldhp,'Position')+[65 0 0 0], ...
                    'String', '>>', 'Tag', 'pushbutton_tom_markGui_findsplit_next'); 
end;

oldhp = hp;
hp = findobj('Style', 'pushbutton', 'Parent', h, 'Tag', 'pushbutton_tom_markGui_findsplit_step');
if (length(hp) ~= 1)
    delete(hp);
    hp = uicontrol('Parent', h, 'Style', 'pushbutton', ...
                    'Callback', 'set(gcf,''userdata'',''step'');uiresume();', 'Position', get(oldhp,'Position')+[65 0 0 0], ...
                    'String', 'next', 'Tag', 'pushbutton_tom_markGui_findsplit_step'); 
end;

oldhp = hp;
hp = findobj('Style', 'pushbutton', 'Parent', h, 'Tag', 'pushbutton_tom_markGui_findsplit_abbort');
if (length(hp) ~= 1)
    delete(hp);
    hp = uicontrol('Parent', h, 'Style', 'pushbutton', ...
                    'Callback', 'set(gcf,''userdata'',''abbort'');uiresume();', 'Position', get(oldhp,'Position')+[65 0 0 0], ...
                    'String', 'cancel', 'Tag', 'pushbutton_tom_markGui_findsplit_abbort'); 
end;



answ = false;

iiim = 1;
while (~answ)
    
    set(h, 'Name', [num2str(iim(iiim)) ': ' num2str(iim) ': ' num2str(length(idxinliers)) ' of ' num2str(size(markerset, 3))]);            
    

    im = tom_emread(filenames{iiim});


    cla; imagesc(im.Value'); colormap gray; axis equal on; hold on;
    
    plot(reshape(markerset(1, :, idxinliers), [length(iim), length(idxinliers)]), reshape(markerset(2, :, idxinliers), [length(iim), length(idxinliers)]), 'bx-');
    plot(squeeze(markerset(1, iiim, idxinliers)), squeeze(markerset(2, iiim, idxinliers)), 'bo');
    
    plot(reshape(markerset(1, :, idxoutliers), [length(iim), length(idxoutliers)]), reshape(markerset(2, :, idxoutliers), [length(iim), length(idxoutliers)]), 'rx-');
    plot(squeeze(markerset(1, iiim, idxoutliers)), squeeze(markerset(2, iiim, idxoutliers)), 'ro');
    

    uiwait(h);
    if (~ishandle(h) || strcmp(get(h, 'Userdata'), 'abbort'))
        if (ishandle(h))
            set(h, 'Visible', 'off');
        end;
        return;
    end;
    switch (get(h, 'Userdata'))
        case 'next'
            iiim = mod(iiim, size(markerset,2))+1;
        case 'prev'
            iiim = mod(iiim-2, size(markerset,2))+1;
        otherwise
            answ = true;
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% l being a Nx3 matrix with homogen lines
% r being a window [x1, x2; y1, y2]
% x a 3x2xN matrix of homogen endpoints with last coordinate 1
% nan if outside of window.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = getLineEndPoints(l, rrect)

lrect(1, 1:3) = cross([rrect(1,1), rrect(2,1), 1]', [rrect(1,2), rrect(2,1), 1]')'; %top
lrect(2, 1:3) = cross([rrect(1,1), rrect(2,2), 1]', [rrect(1,2), rrect(2,2), 1]')'; %bottom
lrect(3, 1:3) = cross([rrect(1,1), rrect(2,1), 1]', [rrect(1,1), rrect(2,2), 1]')'; %left
lrect(4, 1:3) = cross([rrect(1,2), rrect(2,1), 1]', [rrect(1,2), rrect(2,2), 1]')'; %right

x = nan(3, 2, size(l, 1));

for (i=1:size(l,1))
    for (j=1:4)
        prect(j,1:3) = cross(lrect(j, 1:3), l(i, 1:3));
        if (abs(prect(j,3)) > 1e3*eps)
            prect(j,:) = prect(j,:) ./ prect(j,3);
        else
            prect(j,:) = nan;
        end;            
    end;
    pp = zeros(3, 0);
    
    if (prect(1,1)>=rrect(1,1) && prect(1,1)<=rrect(1, 2))
        pp(1:3, end+1) = prect(1,:);
    end;
    if (prect(3,2)>=rrect(2,1) && prect(3,2)<=rrect(2, 2))
        pp(1:3, end+1) = prect(3,:);
    end;
    if (prect(2,1)>=rrect(1,1) && prect(2,1)<=rrect(1, 2))
        pp(1:3, end+1) = prect(2,:);
    end;
    if (prect(4,2)>=rrect(2,1) && prect(4,2)<=rrect(2, 2))
        pp(1:3, end+1) = prect(4,:);
    end;
    if (~isempty(pp))
        if (size(pp, 2) == 2)
            x(1:3,1:2,i) = pp;
        else
            warning('todo: check corner behavior...');
        end;            
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx = getRecIdx(baseset, idx)

mask = ismember(baseset, idx);
fmask = find(mask);
if (isempty(fmask))
    return;
end;
baseset(mask) = nan;
idx = [idx, getRecIdx(baseset, fmask)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Returns an array of tiltangles from the cellarray emheader.
% emheader has for every image, its emheader.
% Missing tiltangles are interpolated!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tiltangle, imsize] = getTiltangles(emheader, interpolate) 

ni = length(emheader);

% Get the tiltangles.
Tiltangle = nan(1, ni);
imsize = nan(2, ni);
for (i=1:ni)
    if (~isempty(emheader{i}))
        imsize(1:2, i) = emheader{i}.Size(1:2);
        Tiltangle(i) = emheader{i}.Tiltangle;
    end;    
end;

if (~exist('interpolate', 'var') || interpolate)
    % For missing images, interpolate the tiltangles and imagesizes.
    idxmissing_image = find(~isfinite(Tiltangle));
    if (~isempty(idxmissing_image))
        if (length(idxmissing_image) == length(emheader))
            imsize = ones(2,ni) * 1024;
            Tiltangle = 120 * (0:1/ni:1) - 60;
        else
            imsize(1:2, idxmissing_image) = repmat(mean(imsize(1:2, setdiff(1:ni, idxmissing_image)), 2), [1, length(idxmissing_image)]);
            Tiltangle(idxmissing_image) = interp1(setdiff(1:ni, idxmissing_image), Tiltangle(isfinite(Tiltangle)), idxmissing_image, 'spline');
        end;
    end;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_align_affineitr(hObject, eventdata, handles) 
if (isempty(handles.data.markersets))
    warndlg('No markerset available.', 'Alignment', 'modal');
    return;
end;

warndlg('Not yet implmented.', 'Alignment', 'modal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_mainmenu_markersets_align_rigbody(hObject, eventdata, handles) 
if (isempty(handles.data.markersets))
    warndlg('No markerset available.', 'Alignment', 'modal');
    return;
end;

sel_markerset = handles.data.sel_markerset;
markerset = handles.data.markersets(sel_markerset).markerset;

[tiltangles, imsize] = getTiltangles(handles.data.emheader); 



status = figureDisableAll(handles.figure_markGui);
figureEnableAll(status); status = cell(0, 3);

try
    switch (get(hObject, 'Tag'))
        case 'mainmenu_markersets_align_rigbodyitr'
            [answ, aligncfg, aligndata] = tom_markGui_Alignment3dItr();
        otherwise
            if (isfield(handles.data.markersets(sel_markerset), 'alignment3d') && isstruct(handles.data.markersets(sel_markerset).alignment3d) && numel(handles.data.markersets(sel_markerset).alignment3d)==1)
                aligncfg = handles.data.markersets(sel_markerset).alignment3d;
            else
                aligncfg.tiltangles = tiltangles;
                aligncfg.imsize = mean(imsize(:));
                sel_markers = getListboxValue(handles.listbox_markers);
                if (length(sel_markers) == 1)
                    aligncfg.refmark = sel_markers;
                end;
            end;
            [answ, aligncfg, aligndata] = tom_markGui_Alignment3d(markerset, handles.data.markersets(sel_markerset).name, { handles.data.fileinfo(:).filename }, aligncfg);
    end;
        
catch    
    figureEnableAll(status);
    rethrow(lasterror);
end;
figureEnableAll(status);

switch ([get(hObject, 'Tag') '_' answ])
    case {'mainmenu_markersets_align_rigbodyitr_cancel', ...
          'mainmenu_markersets_align_rigbody_cancel'}
    case {'mainmenu_markersets_align_rigbody_ok'}
        handles.data.markersets(sel_markerset).alignment3d = aligncfg;
        guidata(handles.figure_markGui, handles);
    case {'mainmenu_markersets_align_rigbodyitr_ok'}
        handles.data.markersets(sel_markerset).alignment3ditr = aligncfg;
        guidata(handles.figure_markGui, handles);
    case {'mainmenu_markersets_align_rigbody_okimport'} 
        handles.data.markersets(sel_markerset).alignment3d = aligncfg;
        handles.data.sel_markerset = length(handles.data.markersets)+1;
        handles.data.markersets(handles.data.sel_markerset).name = [handles.data.markersets(sel_markerset).name ' (rb)'];
        handles.data.markersets(handles.data.sel_markerset).filename = 'markerfile.em';
        handles.data.markersets(handles.data.sel_markerset).markerset = aligndata.x;
        handles.data.markersets(handles.data.sel_markerset).findmatchesparam.matchchain = nan(1, size(aligndata.x, 2), size(aligndata.x, 3));
        handles.data.markersets(handles.data.sel_markerset).findmatchesparam.comparelist = zeros(0, 3);

        sel_images = getListboxValue(handles.listbox_images);
        sel_markers = getListboxValue(handles.listbox_markers);

        updateListboxDisplay(handles);
        updateSelection(handles, sel_images, sel_markers);
    case {'mainmenu_markersets_align_rigbodyitr_okimport'}
        handles.data.markersets(sel_markerset).alignment3ditr = aligncfg;
        handles.data.sel_markerset = length(handles.data.markersets)+1;
        handles.data.markersets(handles.data.sel_markerset).name = [handles.data.markersets(sel_markerset).name ' (rbi)'];
        handles.data.markersets(handles.data.sel_markerset).filename = 'markerfile.em';
        handles.data.markersets(handles.data.sel_markerset).markerset = aligndata.x;
        handles.data.markersets(handles.data.sel_markerset).findmatchesparam.matchchain = nan(1, size(aligndata.x, 2), size(aligndata.x, 3));
        handles.data.markersets(handles.data.sel_markerset).findmatchesparam.comparelist = zeros(0, 3);

        sel_images = getListboxValue(handles.listbox_images);
        sel_markers = getListboxValue(handles.listbox_markers);

        updateListboxDisplay(handles);
        updateSelection(handles, sel_images, sel_markers);

    otherwise
        warning('Unexpected answer!');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = figureDisableAll(hObject)

hChildren = findall(hObject);
nc = length(hChildren);
status = cell(0, 3);


for (i=1:nc)
    if (ishandle(hChildren(i)))
        hGet = get(hChildren(i));
        switch (hGet.Type)
            case {'uicontrol'}
                status{end+1, 1} = hChildren(i);
                status{end, 2} = 'Enable';
                status{end, 3} = hGet.Enable;
                set(hChildren(i), 'Enable', 'inactive');
                status{end+1, 1} = hChildren(i);
                status{end, 2} = 'ButtonDownFcn'; 
                status{end, 3} = hGet.ButtonDownFcn;
                set(hChildren(i), 'ButtonDownFcn', '');
            case {'uimenu'}
                status{end+1, 1} = hChildren(i);
                status{end, 2} = 'Enable';
                status{end, 3} = hGet.Enable;
                set(hChildren(i), 'Enable', 'off');
            case {'text', 'line', 'image', 'uipanel'}
                status{end+1, 1} = hChildren(i);
                status{end, 2} = 'ButtonDownFcn';
                status{end, 3} = hGet.ButtonDownFcn;
                set(hChildren(i), 'ButtonDownFcn', '');
        end;
    end;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figureEnableAll(status)

nc = size(status, 1);
for (i=1:nc)
    if (ishandle(status{i, 1}))
        set(status{i, 1}, status{i,2}, status{i, 3});
    end;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xxx



figure(1);
axis equal;

ms_diff = (ms - msr);
dist = sqrt(sum(ms_diff.^2,1));

plot(reshape(ms_diff(1,:,:), [1,numel(ms_diff)/2]), reshape(ms_diff(2,:,:), [1,numel(ms_diff)/2]), 'b.');
rectangle('Curvature', [1 1], 'Position', max(dist(:)) * [-0.5 -0.5 1 1]*2);
plot(sort(reshape(dist,[1,numel(dist)])));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Movie player
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlayMovie_Callback(hObject, eventdata, handles) 

% Get data
data = handles.data;

% Get projections
projnames = data.filenames;

% Get markerfile
mf = data.markersets(1,data.sel_markerset).markerset;

% Get feature selection
T1val = get(handles.listbox_markers, 'Value');
featureindex = T1val(1,1);

% Get movie size
moviesize = str2num(get(handles.editmoviesize, 'String'));

% Play movie
tom_mark_movie(projnames, mf, featureindex, str2num(get(handles.editmoviesize, 'String')), 0, ...
                                                                              str2num(get(handles.editdelay, 'String')));
% -------------------------------------------------------------------------




