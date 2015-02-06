function varargout = tom_markFindMatchesParam(varargin)
% TOM_MARKFINDMATCHESPARAM M-file for tom_markFindMatchesParam.fig
%      TOM_MARKFINDMATCHESPARAM, by itself, creates a new TOM_MARKFINDMATCHESPARAM or raises the existing
%      singleton*.
%
%      H = TOM_MARKFINDMATCHESPARAM returns the handle to a new TOM_MARKFINDMATCHESPARAM or the handle to
%      the existing singleton*.
%
%      TOM_MARKFINDMATCHESPARAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_MARKFINDMATCHESPARAM.M with the given input arguments.
%
%      TOM_MARKFINDMATCHESPARAM('Property','Value',...) creates a new TOM_MARKFINDMATCHESPARAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_markFindMatchesParam_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_markFindMatchesParam_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_markFindMatchesParam

% Last Modified by GUIDE v2.5 10-Jun-2007 14:00:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_markFindMatchesParam_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_markFindMatchesParam_OutputFcn, ...
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tom_markFindMatchesParam_OpeningFcn(hObject, eventdata, handles, varargin)

% Check and read inputparameters.


% First parameter: number of images.
if (length(varargin) < 1 || ~isnumeric(varargin{1}) || numel(varargin{1})~=1 || round(varargin{1})~=varargin{1} || varargin{1}<2)
    error('Needs number of images as first parameter.');
end;
handles.data.nimages = double(varargin{1});

% Second paramter: stucture with the configuration.
if (length(varargin)<2 || ~(isstruct(varargin{2})||isempty(varargin{2})))
    error('Parameter 2 must be a structure containing the configuration fields.');
end;
handles.data.config0 = varargin{2};

i=3;
handles.data.markersets = cell(0,2);
parse_only = false;

while (i<=length(varargin))
    if (ischar(varargin{i}))
        switch(varargin{i})
            case 'markersets'
                if (length(varargin) < i+1)
                    error('Parameter markersets needs argument');
                else
                    if (~iscell(varargin{i+1}))
                        error('Argument of "markersets" must be a Nx2 cell array with the names and the positions of existing markersets.');
                    end;
                    if (isempty(varargin{i+1}))
                        handles.data.markersets = cell(0, 2);
                    else
                        handles.data.markersets = varargin{i+1};
                        if (size(handles.data.markersets, 2) ~= 2)
                            error('Argument of "markersets" must be a Nx2 cell array with the names and the positions of existing markersets.');
                        end;
                        for (ii=1:size(handles.data.markersets, 1))
                            if (~ischar(handles.data.markersets{ii, 1}) || ~isnumeric(handles.data.markersets{ii, 2}))
                                error('Argument of "markersets" must be a Nx2 cell array with the names and the positions of existing markersets.');
                            end;
                            if (~isempty(handles.data.markersets{ii, 2}) && size(handles.data.markersets{ii,2},2)~=handles.data.nimages)
                                error('The markersets don''t have the same number of images.');
                            end;
                        end;
                    end;
                end;
                i = i+1;
            case 'parse_only'
                parse_only = true;
            otherwise 
                error('Unknown argument');
        end;                    
    else
        error('Unknown argument');
    end;
    i = i+1;
end;


handles.data.configdefault = tom_markFindMatchesParamHelperFcn('getDefaultConfig', handles.data.nimages);


handles.data.config = tom_markFindMatchesParamHelperFcn('parseConfig', handles.data.nimages, handles.data.config0, handles.data.configdefault);


if (parse_only)
    handles.data.closecommand = 'ok';
    guidata(hObject, handles);
else
    
    loadGui(handles);


    set(handles.figure_markFindMatchesParam, 'CloseRequestFcn', '');

    guidata(hObject, handles);

    uiwait(handles.figure_markFindMatchesParam);
end;


    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_markFindMatchesParam_OutputFcn(hObject, eventdata, handles) 


if (~isfield(handles.data, 'closecommand') || ~strcmp(handles.data.closecommand, 'ok'))
    varargout{1} = handles.data.config0;
    varargout{2} = false;
else
    handles = unloadGui(handles);
    varargout{1} = handles.data.config;
    varargout{2} = true;
end;

delete(handles.figure_markFindMatchesParam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = getOnOff(status)

if (status)
    s = 'on';
else
    s = 'off';
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imshift = getImShift(handles)

if (get(handles.radiobutton_im_shift_none, 'Value'))
    imshift = tom_mark_getAvgShift(handles.data.nimages);
elseif (get(handles.radiobutton_im_shift_set, 'Value'))
    imshift = handles.data.config.im_shift;
elseif (get(handles.radiobutton_im_shift_zero, 'Value'))
    imshift = zeros(handles.data.nimages, handles.data.nimages, 2);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = unloadGui(handles)

nval = str2num(get(handles.edit_ncorrelate, 'String'));
if (numel(nval)~=1 || round(nval)~=nval || ~(nval >= 1))
    handles.data.config.ncorrelate = handles.data.configdefault.ncorrelate;
else
    handles.data.config.ncorrelate = nval;
end;

nval = str2num(get(handles.edit_maxdisparity, 'String'));
if (numel(nval)~=1 || ~(nval >= 1))
    handles.data.config.max_disparity = handles.data.configdefault.max_disparity;
else
    handles.data.config.max_disparity = nval;
end;

nval = str2num(get(handles.edit_min_val, 'String'));
if (numel(nval)~=1 || isnan(nval))
    handles.data.config.min_val = handles.data.configdefault.min_val;
else
    handles.data.config.min_val = nval;
end;

nval = str2num(get(handles.edit_w, 'String'));
if (numel(nval)~=1 || ~(nval>0))
    handles.data.config.w = handles.data.configdefault.w;
else
    handles.data.config.w = nval;
end;

nval = str2num(get(handles.edit_imreadbinning, 'String'));
if (numel(nval)~=1 || round(nval)~=nval || ~(nval>=0))
    handles.data.config.imreadbinning = handles.data.configdefault.imreadbinning;
else
    handles.data.config.imreadbinning = nval;
end;


handles.data.config.comparelist_use = get(handles.checkbox_comparelist, 'Value');
handles.data.config.mutual = get(handles.checkbox_mutual, 'Value');


if (get(handles.radiobutton_im_shift_zero, 'Value'))
    handles.data.config.im_shift_mode = 'zero';
elseif (get(handles.radiobutton_im_shift_set, 'Value'))
    handles.data.config.im_shift_mode = 'set';
else %if (get(handles.radiobutton_im_shift_none, 'Value'))
    handles.data.config.im_shift_mode = 'none';
end;
    
handles.data.config.im_shift = getImShift(handles);
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadGui(handles)


set(handles.edit_ncorrelate, 'String', num2str(handles.data.config.ncorrelate));
set(handles.edit_maxdisparity, 'String', num2str(handles.data.config.max_disparity));
set(handles.edit_min_val, 'String', num2str(handles.data.config.min_val));
set(handles.edit_w, 'String', num2str(handles.data.config.w));
set(handles.edit_imreadbinning, 'String', num2str(handles.data.config.imreadbinning));

set(handles.checkbox_comparelist, 'Value', 1, 'Enable', getOnOff(~isempty(handles.data.config.comparelist)));
set(handles.checkbox_mutual, 'Value', handles.data.config.mutual);

set(handles.(['radiobutton_im_shift_' handles.data.config.im_shift_mode]), 'Value', 1);
set(handles.listbox_im_shift_markersets, 'String', handles.data.markersets(:, 1)', 'Value', 1:size(handles.data.markersets,1));

Callback_radiobutton_im_shift(handles.(['radiobutton_im_shift_' handles.data.config.im_shift_mode]), [], handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_im_shift_compute(hObject, eventdata, handles) 

ms = [];
lvalue = get(handles.listbox_im_shift_markersets, 'Value');
for (i=lvalue)
    if (~isempty(handles.data.markersets{ i, 2 }))
        ms = cat(3, ms, handles.data.markersets{ i, 2 });
    end;
end;


handles.data.config.im_shift = tom_mark_getAvgShift(ms);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CloseRequestFcn_figure(hObject, eventdata, handles) 


closereq;
%delete(handles.figure_markFindMatchesParam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_im_shift_tows(hObject, eventdata, handles) 

im_shift = getImShift(handles);

assignin('base', 'im_shift_FindMatchesParam', im_shift);

evalin('base', 'open(''im_shift_FindMatchesParam'');');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_im_shift_fromws(hObject, eventdata, handles) 


vwhos = evalin('base', 'whos();');
vnames = {};
for (i=1:length(vwhos))
    if (any(strcmp(vwhos(i).class, {'double', 'single', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'})) && ...
        length(vwhos(i).size) == 3 && all(vwhos(i).size==[handles.data.nimages*[1, 1], 2]))
        vnames{end+1} = vwhos(i).name;        
    end;
end;

if (isempty(vnames))
    errordlg(['No suitable array of size ' num2str(handles.data.nimages) 'x' num2str(handles.data.nimages) 'x2 in base workspace found'], 'Load imshifts from workspace', 'modal');
    return;
end;

[selection, answ] = listdlg('PromptString', 'Select the workspace variable', 'SelectionMode', 'single', 'ListString' , vnames, 'InitialValue', 1, 'Name', 'Import from workspace', 'ListSize', [300, 200]);
if (~answ)
    return;
end;

handles.data.config.im_shift = double(evalin('base', vnames{selection}));

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_ok(hObject, eventdata, handles) 

handles.data.closecommand = 'ok';
guidata(hObject, handles);
uiresume;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_cancel(hObject, eventdata, handles) 

handles.data.closecommand = 'cancel';
guidata(hObject, handles);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_reset(hObject, eventdata, handles) 


handles.data.config = tom_markFindMatchesParamHelperFcn('parseConfig', handles.data.nimages, handles.data.config0, handles.data.configdefault);

loadGui(handles);
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_default(hObject, eventdata, handles) 


handles.data.config = handles.data.configdefault;

loadGui(handles);
guidata(hObject, handles);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_im_shift_showdisp(hObject, eventdata, handles) 

im_shift = getImShift(handles);


ms = [];

lvalue = get(handles.listbox_im_shift_markersets, 'Value');
for (i=lvalue)
    if (~isempty(handles.data.markersets{ i, 2 }))
        ms = cat(3, ms, handles.data.markersets{ i, 2 });
    end;
end;

ncorrelate = str2num(get(handles.edit_ncorrelate, 'String'));
if (numel(ncorrelate)~=1 || isnan(ncorrelate))
    ncorrelate = handles.data.nimages;
end;
disparity = cell(1, ncorrelate);

maxlength = 0;
msize3 = size(ms, 3);
for (i=1:ncorrelate)
    dispi = [];
    for (j=1:(handles.data.nimages-i))
        dispj = squeeze(sqrt(sum((ms(:, [j],:) -ms(:,i+j,:) + repmat(squeeze(im_shift(j, i+j, 1:2)), [1 1 msize3])) .^ 2, 1)));
        dispj(isnan(dispj)) = [];
        dispi = [dispi; dispj];
    end;
    maxlength = max(maxlength, length(dispi));
    disparity{i} = dispi;
end;

h = findobj('Tag', 'figure_tom_markFindMatchesParam_showdisparity', 'Type', 'figure');
if (isempty(h))
    h = figure('Tag', 'figure_tom_markFindMatchesParam_showdisparity');
else
    figure(h);
end;


cla;
axis on;
hold on;
chsv = hsv(ncorrelate);
for (i=1:ncorrelate)
    j = length(disparity{i});
    if (j>0)
        plot((1:j)*maxlength/j, sort(disparity{i}), '.-', 'Color', chsv(i,:));
    end;
end;
hold off;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_radiobutton_im_shift(hObject, eventdata, handles) 

sTag = get(hObject, 'Tag');


set(handles.edit_maxdisparity, 'Enable', getOnOff(strcmp(sTag, 'radiobutton_im_shift_zero') || strcmp(sTag, 'radiobutton_im_shift_set')));
set(handles.pushbutton_im_shift_showdisp, 'Enable', getOnOff(strcmp(sTag, 'radiobutton_im_shift_zero') || strcmp(sTag, 'radiobutton_im_shift_set')));
set(handles.pushbutton_im_shift_compute, 'Enable', getOnOff(size(handles.data.markersets,1)>0 && strcmp(sTag, 'radiobutton_im_shift_set')));
set(handles.pushbutton_im_shift_fromws, 'Enable', getOnOff(strcmp(sTag, 'radiobutton_im_shift_set')));
%set(handles.pushbutton_im_shift_tows, 'Enable', getOnOff(status));
set(handles.listbox_im_shift_markersets, 'Enable', getOnOff(strcmp(sTag, 'radiobutton_im_shift_zero') || strcmp(sTag, 'radiobutton_im_shift_set')));






















