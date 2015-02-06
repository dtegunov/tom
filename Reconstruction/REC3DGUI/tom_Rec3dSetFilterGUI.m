function varargout = tom_Rec3dSetFilterGUI(varargin)
%TOM_REC3DSETFILTERGUI is a module of TOM_REC3DGUI.
%
%   FILTER = tom_Rec3dSetFilterGUI(FILTER)
%
%   It works on a structure which contains parameter for filtering.
%
%PARAMETERS
%
%  INPUT
%   FILTER              sub-structure of the main-structure of Rec3d
%                       FILTER = Rec3dProject.RECONSTRUCTION.FILTER
%  
%  OUTPUT
%   FILTER              aktualized sub-structure
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%   TOM_APPLY_FILTER   TOM_FILTER   TOM_BANDPASS
%
%   01/01/07 ME
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


% -------------------------------------------------------------------------
% Begin initialization code - DO NOT EDIT
% -------------------------------------------------------------------------
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_Rec3dSetFilterGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_Rec3dSetFilterGUI_OutputFcn, ...
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
% -------------------------------------------------------------------------
% End initialization code - DO NOT EDIT
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% OpeningFunction
% -------------------------------------------------------------------------
function tom_Rec3dSetFilterGUI_OpeningFcn(hObject, eventdata, handles, varargin)

if (size(varargin, 2) < 1)
    
    % Input FILTER
    handles.inputFILTER.Apply = 0;
    handles.inputFILTER.Value.times = 1;
    handles.inputFILTER.Value.low = 0;
    handles.inputFILTER.Value.high = 0;
    handles.inputFILTER.Value.smooth = 0;
    handles.inputFILTER.Value.space = 'real';
    handles.inputFILTER.Value.method = 'quadr';
    handles.inputFILTER.Value.radius = 0;
    handles.inputFILTER.Type = 'bandpass';
    
else
    
    % Input FILTER
    handles.inputFILTER = varargin{1};
    
end

% handles.inputFILTER -> handles.outputFILTER
handles.outputFILTER = handles.inputFILTER;

% Display FILTER
setRec3dSetFilterGUIValues(handles);

% Update handles structure
guidata(hObject, handles);
uiwait;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% OutputFunction
% -------------------------------------------------------------------------
function varargout = tom_Rec3dSetFilterGUI_OutputFcn(hObject, eventdata, handles)

% Output FILTER
varargout{1} = handles.outputFILTER;

% Delete Rec3dSetFilterGUI
delete(handles.Rec3dSetFilterGUI);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Cancel
% -------------------------------------------------------------------------
function pushbuttonCancel_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);
uiresume;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% OK
% -------------------------------------------------------------------------
function pushbuttonOK_Callback(hObject, eventdata, handles)

% Aktualize outputFILTER
handles = getRec3dSetFilterGUIValues(handles);

% Update handles structure
guidata(hObject, handles);
uiresume;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% setRec3dSetFilterGUIValues
% -------------------------------------------------------------------------
function setRec3dSetFilterGUIValues(handles);

% Simplify
% RECONSTRUCTION / Filter
Apply = handles.inputFILTER.Apply;
times = handles.inputFILTER.Value.times;
low = handles.inputFILTER.Value.low;
high = handles.inputFILTER.Value.high;
smooth = handles.inputFILTER.Value.smooth;
space = handles.inputFILTER.Value.space;
method = handles.inputFILTER.Value.method;
radius = handles.inputFILTER.Value.radius;
Type = handles.inputFILTER.Type;

% Type
if strcmp(Type, 'bandpass') == 1
    set(handles.popupType, 'Value', 1);
elseif strcmp(Type, 'kernel') == 1
    set(handles.popupType, 'Value', 2);
end

% Apply
if isequal(Apply,1)
    set(handles.popupApply, 'Value', 1); %on
elseif isequal(Apply,2)
    set(handles.popupApply, 'Value', 2); %default
elseif isequal(Apply,0)
    set(handles.popupApply, 'Value', 3); %off
end

% Times
set(handles.editTimes, 'String', num2str(times));

% Low
set(handles.editLow, 'String', num2str(low));

% High
set(handles.editHigh, 'String', num2str(high));

% Smooth
set(handles.editSmooth, 'String', num2str(smooth));

% Space
if strcmp(space, 'real') == 1
    set(handles.popupSpace, 'Value', 1);
elseif strcmp(space, 'fourier') == 1
    set(handles.popupSpace, 'Value', 2);
end

% Method
if strcmp(method, 'quadr') == 1
    set(handles.popupMethod, 'Value', 1);
elseif strcmp(method, 'circ') == 1
    set(handles.popupMethod, 'Value', 2);
end

% Radius
set(handles.editRadius, 'String', num2str(radius));
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% getRec3dSetFilterGUIValues
% -------------------------------------------------------------------------
function [handles] = getRec3dSetFilterGUIValues(handles)

% Apply
str = get(handles.popupApply, 'String');
val = get(handles.popupApply, 'Value');
ApplyNewWord = str{val};

if strcmp(ApplyNewWord, 'on')
    ApplyNew = 1;
elseif strcmp(ApplyNewWord, 'default')
    ApplyNew = 2;
elseif strcmp(ApplyNewWord, 'off')
    ApplyNew = 0;
end

% Times
timesNew = str2num(get(handles.editTimes, 'String'));

% Low
lowNew = str2num(get(handles.editLow, 'String'));

% High
highNew = str2num(get(handles.editHigh, 'String'));

% Smooth
smoothNew = str2num(get(handles.editSmooth, 'String'));

% Space
str = get(handles.popupSpace, 'String');
val = get(handles.popupSpace, 'Value');
spaceNew = str{val};

% Method
str = get(handles.popupMethod, 'String');
val = get(handles.popupMethod, 'Value');
methodNew = str{val};

% Radius
radiusNew = str2num(get(handles.editRadius, 'String'));

% Type
str = get(handles.popupType, 'String');
val = get(handles.popupType, 'Value');
TypeNew = str{val};

% Aktualize outputFILTER
handles.outputFILTER.Apply = ApplyNew;
handles.outputFILTER.Value.times = timesNew;
handles.outputFILTER.Value.low = lowNew;
handles.outputFILTER.Value.high = highNew;
handles.outputFILTER.Value.smooth = smoothNew;
handles.outputFILTER.Value.space = spaceNew;
handles.outputFILTER.Value.method = methodNew;
handles.outputFILTER.Value.radius = radiusNew;
handles.outputFILTER.Type = TypeNew;
% -------------------------------------------------------------------------

