function varargout = tom_Rec3dSetDetailGUI(varargin)
%TOM_REC3DSETDETAILGUI is a module of TOM_REC3DGUI.
%
%   Rec3dProject = tom_Rec3dSetDetailGUI(Rec3dProject,Rec3dReconstruction)
%
%   It works on the main-structure for detail-reconstruction.
%
%PARAMETERS
%
%  INPUT
%   Rec3dProject              main-structure of Rec3d
%   Rec3dReconstruction       overview reconstruction
%  
%  OUTPUT
%   Rec3dProject              aktualized main-structure of Rec3d
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%   TOM_VOLXYZ    
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
                   'gui_OpeningFcn', @tom_Rec3dSetDetailGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_Rec3dSetDetailGUI_OutputFcn, ...
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
function tom_Rec3dSetDetailGUI_OpeningFcn(hObject, eventdata, handles, varargin)

if (size(varargin, 2) < 1)
    
    % Input REC3DPROJECT
    handles.inputREC3DPROJECT = [];
    % Input REC3DRECONSTRUCTION
    handles.inputREC3DRECONSTRUCTION = [];
    
else
    
    % Input REC3DPROJECT
    handles.inputREC3DPROJECT = varargin{1};
    % Input REC3DRECONSTRUCTION
    handles.inputREC3DRECONSTRUCTION = varargin{2};
    
end

% handles.inputREC3DPROJECT -> handles.outputREC3DPROJECT 
handles.outputREC3DPROJECT = handles.inputREC3DPROJECT;

% Display REC3DPROJECT
setRec3dSetDetailGUIValues(handles);

% Update handles structure
guidata(hObject, handles);
uiwait;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% OutputFunction
% -------------------------------------------------------------------------
function varargout = tom_Rec3dSetDetailGUI_OutputFcn(hObject, eventdata, handles)

% Output METHOD
varargout{1} = handles.outputREC3DPROJECT;

% Delete Rec3dSetDetailGUI
delete(handles.Rec3dSetDetailGUI);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% pushbuttonCancel
% -------------------------------------------------------------------------
function pushbuttonCancel_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);
uiresume;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% pushbuttonOK
% -------------------------------------------------------------------------
function pushbuttonOK_Callback(hObject, eventdata, handles)

% Aktualize outputREC3DPROJECT
handles = getRec3dSetDetailGUIValues(handles);

% Update handles structure
guidata(hObject, handles);
uiresume;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% setRec3dSetDetailGUIValues
% -------------------------------------------------------------------------
function setRec3dSetDetailGUIValues(handles)

% Default initialization if Rec3dProject is empty
if isempty(handles.inputREC3DPROJECT) == 1
    set(handles.popupDetailMode, 'Value', 3);
    set(handles.popupDetailMode, 'Enable', 'off');
    set(handles.editOverviewSizeZ, 'String', num2str(0));
    set(handles.editOverviewSizeZ, 'Enable', 'off');
    set(handles.editOverviewPreBinning, 'String', num2str(0));
    set(handles.editOverviewPreBinning, 'Enable', 'off');
    set(handles.editOverviewPostBinning, 'String', num2str(0));
    set(handles.editOverviewPostBinning, 'Enable', 'off');
    set(handles.editNumberOfDetails, 'String', num2str(0));
    set(handles.editNumberOfDetails, 'Enable', 'off');
    set(handles.listboxDetailCoordinates, 'Enable', 'off');
    set(handles.pushbuttonGet, 'Enable', 'off');
    set(handles.pushbuttonClick, 'Enable', 'off');
    set(handles.pushbuttonClear, 'Enable', 'off');
    return;
end

% Deactivate Click if Rec3dReconstruction is empty
if isempty(handles.inputREC3DRECONSTRUCTION) == 1
    set(handles.pushbuttonClick, 'Enable', 'off');
end

% Simplify
% RECONSTRUCTION / Detail
DetailMode = handles.inputREC3DPROJECT.DETAIL.DetailMode;
OverviewSizeZ = handles.inputREC3DPROJECT.DETAIL.OverviewSizeZ;
OverviewPreBinning = handles.inputREC3DPROJECT.DETAIL.OverviewPreBinning;
OverviewPostBinning = handles.inputREC3DPROJECT.DETAIL.OverviewPostBinning;
NumberOfDetails = handles.inputREC3DPROJECT.DETAIL.NumberOfDetails;
DetailCoordinates = handles.inputREC3DPROJECT.DETAIL.DetailCoordinates;

% DetailMode
if strcmp(DetailMode, 'on') == 1
    set(handles.popupDetailMode, 'Value', 1);
elseif strcmp(DetailMode, 'av3') == 1
    set(handles.popupDetailMode, 'Value', 2);
elseif strcmp(DetailMode, 'off') == 1
    set(handles.popupDetailMode, 'Value', 3);
end

% OverviewSizeZ
set(handles.editOverviewSizeZ, 'String', num2str(OverviewSizeZ));

% OverviewPreBinning
set(handles.editOverviewPreBinning, 'String', num2str(OverviewPreBinning));

% OverviewPostBinning
set(handles.editOverviewPostBinning, 'String', num2str(OverviewPostBinning));

% NumberOfDetails
set(handles.editNumberOfDetails, 'String', num2str(NumberOfDetails));

% DetailCoordinates
set(handles.listboxDetailCoordinates, 'String', num2str(DetailCoordinates));
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% getRec3dSetDetailGUIValues
% -------------------------------------------------------------------------
function [handles] = getRec3dSetDetailGUIValues(handles)

% DetailMode
str = get(handles.popupDetailMode, 'String');
val = get(handles.popupDetailMode, 'Value');
DetailModeNew = str{val};

% OverviewSizeZ
OverviewSizeZNew = str2num(get(handles.editOverviewSizeZ, 'String'));

% OverviewPreBinning
OverviewPreBinningNew = str2num(get(handles.editOverviewPreBinning, 'String'));

% OverviewPostBinning
OverviewPostBinningNew = str2num(get(handles.editOverviewPostBinning, 'String'));

% NumberOfDetails
NumberOfDetailsNew = str2num(get(handles.editNumberOfDetails, 'String'));

% DetailCoordinates
DetailCoordinatesNew = str2num(get(handles.listboxDetailCoordinates, 'String'));

% Aktualize outputREC3DPROJECT
% RECONSTRUCTION / Detail
handles.outputREC3DPROJECT.DETAIL.DetailMode = DetailModeNew;
handles.outputREC3DPROJECT.DETAIL.OverviewSizeZ = OverviewSizeZNew;
handles.outputREC3DPROJECT.DETAIL.OverviewPreBinning = OverviewPreBinningNew;
handles.outputREC3DPROJECT.DETAIL.OverviewPostBinning = OverviewPostBinningNew;
handles.outputREC3DPROJECT.DETAIL.NumberOfDetails = NumberOfDetailsNew;
handles.outputREC3DPROJECT.DETAIL.DetailCoordinates = DetailCoordinatesNew;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% pushbuttonGet
% -------------------------------------------------------------------------
function pushbuttonGet_Callback(hObject, eventdata, handles)

% Simplify
% RECONSTRUCTION / Volume
SizeZ = handles.inputREC3DPROJECT.VOLUME.SizeZ;
PreBinning = handles.inputREC3DPROJECT.VOLUME.PreBinning;
PostBinning = handles.inputREC3DPROJECT.VOLUME.PostBinning;

% OverviewSizeZ
set(handles.editOverviewSizeZ, 'String', num2str(SizeZ));

% OverviewPreBinning
set(handles.editOverviewPreBinning, 'String', num2str(PreBinning));

% OverviewPostBinning
set(handles.editOverviewPostBinning, 'String', num2str(PostBinning));

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% pushbuttonClick
% -------------------------------------------------------------------------
function pushbuttonClick_Callback(hObject, eventdata, handles)

% Overview Reconstruction
OverviewReconstruction = handles.inputREC3DRECONSTRUCTION;

% Create av3_alignstruct
av3_alignstruct = tom_av3_createnewalignstruct;

% Click Details
av3_alignstructNew = tom_volxyz(OverviewReconstruction, av3_alignstruct);

% NumberOfDetailsNew
NumberOfDetailsNew = size(av3_alignstructNew, 2);

% DetailCoordinatesNew
if NumberOfDetailsNew == 0
    
    DetailCoordinatesNew = [];
    
else
    
    for k=1:NumberOfDetailsNew

         DetailCoordinatesNew(k,1) = av3_alignstructNew(1,k).Tomogram.Position.X;
         DetailCoordinatesNew(k,2) = av3_alignstructNew(1,k).Tomogram.Position.Y;
         DetailCoordinatesNew(k,3) = av3_alignstructNew(1,k).Tomogram.Position.Z;

    end

end

% Set NumberOfDetailsNew
set(handles.editNumberOfDetails, 'String', num2str(NumberOfDetailsNew));

% Set DetailCoordinatesNew
set(handles.listboxDetailCoordinates, 'String', num2str(DetailCoordinatesNew));

% Aktualize outputREC3DPROJECT
% RECONSTRUCTION / Detail
handles.outputREC3DPROJECT.DETAIL.av3_alignstruct = av3_alignstructNew;
                                                              
% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% pushbuttonClear
% -------------------------------------------------------------------------
function pushbuttonClear_Callback(hObject, eventdata, handles)

% DetailMode
set(handles.popupDetailMode, 'Value', 3);

% OverviewSizeZ
set(handles.editOverviewSizeZ, 'String', num2str(0));

% OverviewPreBinning
set(handles.editOverviewPreBinning, 'String', num2str(0));

% OverviewPostBinning
set(handles.editOverviewPostBinning, 'String', num2str(0));

% NumberOfDetails
set(handles.editNumberOfDetails, 'String', num2str(0));

% DetailCoordinates
set(handles.listboxDetailCoordinates, 'String', num2str([]));

% av3_alignstruct
handles.inputREC3DPROJECT.DETAIL.av3_alignstruct = [];
handles.outputREC3DPROJECT.DETAIL.av3_alignstruct = [];

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------










% -------------------------------------------------------------------------
% tom_av3_createnewalignstruct
% -------------------------------------------------------------------------
function Align = tom_av3_createnewalignstruct()

% Create new align structure
run = 1;
i = 1;
Align(run,i).Filename = '';
Align(run,i).Tomogram.Filename = '';
Align(run,i).Tomogram.Header = struct();
Align(run,i).Tomogram.Position.X = 0;
Align(run,i).Tomogram.Position.Y = 0;
Align(run,i).Tomogram.Position.Z = 0;
Align(run,i).Tomogram.Regfile = '';
Align(run,i).Tomogram.Offset = 0;
Align(run,i).Tomogram.Binning = 0;
Align(run,i).Tomogram.AngleMin = 0;
Align(run,i).Tomogram.AngleMax = 0;
Align(run,i).Shift.X = 0;
Align(run,i).Shift.Y = 0;
Align(run,i).Shift.Z = 0;
Align(run,i).Angle.Phi = 0;
Align(run,i).Angle.Psi = 0;
Align(run,i).Angle.Theta = 0;
Align(run,i).Angle.Rotmatrix = [];
Align(run,i).CCC = 0;
Align(run,i).Class = 0;
Align(run,i).ProjectionClass = 0;
Align(run,i).NormFlag = 0;
Align(run,i).Filter = [0 0];
% -------------------------------------------------------------------------

