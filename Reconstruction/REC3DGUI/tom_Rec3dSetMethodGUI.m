function varargout = tom_Rec3dSetMethodGUI(varargin)
%TOM_REC3DSETMETHODGUI is a module of TOM_REC3DGUI.
%
%   METHOD = tom_Rec3dSetMethodGUI(METHOD)
%
%   It works on a structure which contains parameter for reconstruction.
%
%PARAMETERS
%
%  INPUT
%   METHOD              sub-structure of the main-structure of Rec3d
%                                METHOD = Rec3dProject.METHOD
%  
%  OUTPUT
%   METHOD              aktualized sub-structure
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
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
                   'gui_OpeningFcn', @tom_Rec3dSetMethodGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_Rec3dSetMethodGUI_OutputFcn, ...
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
function tom_Rec3dSetMethodGUI_OpeningFcn(hObject, eventdata, handles, varargin)

if (size(varargin, 2) < 1)
    
    % Input METHOD
    handles.inputMETHOD.ReconstructionMethod = 'WBP';
    handles.inputMETHOD.Normalization = 'phase';
    handles.inputMETHOD.SmoothBorders = 0;
    handles.inputMETHOD.AlignTiltaxis = 'ProjDir';
    handles.inputMETHOD.Handedness = 0;
    handles.inputMETHOD.ApplyWeighting = 'on';
    handles.inputMETHOD.WeightingMethod = 'exact';
    handles.inputMETHOD.ObjectThickness = 0;
    handles.inputMETHOD.Taper = 'off';
    handles.inputMETHOD.Iterations = 0;
    handles.inputMETHOD.Relaxation = 0;
    handles.inputMETHOD.Pathlength = 'ones';
    
else
    
    % Input METHOD
    handles.inputMETHOD = varargin{1};
    
end

% handles.inputMETHOD -> handles.outputMETHOD 
handles.outputMETHOD = handles.inputMETHOD;

% Display METHOD
setRec3dSetMethodGUIValues(handles);

% Update handles structure
guidata(hObject, handles);
uiwait;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% OutputFunction
% -------------------------------------------------------------------------
function varargout = tom_Rec3dSetMethodGUI_OutputFcn(hObject, eventdata, handles)

% Output METHOD
varargout{1} = handles.outputMETHOD;

% Delete Rec3dSetMethodGUI
delete(handles.Rec3dSetMethodGUI);
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

% Aktualize outputMETHOD
handles = getRec3dSetMethodGUIValues(handles);

% Update handles structure
guidata(hObject, handles);
uiresume;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% setRec3dSetMethodGUIValues
% -------------------------------------------------------------------------
function setRec3dSetMethodGUIValues(handles);

% Simplify
% RECONSTRUCTION / Method
ReconstructionMethod = handles.inputMETHOD.ReconstructionMethod;
Normalization = handles.inputMETHOD.Normalization;
SmoothBorders = handles.inputMETHOD.SmoothBorders;
AlignTiltaxis = handles.inputMETHOD.AlignTiltaxis;
Handedness = handles.inputMETHOD.Handedness;
ApplyWeighting = handles.inputMETHOD.ApplyWeighting;
WeightingMethod = handles.inputMETHOD.WeightingMethod;
ObjectThickness = handles.inputMETHOD.ObjectThickness;
Taper = handles.inputMETHOD.Taper;
Iterations = handles.inputMETHOD.Iterations;
Relaxation = handles.inputMETHOD.Relaxation;
Pathlength = handles.inputMETHOD.Pathlength;

% ReconstructionMethod
if strcmp(ReconstructionMethod, 'WBP') == 1
    set(handles.popupReconstructionMethod, 'Value', 1);
elseif strcmp(ReconstructionMethod, 'SIRT') == 1
    set(handles.popupReconstructionMethod, 'Value', 2);
end

% Normalization
if strcmp(Normalization, 'phase') == 1
    set(handles.popupNormalization, 'Value', 1);
elseif strcmp(Normalization, 'off') == 1
    set(handles.popupNormalization, 'Value', 2);
end

% SmoothBorders
set(handles.editSmoothBorders, 'String', num2str(SmoothBorders));

% AlignTiltaxis
if strcmp(AlignTiltaxis, 'ProjDir') == 1
    set(handles.popupAlignTiltaxis, 'Value', 1);
elseif strcmp(AlignTiltaxis, 'Y-axis') == 1
    set(handles.popupAlignTiltaxis, 'Value', 2);
end

% Handedness
if isequal(Handedness, 0) == 1
    set(handles.popupHandedness, 'Value', 1);
elseif isequal(Handedness, 180) == 1
    set(handles.popupHandedness, 'Value', 2);
end

% ApplyWeighting
if strcmp(ApplyWeighting, 'on') == 1
    set(handles.popupApplyWeighting, 'Value', 1);
elseif strcmp(ApplyWeighting, 'off') == 1
    set(handles.popupApplyWeighting, 'Value', 2);
end

% WeightingMethod
if strcmp(WeightingMethod, 'exact') == 1
    set(handles.popupWeightingMethod, 'Value', 1);
elseif strcmp(WeightingMethod, 'analytical') == 1
    set(handles.popupWeightingMethod, 'Value', 2);
elseif strcmp(WeightingMethod, '3d') == 1
    set(handles.popupWeightingMethod, 'Value', 3);
end

% ObjectThickness
set(handles.editObjectThickness, 'String', num2str(ObjectThickness));

% Taper
if strcmp(Taper, 'off') == 1
    set(handles.checkboxTaper, 'Value', 0);
elseif strcmp(Taper, 'on') == 1
    set(handles.checkboxTaper, 'Value', 1);
end

% Iterations
set(handles.editIterations, 'String', num2str(Iterations));

% Relaxation
set(handles.editRelaxation, 'String', num2str(Relaxation));

% Pathlength
if strcmp(Pathlength, 'ones') == 1
    set(handles.popupPathlength, 'Value', 1);
elseif strcmp(Pathlength, 'slab') == 1
    set(handles.popupPathlength, 'Value', 2);
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% getRec3dSetMethodGUIValues
% -------------------------------------------------------------------------
function [handles] = getRec3dSetMethodGUIValues(handles)

% ReconstructionMethod
str = get(handles.popupReconstructionMethod, 'String');
val = get(handles.popupReconstructionMethod, 'Value');
ReconstructionMethodNew = str{val};

% Normalization
str = get(handles.popupNormalization, 'String');
val = get(handles.popupNormalization, 'Value');
NormalizationNew = str{val};

% SmoothBorders
SmoothBordersNew = str2num(get(handles.editSmoothBorders, 'String'));

% AlignTiltaxis
str = get(handles.popupAlignTiltaxis, 'String');
val = get(handles.popupAlignTiltaxis, 'Value');
AlignTiltaxisNew = str{val};

% Handedness
str = get(handles.popupHandedness, 'String');
val = get(handles.popupHandedness, 'Value');
HandednessNew = str2num(str{val});

% ApplyWeighting
str = get(handles.popupApplyWeighting, 'String');
val = get(handles.popupApplyWeighting, 'Value');
ApplyWeightingNew = str{val};

% WeightingMethod
str = get(handles.popupWeightingMethod, 'String');
val = get(handles.popupWeightingMethod, 'Value');
WeightingMethodNew = str{val};

% ObjectThickness
ObjectThicknessNew = str2num(get(handles.editObjectThickness, 'String'));

% Taper
TaperNew = get(handles.checkboxTaper, 'Value');
if isequal(TaperNew, 0) == 1
    TaperNew = 'off';
elseif  isequal(TaperNew, 1) == 1
    TaperNew = 'on';
end

% Iterations
IterationsNew = str2num(get(handles.editIterations, 'String'));

% Relaxation
RelaxationNew = str2num(get(handles.editRelaxation, 'String'));

% Pathlength
str = get(handles.popupPathlength, 'String');
val = get(handles.popupPathlength, 'Value');
PathlengthNew = str{val};

% Aktualize outputMETHOD
handles.outputMETHOD.ReconstructionMethod = ReconstructionMethodNew;
handles.outputMETHOD.Normalization = NormalizationNew;
handles.outputMETHOD.SmoothBorders = SmoothBordersNew;
handles.outputMETHOD.AlignTiltaxis = AlignTiltaxisNew;
handles.outputMETHOD.Handedness = HandednessNew;
handles.outputMETHOD.ApplyWeighting = ApplyWeightingNew;
handles.outputMETHOD.WeightingMethod = WeightingMethodNew;
handles.outputMETHOD.ObjectThickness = ObjectThicknessNew;
handles.outputMETHOD.Taper = TaperNew;
handles.outputMETHOD.Iterations = IterationsNew;
handles.outputMETHOD.Relaxation = RelaxationNew;
handles.outputMETHOD.Pathlength = PathlengthNew;
% -------------------------------------------------------------------------

