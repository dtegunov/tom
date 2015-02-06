function varargout = tom_Rec3dShowAlignmentGUI(varargin)
%TOM_REC3DSHOWALIGNMENTGUI is a module of TOM_REC3DGUI.
%
%   tom_Rec3dShowAlignmentGUI(Rec3dProject)
%
%   It displays the alignment.
%
%PARAMETERS
%
%  INPUT
%   Rec3dProject              main-structure of Rec3d
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


%--------------------------------------------------------------------------
% Begin initialization code - DO NOT EDIT
%--------------------------------------------------------------------------
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_Rec3dShowAlignmentGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_Rec3dShowAlignmentGUI_OutputFcn, ...
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
%--------------------------------------------------------------------------
% End initialization code - DO NOT EDIT
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% OpeningFunction
%--------------------------------------------------------------------------
function tom_Rec3dShowAlignmentGUI_OpeningFcn(hObject,eventdata,handles,varargin)

if (size(varargin, 2) < 1)
    
    % Input REC3DPROJECT
    handles.inputREC3DPROJECT = [];
    
else
    
    % Input REC3DPROJECT
    handles.inputREC3DPROJECT = varargin{1};
    
end

% Display alignment
setRec3dShowAlignmentGUIValues(handles);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% OutputFunction
%--------------------------------------------------------------------------
function varargout = tom_Rec3dShowAlignmentGUI_OutputFcn(hObject,eventdata,handles) 
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% setRec3dShowAlignmentGUIValues
% -------------------------------------------------------------------------
function setRec3dShowAlignmentGUIValues(handles)

% Return if inputREC3DPROJECT is empty
% Use GUI initialisation from guide
if isempty(handles.inputREC3DPROJECT) == 1
    return;
end

% Simplify
ProjectName = handles.inputREC3DPROJECT.ProjectName;
AlignmentMethod = handles.inputREC3DPROJECT.ALIGNMENT.AlignmentMethod;
ReferenceMarker = handles.inputREC3DPROJECT.ALIGNMENT.ReferenceMarker;
TiltingGeometry = handles.inputREC3DPROJECT.PROJECT.TiltingGeometry;
NumOfMarkers = handles.inputREC3DPROJECT.PROJECT.NumOfMarkers;
Imdim = handles.inputREC3DPROJECT.PROJECT.Imdim;

tx = handles.inputREC3DPROJECT.PARAMETER.tx;
ty = handles.inputREC3DPROJECT.PARAMETER.ty;
Tiltaxis = handles.inputREC3DPROJECT.PARAMETER.Tiltaxis;
isoscale = handles.inputREC3DPROJECT.PARAMETER.isoscale;

% Display ProjectName
set(handles.dispProjectName, 'String', ProjectName);
set(handles.dispAlignmentMethod, 'String', AlignmentMethod);
set(handles.dispReferenceMarker, 'String', num2str(ReferenceMarker));
set(handles.dispTiltingGeometry, 'String', TiltingGeometry);
set(handles.dispNumOfMarkers, 'String', NumOfMarkers);
set(handles.dispImdim, 'String', num2str(Imdim));

% Projection translation in x direction
axes(handles.plottx);plot(tx, 'g-');

% Projection translation in y direction
axes(handles.plotty);plot(ty, 'g-');

% Tiltaxis
axes(handles.plotTiltaxis);plot(Tiltaxis, 'r-');

% isoscale
axes(handles.plotisoscale);plot(isoscale, 'b-');
% -------------------------------------------------------------------------

