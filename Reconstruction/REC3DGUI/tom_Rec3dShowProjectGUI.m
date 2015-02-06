function varargout = tom_Rec3dShowProjectGUI(varargin)
%TOM_REC3DSHOWPROJECTGUI is a module of TOM_REC3DGUI.
%
%   tom_Rec3dShowProjectGUI(Rec3dProject)
%
%   It displays project variables in the MATLAB Command Window.
%
%PARAMETERS
%
%  INPUT
%   Rec3dProject     main-structure of Rec3d
%  
%  OUTPUT
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
                   'gui_OpeningFcn', @tom_Rec3dShowProjectGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_Rec3dShowProjectGUI_OutputFcn, ...
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
function tom_Rec3dShowProjectGUI_OpeningFcn(hObject, eventdata, handles, varargin)

if (size(varargin, 2) < 1)
    
    handles.ProjectName = '';
    handles.ProjectStatus = '';
    handles.ProjectVersion = '';
    handles.PROJECT = [];
    handles.ALIGNMENT = [];
    handles.ALGRESIDUALS = [];
    handles.NAMEPATHEXT = [];
    handles.PARAMETER = [];
    handles.VOLUME = [];
    handles.METHOD = [];
    handles.FILTER = [];
    handles.DETAIL = [];
    handles.PARALLEL = [];
    handles.SIRTRESIDUALS = [];
    
else
    
    % Input REC3DPROJECT
    inputREC3DPROJECT = varargin{1};

    handles.ProjectName = inputREC3DPROJECT.ProjectName;
    handles.ProjectStatus = inputREC3DPROJECT.ProjectStatus;
    handles.ProjectVersion = inputREC3DPROJECT.ProjectVersion;
    handles.PROJECT = inputREC3DPROJECT.PROJECT;
    handles.ALIGNMENT = inputREC3DPROJECT.ALIGNMENT;
    handles.ALGRESIDUALS = inputREC3DPROJECT.ALGRESIDUALS;
    handles.NAMEPATHEXT = inputREC3DPROJECT.NAMEPATHEXT;
    handles.PARAMETER = inputREC3DPROJECT.PARAMETER;
    handles.VOLUME = inputREC3DPROJECT.VOLUME;
    handles.METHOD = inputREC3DPROJECT.METHOD;
    handles.FILTER = inputREC3DPROJECT.FILTER;
    handles.DETAIL = inputREC3DPROJECT.DETAIL;
    handles.PARALLEL = inputREC3DPROJECT.PARALLEL;
    handles.SIRTRESIDUALS = inputREC3DPROJECT.SIRTRESIDUALS;

end

% Update handles structure
guidata(hObject, handles);
uiwait;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% OutputFunction
% -------------------------------------------------------------------------
function varargout = tom_Rec3dShowProjectGUI_OutputFcn(hObject, eventdata, handles)

% Delete Rec3dShowProjectGUI
delete(handles.Rec3dShowProjectGUI);
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

% ShowProject
Rec3dShowProject(handles);

% Update handles structure
guidata(hObject, handles);
uiresume;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Rec3dShowProject
% -------------------------------------------------------------------------
function Rec3dShowProject(handles)

% Simplify
ProjectName = handles.ProjectName;
ProjectStatus = handles.ProjectStatus;
ProjectVersion = handles.ProjectVersion;
PROJECT = handles.PROJECT;
ALIGNMENT = handles.ALIGNMENT;
ALGRESIDUALS = handles.ALGRESIDUALS;
NAMEPATHEXT = handles.NAMEPATHEXT;
PARAMETER = handles.PARAMETER;
VOLUME = handles.VOLUME;
METHOD = handles.METHOD;
FILTER = handles.FILTER;
DETAIL = handles.DETAIL;
PARALLEL = handles.PARALLEL;
SIRTRESIDUALS = handles.SIRTRESIDUALS;

% Return if PROJECT is empty
if isempty(PROJECT) == 1
    return;
end

% Show PROJECT
valPROJECT = get(handles.checkboxPROJECT, 'Value');

if valPROJECT == 1
    
    % Show HEAD
    disp(' ');
    disp('ProjectName');
    disp(ProjectName);
    disp(' ');
    disp('ProjectStatus')
    disp(ProjectStatus);
    disp(' ');
    disp('ProjectVersion');
    disp(ProjectVersion);
    
    disp(' ');
    disp(' ');
    disp(' ');
    disp('PROJECT');
    
    % FieldNames
    FieldNames = fieldnames(PROJECT);
    % FieldValues
    for k = 1:size(FieldNames,1)
         FieldValues{k,1} = PROJECT.(FieldNames{k});
    end
    % Show FieldNames and FieldValues
    for k = 1:size(FieldNames,1)
         disp(' ');
         disp(FieldNames{k,1});
         disp(FieldValues{k,1});
    end

end
clear FieldNames;
clear FieldValues;

% Show ALIGNMENT
valALIGNMENT = get(handles.checkboxALIGNMENT, 'Value');

if valALIGNMENT == 1
    
    disp(' ');
    disp(' ');
    disp(' ');
    disp('ALIGNMENT');
    
    % FieldNames
    FieldNames = fieldnames(ALIGNMENT);
    % FieldValues
    for k = 1:size(FieldNames,1)
         FieldValues{k,1} = ALIGNMENT.(FieldNames{k});
    end
    % Show FieldNames and FieldValues
    for k = 1:size(FieldNames,1)
         disp(' ');
         disp(FieldNames{k,1});
         disp(FieldValues{k,1});
    end

end
clear FieldNames;
clear FieldValues;

% Show ALGRESIDUALS
valALGRESIDUALS = get(handles.checkboxALGRESIDUALS, 'Value');

if valALGRESIDUALS == 1
    
    disp(' ');
    disp(' ');
    disp(' ');
    disp('ALGRESIDUALS');
    
    % FieldNames
    FieldNames = fieldnames(ALGRESIDUALS);
    % FieldValues
    for k = 1:size(FieldNames,1)
         FieldValues{k,1} = ALGRESIDUALS.(FieldNames{k});
    end
    % Show FieldNames and FieldValues
    for k = 1:size(FieldNames,1)
         disp(' ');
         disp(FieldNames{k,1});
         disp(FieldValues{k,1});
    end

end
clear FieldNames;
clear FieldValues;

% Show NAMEPATHEXT
valNAMEPATHEXT = get(handles.checkboxNAMEPATHEXT, 'Value');

if valNAMEPATHEXT == 1
    
    disp(' ');
    disp(' ');
    disp(' ');
    disp('NAMEPATHEXT');
    
    % FieldNames
    FieldNames = fieldnames(NAMEPATHEXT);
    % FieldValues
    for k = 1:size(FieldNames,1)
         FieldValues{k,1} = NAMEPATHEXT.(FieldNames{k});
    end
    % Show FieldNames and FieldValues
    for k = 1:size(FieldNames,1)
         disp(' ');
         disp(FieldNames{k,1});
         disp(FieldValues{k,1});
    end

end
clear FieldNames;
clear FieldValues;

% Show PARAMETER
valPARAMETER = get(handles.checkboxPARAMETER, 'Value');

if valPARAMETER == 1
    
    disp(' ');
    disp(' ');
    disp(' ');
    disp('PARAMETER');
    
    % FieldNames
    FieldNames = fieldnames(PARAMETER);
    % FieldValues
    for k = 1:size(FieldNames,1)
         FieldValues{k,1} = PARAMETER.(FieldNames{k});
    end
    % Show FieldNames and FieldValues
    for k = 1:size(FieldNames,1)
         disp(' ');
         disp(FieldNames{k,1});
         disp(FieldValues{k,1});
    end

end
clear FieldNames;
clear FieldValues;

% Show VOLUME
valVOLUME = get(handles.checkboxVOLUME, 'Value');

if valVOLUME == 1
    
    disp(' ');
    disp(' ');
    disp(' ');
    disp('VOLUME');
    
    % FieldNames
    FieldNames = fieldnames(VOLUME);
    % FieldValues
    for k = 1:size(FieldNames,1)
         FieldValues{k,1} = VOLUME.(FieldNames{k});
    end
    % Show FieldNames and FieldValues
    for k = 1:size(FieldNames,1)
         disp(' ');
         disp(FieldNames{k,1});
         disp(FieldValues{k,1});
    end

end
clear FieldNames;
clear FieldValues;

% Show METHOD
valMETHOD = get(handles.checkboxMETHOD, 'Value');

if valMETHOD == 1
    
    disp(' ');
    disp(' ');
    disp(' ');
    disp('METHOD');
    
    % FieldNames
    FieldNames = fieldnames(METHOD);
    % FieldValues
    for k = 1:size(FieldNames,1)
         FieldValues{k,1} = METHOD.(FieldNames{k});
    end
    % Show FieldNames and FieldValues
    for k = 1:size(FieldNames,1)
         disp(' ');
         disp(FieldNames{k,1});
         disp(FieldValues{k,1});
    end

end
clear FieldNames;
clear FieldValues;

% Show FILTER
valFILTER = get(handles.checkboxFILTER, 'Value');

if valFILTER == 1
    
    disp(' ');
    disp(' ');
    disp(' ');
    disp('FILTER');
    
    % FieldNames
    FieldNames = fieldnames(FILTER);
    % FieldValues
    for k = 1:size(FieldNames,1)
         FieldValues{k,1} = FILTER.(FieldNames{k});
    end
    % Show FieldNames and FieldValues
    for k = 1:size(FieldNames,1)
         disp(' ');
         disp(FieldNames{k,1});
         disp(FieldValues{k,1});
    end

end
clear FieldNames;
clear FieldValues;

% Show DETAIL
valDETAIL = get(handles.checkboxDETAIL, 'Value');

if valDETAIL == 1
    
    disp(' ');
    disp(' ');
    disp(' ');
    disp('DETAIL');
    
    % FieldNames
    FieldNames = fieldnames(DETAIL);
    % FieldValues
    for k = 1:size(FieldNames,1)
         FieldValues{k,1} = DETAIL.(FieldNames{k});
    end
    % Show FieldNames and FieldValues
    for k = 1:size(FieldNames,1)
         disp(' ');
         disp(FieldNames{k,1});
         disp(FieldValues{k,1});
    end

end
clear FieldNames;
clear FieldValues;

% Show PARALLEL
valPARALLEL = get(handles.checkboxPARALLEL, 'Value');

if valPARALLEL == 1
    
    disp(' ');
    disp(' ');
    disp(' ');
    disp('PARALLEL');
    
    % FieldNames
    FieldNames = fieldnames(PARALLEL);
    % FieldValues
    for k = 1:size(FieldNames,1)
         FieldValues{k,1} = PARALLEL.(FieldNames{k});
    end
    % Show FieldNames and FieldValues
    for k = 1:size(FieldNames,1)
         disp(' ');
         disp(FieldNames{k,1});
         disp(FieldValues{k,1});
    end

end
clear FieldNames;
clear FieldValues;

% Show SIRTRESIDUALS
valSIRTRESIDUALS = get(handles.checkboxSIRTRESIDUALS, 'Value');

if valSIRTRESIDUALS == 1
    
    disp(' ');
    disp(' ');
    disp(' ');
    disp('SIRTRESIDUALS');
    
    % FieldNames
    FieldNames = fieldnames(SIRTRESIDUALS);
    % FieldValues
    for k = 1:size(FieldNames,1)
         FieldValues{k,1} = SIRTRESIDUALS.(FieldNames{k});
    end
    % Show FieldNames and FieldValues
    for k = 1:size(FieldNames,1)
         disp(' ');
         disp(FieldNames{k,1});
         disp(FieldValues{k,1});
    end

end
clear FieldNames;
clear FieldValues;
% -------------------------------------------------------------------------

