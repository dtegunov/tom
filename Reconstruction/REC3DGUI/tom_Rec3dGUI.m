function varargout = tom_Rec3dGUI(varargin)
%TOM_REC3DGUI is a tool for tomographic reconstruction.
%
%   tom_Rec3dGUI
%
%   It is naturally designed for both >>singleaxis<< and 
%   >>dualaxis<< Tilting Geometry and for both WBP and 
%   SIRT reconstruction algorithms.
%
%PARAMETERS
%
%  INPUT
%   TILTSERIES                 Select a singleaxis or a dualaxis tiltseries
%                                      for reconstruction with the Generate Project dialog.
%   MARKERFILE              Select a singleaxis or a dualaxis markerfile
%                                      for alignment with the Generate Project dialog.
%
%  OUTPUT
%   RECONSTRUCTION   3-d reconstruction
%   REC3DPROJECT         structure with all reconstruction settings
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%
%   created by ME 01/01/07
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
                   'gui_OpeningFcn', @tom_Rec3dGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_Rec3dGUI_OutputFcn, ...
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
function tom_Rec3dGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Specify application-defined data
setappdata(0,'UseNativeSystemDialogs',0);
%rmappdata(0,'UseNativeSystemDialogs');

% Initialize Rec3dControl
    %%%%%% Aktualize Rec3dControl here
    Rec3dControl.WaitbarMode = 'on';
    Rec3dControl.ProtocolMode = 'off';
    Rec3dControl.DemoMode = 'off';
    Rec3dControl.DelMode = 'on';
    Rec3dControl.ParallelMode = 'off';
    Rec3dControl.DistributedComputingToolbox = ...
                 CheckDistributedComputingToolbox();
    Rec3dControl.RecInMemory = 'no';
    Rec3dControl.Rec3dVersion = '1.0';

% Initialize Rec3dProject
Rec3dProject = [];

% Initialize Rec3dReconstruction
Rec3dReconstruction = [];

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles...
         (handles, Rec3dControl, Rec3dProject, Rec3dReconstruction);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% OutputFunction
% -------------------------------------------------------------------------
function varargout = tom_Rec3dGUI_OutputFcn(hObject, eventdata, handles) 
% -------------------------------------------------------------------------










% -------------------------------------------------------------------------
% PROJECT
% -------------------------------------------------------------------------
function Project_Callback(hObject, eventdata, handles)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% LoadProject
% -------------------------------------------------------------------------
function LoadProject_Callback(hObject, eventdata, handles)

% handles -> Rec3d
Rec3dControl = handles2Rec3d(handles);

% Get ProjectFileName and ProjectFilePath              
[ProjectFileName, ProjectFilePath] = uigetfile({'*.mat', 'MAT-files (*.mat)';}, 'Load Project'); 

if isequal(ProjectFileName, 0) || isequal(ProjectFilePath, 0)
    return;
end

% Load Project
workingdir = pwd;
cd (ProjectFilePath);
Rec3dProject = load(ProjectFileName);
cd (workingdir);

% Rec3dProject.Rec3dProject -> Rec3dProject
Rec3dProject = Rec3dProject.Rec3dProject;

% Check ProjectVersion
% Only needed for MPI internal Rec3dGUI version
if strcmp(Rec3dProject.ProjectVersion, '1.0') ~= 1
    % ERROR
    msgbox('This Rec3dProject can not be loaded in version 1.0', ...
                   'Load Project', 'error');
    return;
end

% Clear Rec3dReconstruction
Rec3dReconstruction = [];

% Aktualize 'RecInMemory'
Rec3dControl.RecInMemory = 'no';

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles...
          (handles, Rec3dControl, Rec3dProject, Rec3dReconstruction);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% GenerateProject
% -------------------------------------------------------------------------
function GenerateProject_Callback(hObject, eventdata, handles)

% handles -> Rec3d
Rec3dControl = handles2Rec3d(handles);

% Rec3dGenerateProject
Rec3dProject = Rec3dGenerateProject();

% Error return
if isempty(Rec3dProject) == 1
    return;
end

% Clear Rec3dReconstruction
Rec3dReconstruction = [];

% Aktualize 'RecInMemory'
Rec3dControl.RecInMemory = 'no';

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles...
          (handles, Rec3dControl, Rec3dProject, Rec3dReconstruction);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ModifyProject
% -------------------------------------------------------------------------
function ModifyProject_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% ModifyProject
Rec3dProject = tom_Rec3dModifyProjectGUI(Rec3dProject);

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SaveProject
% -------------------------------------------------------------------------
function SaveProject_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Put ProjectFileName and ProjectFilePath              
[ProjectFileName, ProjectFilePath] = uiputfile({'*.mat', 'MAT-files (*.mat)';}, ...
                            'Save Project', [Rec3dProject.ProjectName '.mat']); 

if isequal(ProjectFileName, 0) || isequal(ProjectFilePath, 0)
    return;
end

% Save Project
workingdir = pwd;
cd (ProjectFilePath);
save(ProjectFileName, 'Rec3dProject');
cd (workingdir);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% WaitbarMode
% -------------------------------------------------------------------------
function WaitbarMode_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Change checkbox
if strcmp(Rec3dControl.WaitbarMode, 'off')
    Rec3dControl.WaitbarMode = 'on';
else
    Rec3dControl.WaitbarMode = 'off';
end

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ProtocolMode
% -------------------------------------------------------------------------
function ProtocolMode_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Change checkbox
if strcmp(Rec3dControl.ProtocolMode, 'off')
    Rec3dControl.ProtocolMode = 'on';
else
    Rec3dControl.ProtocolMode = 'off';
end

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% DemoMode
% -------------------------------------------------------------------------
function DemoMode_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Change checkbox
if strcmp(Rec3dControl.DemoMode, 'off')
    Rec3dControl.DemoMode = 'on';
else
    Rec3dControl.DemoMode = 'off';
end

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% DelMode
% -------------------------------------------------------------------------
function DelMode_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Change checkbox
if strcmp(Rec3dControl.DelMode, 'off')
    Rec3dControl.DelMode = 'on';
else
    Rec3dControl.DelMode = 'off';
end

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ParallelMode
% -------------------------------------------------------------------------
function ParallelMode_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Change checkbox
if strcmp(Rec3dControl.ParallelMode, 'off')
    Rec3dControl.ParallelMode = 'on';
else
    Rec3dControl.ParallelMode = 'off';
end

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ShowProjections
% -------------------------------------------------------------------------
function ShowProjections_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% No Rec3dProject loaded
if isempty(Rec3dProject) == 1
    tom_embrowse;
    return;
end

% Simplify PROJECT
TiltingGeometry  = Rec3dProject.PROJECT.TiltingGeometry;
ProjectionsPath1 = Rec3dProject.PROJECT.ProjectionsPath1;
ProjectionsPath2 = Rec3dProject.PROJECT.ProjectionsPath2;

% Show projections of your 'singleaxis' project
if strcmp(TiltingGeometry, 'singleaxis') == 1 
    workingdir = pwd;
    cd (ProjectionsPath1);
    tom_embrowse;
    cd (workingdir);
    return;
end

% Show projections of your 'dualaxis' project
if strcmp(TiltingGeometry, 'dualaxis') == 1 
    
    % Which Tiltseries
    Tiltseries = questdlg('Please choose tiltseries:', ...
    'Show Projections', 'Tiltseries1', 'Tiltseries2', 'Cancel', 'Tiltseries1');
                 
    if strcmp(Tiltseries, 'Cancel') == 1 || strcmp(Tiltseries, '') == 1
        return;
    end
    
    if strcmp(Tiltseries, 'Tiltseries1') == 1
         workingdir = pwd;
         cd (ProjectionsPath1);
         tom_embrowse;
         cd (workingdir);
         return;
    else
         workingdir = pwd;
         cd (ProjectionsPath2);
         tom_embrowse;
         cd (workingdir);
         return;
    end

end

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ShowProject
% -------------------------------------------------------------------------
function ShowProject_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% ShowProject
tom_Rec3dShowProjectGUI(Rec3dProject);

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% CopyProject
% -------------------------------------------------------------------------
function CopyProject_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Copy Control
assignin('base', 'Rec3dControl', Rec3dControl);

% Copy Project
assignin('base', 'Rec3dProject', Rec3dProject);

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Quit
% -------------------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)

QuitValue = questdlg('Do you really want to quit?', ...
                     'Quit','Yes','No','Yes');
                 
if strcmp(QuitValue, 'No') == 1 || strcmp(QuitValue, '') == 1
    return;
else
    close(handles.Rec3dGUI);
end
% -------------------------------------------------------------------------










% -------------------------------------------------------------------------
% ALIGNMENT
% -------------------------------------------------------------------------
function Alignment_Callback(hObject, eventdata, handles)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% DoAlignment
% -------------------------------------------------------------------------
function DoAlignment_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% AlignmentEngine
Rec3dProject = tom_Rec3dAlignmentEngine(Rec3dControl, Rec3dProject);

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ClearAlignment
% -------------------------------------------------------------------------
function ClearAlignment_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Simplify ALIGNMENT
AlignmentMethod = Rec3dProject.ALIGNMENT.AlignmentMethod;
ReferenceMarker = Rec3dProject.ALIGNMENT.ReferenceMarker;

% Clear
% Clear HEAD
Rec3dProject.ProjectStatus = 'loaded';
% Clear ALIGNMENT
Rec3dProject.ALIGNMENT.AlignmentMethod = AlignmentMethod;
Rec3dProject.ALIGNMENT.ReferenceMarker = ReferenceMarker;
Rec3dProject.ALIGNMENT.Origin1 = [0 0 0];
Rec3dProject.ALIGNMENT.m3d1 = 0;
Rec3dProject.ALIGNMENT.Tiltaxis1 = 0;
Rec3dProject.ALIGNMENT.tx1 = 0;
Rec3dProject.ALIGNMENT.ty1 = 0;
Rec3dProject.ALIGNMENT.isoscale1 = 1;
Rec3dProject.ALIGNMENT.WarpDone1 = 'no';
Rec3dProject.ALIGNMENT.WarpAlignment1 = 0;
Rec3dProject.ALIGNMENT.Origin2 = [0 0 0];
Rec3dProject.ALIGNMENT.m3d2 = 0;
Rec3dProject.ALIGNMENT.Tiltaxis2 = 0;
Rec3dProject.ALIGNMENT.tx2 = 0;
Rec3dProject.ALIGNMENT.ty2 = 0;
Rec3dProject.ALIGNMENT.isoscale2 = 1;
Rec3dProject.ALIGNMENT.WarpDone2 = 'no';
Rec3dProject.ALIGNMENT.WarpAlignment2 = 0;
Rec3dProject.ALIGNMENT.RotMatrix = 0;
Rec3dProject.ALIGNMENT.Psi = 0;
Rec3dProject.ALIGNMENT.Theta = 0;
Rec3dProject.ALIGNMENT.Phi = 0;
% Clear ALGRESIDUALS
Rec3dProject.ALGRESIDUALS.ResidualMatrix1 = 0;
Rec3dProject.ALGRESIDUALS.Sigma1 = 0;
Rec3dProject.ALGRESIDUALS.AveragePerProjection1 = 0;
Rec3dProject.ALGRESIDUALS.AveragePerMarker1 = 0;
Rec3dProject.ALGRESIDUALS.ResidualMatrix2 = 0;
Rec3dProject.ALGRESIDUALS.Sigma2 = 0;
Rec3dProject.ALGRESIDUALS.AveragePerProjection2 = 0;
Rec3dProject.ALGRESIDUALS.AveragePerMarker2 = 0;
Rec3dProject.ALGRESIDUALS.EulerAnglesResidual = 0;
Rec3dProject.ALGRESIDUALS.MaximumResidual = 0;
Rec3dProject.ALGRESIDUALS.ResidualSpheres = 0;
Rec3dProject.ALGRESIDUALS.AverageResidualSphere = 0;
% Clear PARAMETER
Rec3dProject.PARAMETER.Tiltangles = 0;
Rec3dProject.PARAMETER.Tiltaxis = 0;
Rec3dProject.PARAMETER.ProjDir = 0;
Rec3dProject.PARAMETER.tx = 0;
Rec3dProject.PARAMETER.ty = 0;
Rec3dProject.PARAMETER.isoscale = 1;

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ShowAlignment
% -------------------------------------------------------------------------
function ShowAlignment_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% ShowAlignmentGUI
tom_Rec3dShowAlignmentGUI(Rec3dProject);

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% DoWarpAlignment
% -------------------------------------------------------------------------
function DoWarpAlignment_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% ERROR
msgbox('Warp Alignment not implemented!', ...
               'Do Warp Alignment', 'error');

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------










% -------------------------------------------------------------------------
% RECONSTRUCTION
% -------------------------------------------------------------------------
function Reconstruction_Callback(hObject, eventdata, handles)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SetMethod
% -------------------------------------------------------------------------
function SetMethod_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Rec3dProject.METHOD -> METHOD
METHOD = Rec3dProject.METHOD;

% SetMethod
METHOD = tom_Rec3dSetMethodGUI(METHOD);

% METHOD -> Rec3dProject.METHOD 
Rec3dProject.METHOD = METHOD;

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SetFilter
% -------------------------------------------------------------------------
function SetFilter_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Rec3dProject.FILTER -> FILTER
FILTER = Rec3dProject.FILTER;

% SetFilter
FILTER = tom_Rec3dSetFilterGUI(FILTER);

% FILTER -> Rec3dProject.FILTER 
Rec3dProject.FILTER = FILTER;

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SetDetail
% -------------------------------------------------------------------------
function SetDetail_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject, Rec3dReconstruction] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% SetDetail
Rec3dProject = tom_Rec3dSetDetailGUI(Rec3dProject, Rec3dReconstruction);

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SetParallel
% -------------------------------------------------------------------------
function SetParallel_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% ERROR
msgbox('Set Parallel not implemented!', ...
               'Set Parallel', 'error');

% % Rec3dProject.PARALLEL -> PARALLEL
% PARALLEL = Rec3dProject.PARALLEL;
% 
% % SetParallel
% PARALLELNEW = tom_parallelsettings(PARALLEL);
% 
% if strcmp(PARALLELNEW, '') == 1 % Cancel clicked in tom_parallelsettings
%     PARALLELNEW = PARALLEL;
% end
% 
% % PARALLELNEW -> Rec3dProject.PARALLEL 
% Rec3dProject.PARALLEL = PARALLELNEW;

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% DoReconstruction
% -------------------------------------------------------------------------
function DoReconstruction_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Enter ReconstructionName and locate ReconstructionPath
[ReconstructionNamelong, ReconstructionPath] = ...
    uiputfile({'*.vol','VOL-files (*.vol)';}, ...
'Enter Reconstruction Name and locate Reconstruction Path', ...
    [Rec3dProject.ProjectName '.vol']);

if isequal(ReconstructionNamelong, 0) || isequal(ReconstructionPath, 0)
    return;
end

% Get ReconstructionExt
[pathstr, ReconstructionName, ReconstructionExt] = ...
                                       fileparts(ReconstructionNamelong);

% Save ReconstructionName, ReconstructionPath, ReconstructionExt
Rec3dProject.NAMEPATHEXT.ReconstructionName = ReconstructionName;
Rec3dProject.NAMEPATHEXT.ReconstructionPath = ReconstructionPath;
Rec3dProject.NAMEPATHEXT.ReconstructionExt = ReconstructionExt;

% Rec3dProject.METHOD.ReconstructionMethod -> ReconstructionMethod
ReconstructionMethod = Rec3dProject.METHOD.ReconstructionMethod;

% Rec3dControl.ParallelMode -> ParallelMode
ParallelMode = Rec3dControl.ParallelMode;

% WBP
if strcmp(ReconstructionMethod, 'WBP') == 1
    
    % WBP / stand-alone computer
    if strcmp(ParallelMode, 'off') == 1
         [Rec3dReconstruction, Rec3dProject] = ...
         tom_Rec3dWBPReconstructionEngine(Rec3dControl, Rec3dProject, 'rec', '001', 'guistyle');
    end
    
    % WBP / parallel computer
    if strcmp(ParallelMode, 'on') == 1
         % ERROR
         msgbox('PARALLEL-WBP not implemented!', ...
              'Do Reconstruction', 'error');
         Rec3dReconstruction = [];   
    end
    
end

% SIRT
if strcmp(ReconstructionMethod, 'SIRT') == 1
    
    % SIRT / stand-alone computer
    if strcmp(ParallelMode, 'off') == 1
         [Rec3dReconstruction, Rec3dProject] = ...
         tom_Rec3dSIRTReconstructionEngine(Rec3dControl, Rec3dProject, '001', 'guistyle');
    end
    
    % SIRT / parallel computer largearray
    if strcmp(ParallelMode, 'on') == 1
         % ERROR
         msgbox('PARALLEL-SIRT not implemented!', ...
              'Do Reconstruction', 'error');
         Rec3dReconstruction = [];   
    end
    
end

% Error return
if isempty(Rec3dReconstruction) == 1
    return;
end

% Rec3dProject.DETAIL.DetailMode -> DetailMode
DetailMode = Rec3dProject.DETAIL.DetailMode;

% Rec3dProject.DETAIL.NumberOfDetails -> NumberOfDetails
NumberOfDetails = Rec3dProject.DETAIL.NumberOfDetails;
    
% Aktualize ReconstructionName for DetailMode
if (strcmp(DetailMode, 'on') == 1 || strcmp(DetailMode, 'av3') == 1) && NumberOfDetails >= 1
    Rec3dProject.NAMEPATHEXT.ReconstructionName = [ReconstructionName '_' num2str(NumberOfDetails)];
end

% Aktualize 'RecInMemory'
Rec3dControl.RecInMemory = 'yes';

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles...
          (handles, Rec3dControl, Rec3dProject, Rec3dReconstruction);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ClearReconstruction
% -------------------------------------------------------------------------
function ClearReconstruction_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Clear
% Clear HEAD
Rec3dProject.ProjectStatus = 'aligned';
% Clear RECONSTRUCTION
Rec3dProject.NAMEPATHEXT.ReconstructionName = '';
Rec3dProject.NAMEPATHEXT.ReconstructionPath = '';
Rec3dProject.NAMEPATHEXT.ReconstructionExt = '';
Rec3dProject.NAMEPATHEXT.TempFilesName = '';
Rec3dProject.NAMEPATHEXT.TempFilesPath = '';
Rec3dProject.NAMEPATHEXT.TempFilesExt = '';
% Clear SIRTRESIDUALS
Rec3dProject.SIRTRESIDUALS.EuclidianDistancePerProjection = 0;
Rec3dProject.SIRTRESIDUALS.EuclidianDistancePerIteration = 0;
Rec3dProject.SIRTRESIDUALS.dProjVal = 0;
Rec3dProject.SIRTRESIDUALS.DifVolVal = 0;
Rec3dProject.SIRTRESIDUALS.RecVolVal = 0;

% Clear Rec3dReconstruction
Rec3dReconstruction = [];

% Aktualize 'RecInMemory'
Rec3dControl.RecInMemory = 'no';

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles...
          (handles, Rec3dControl, Rec3dProject, Rec3dReconstruction);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% LoadReconstruction
% -------------------------------------------------------------------------
function LoadReconstruction_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% Simplify HEAD
ReconstructionName = Rec3dProject.NAMEPATHEXT.ReconstructionName;
ReconstructionPath = Rec3dProject.NAMEPATHEXT.ReconstructionPath;
ReconstructionExt = Rec3dProject.NAMEPATHEXT.ReconstructionExt;

% Reconstruction not saved
if strcmp(ReconstructionPath, '') == 1
    % ERROR
    msgbox('Reconstruction not saved!', 'Load Reconstruction', 'error');
    return;
end

% Check directory already exists
CheckDirectory = isdir(ReconstructionPath);

if CheckDirectory == 0
    % ERROR
    msgbox('Reconstruction directory lost!', 'Load Reconstruction', 'error');
    return;
end

% Check reconstruction already exists
workingdir = pwd;
cd (ReconstructionPath);
CheckReconstruction = tom_isemfile([ReconstructionName ReconstructionExt]);
cd (workingdir);

if CheckReconstruction == 0 || CheckReconstruction == -1
    % ERROR
    msgbox('Reconstruction lost!', 'Load Reconstruction', 'error');
    return;
end

% Load reconstruction
workingdir = pwd;
cd (ReconstructionPath);
Rec3dReconstruction = tom_emread([ReconstructionName ReconstructionExt]);
cd (workingdir);

% Clear reconstruction header
Rec3dReconstruction = Rec3dReconstruction.Value;

% Aktualize 'RecInMemory'
Rec3dControl.RecInMemory = 'yes';

% Reconstruction loaded
% MESSAGE
uiwait(msgbox('Reconstruction loaded!', 'Load Reconstruction', 'warn'));

% setRec3dGUI
setRec3dGUI(handles, Rec3dControl, Rec3dProject);

% Rec3d -> handles
handles = Rec3d2handles...
          (handles, Rec3dControl, Rec3dProject, Rec3dReconstruction);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ShowReconstruction
% -------------------------------------------------------------------------
function ShowReconstruction_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject, Rec3dReconstruction] = handles2Rec3d(handles);

tom_volxyz(Rec3dReconstruction);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SaveReconstruction
% -------------------------------------------------------------------------
function SaveReconstruction_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject, Rec3dReconstruction] = handles2Rec3d(handles);

% Enter ReconstructionName and locate ReconstructionPath
[ReconstructionNamelong, ReconstructionPath] = ...
    uiputfile({'*.vol','VOL-files (*.vol)';}, ...
'Enter Reconstruction Name and locate Reconstruction Path', ...
    [Rec3dProject.ProjectName '.vol']);

if isequal(ReconstructionNamelong, 0) || isequal(ReconstructionPath, 0)
    return;
end

% Get ReconstructionExt
[pathstr, ReconstructionName, ReconstructionExt] = ...
                                       fileparts(ReconstructionNamelong);

% Save ReconstructionName, ReconstructionPath, ReconstructionExt
Rec3dProject.NAMEPATHEXT.ReconstructionName = ReconstructionName;
Rec3dProject.NAMEPATHEXT.ReconstructionPath = ReconstructionPath;
Rec3dProject.NAMEPATHEXT.ReconstructionExt = ReconstructionExt;

% Save reconstruction
workingdir = pwd;
cd (ReconstructionPath);
tom_emwrite([ReconstructionName ReconstructionExt], Rec3dReconstruction);
cd (workingdir);

% Reconstruction saved
% MESSAGE
uiwait(msgbox('Reconstruction saved!', 'Save Reconstruction', 'warn'));

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% CreateTempFiles
% -------------------------------------------------------------------------
function CreateTempFiles_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% getRec3dGUI
[Rec3dControl, Rec3dProject] = getRec3dGUI...
                              (handles, Rec3dControl, Rec3dProject);

% Enter TempFilesName and locate TempFilesPath
[TempFilesNamelong, TempFilesPath] = ...
    uiputfile({'*.em','EM-files (*.em)';}, ...
'Enter TempFiles Name and locate TempFiles Path', ...
    [Rec3dProject.ProjectName '_']);

if isequal(TempFilesNamelong, 0) || isequal(TempFilesPath, 0)
    return;
end

% Get TempFilesName and TempFilesExt
[pathstr, TempFilesName, TempFilesExt] = ...
                                       fileparts(TempFilesNamelong);

% Save TempFilesName, TempFilesPath, TempFilesExt
Rec3dProject.NAMEPATHEXT.TempFilesName = TempFilesName;
Rec3dProject.NAMEPATHEXT.TempFilesPath = TempFilesPath;
Rec3dProject.NAMEPATHEXT.TempFilesExt = TempFilesExt;

% Rec3dControl.ParallelMode -> ParallelMode
ParallelMode = Rec3dControl.ParallelMode;

% Create TempFiles
if strcmp(ParallelMode, 'off') == 1
    tom_Rec3dWBPReconstructionEngine(Rec3dControl, Rec3dProject, 'temp', '001', 'guistyle');
end

% PARALLEL TempFile Engine
if strcmp(ParallelMode, 'on') == 1
    % ERROR
    msgbox('PARALLEL-TempFile Engine not implemented!', ...
         'Create TempFiles', 'error');
end

% Rec3d -> handles
handles = Rec3d2handles(handles, Rec3dControl, Rec3dProject);

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% PreBinning
% -------------------------------------------------------------------------
function popupPreBinning_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% Rec3dProject loaded
if isempty(Rec3dProject) == 0   
    
    % Simplify Rec3dProject.ALPHA.beta -> beta
    % PROJECT
    Imdim = Rec3dProject.PROJECT.Imdim;
    
    % Get current PostBinning
    PostBinning = get(handles.popupPostBinning, 'Value')-1;
    
    val = get(hObject,'Value');
    str = get(hObject, 'String');
    switch str{val};
    case '0'
         PreBinnedImdim = Imdim/(2^(0+PostBinning));
    case '1'
         PreBinnedImdim = Imdim/(2^(1+PostBinning));
    case '2'
         PreBinnedImdim = Imdim/(2^(2+PostBinning));
    case '3'
         PreBinnedImdim = Imdim/(2^(3+PostBinning));
    case '4'
         PreBinnedImdim = Imdim/(2^(4+PostBinning));
    case '5'
         PreBinnedImdim = Imdim/(2^(5+PostBinning));
    case '6'
         PreBinnedImdim = Imdim/(2^(6+PostBinning));
    case '7'
         PreBinnedImdim = Imdim/(2^(7+PostBinning));
    end
    
    set(handles.editSizeX, 'String', num2str(PreBinnedImdim));
    set(handles.editSizeY, 'String', num2str(PreBinnedImdim));
    set(handles.editSizeZ, 'String', num2str(PreBinnedImdim));
    
end

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% PostBinning
% -------------------------------------------------------------------------
function popupPostBinning_Callback(hObject, eventdata, handles)

% handles -> Rec3d
[Rec3dControl, Rec3dProject] = handles2Rec3d(handles);

% Rec3dProject loaded
if isempty(Rec3dProject) == 0
    
    % Simplify Rec3dProject.ALPHA.beta -> beta
    % PROJECT
    Imdim = Rec3dProject.PROJECT.Imdim;
    
    % Get current PreBinning
    PreBinning = get(handles.popupPreBinning, 'Value')-1;
    
    val = get(hObject,'Value');
    str = get(hObject, 'String');
    switch str{val};
    case '0'
         PostBinnedImdim = Imdim/(2^(0+PreBinning));
    case '1'
         PostBinnedImdim = Imdim/(2^(1+PreBinning));
    case '2'
         PostBinnedImdim = Imdim/(2^(2+PreBinning));
    case '3'
         PostBinnedImdim = Imdim/(2^(3+PreBinning));
    case '4'
         PostBinnedImdim = Imdim/(2^(4+PreBinning));
    case '5'
         PostBinnedImdim = Imdim/(2^(5+PreBinning));
    case '6'
         PostBinnedImdim = Imdim/(2^(6+PreBinning));
    case '7'
         PostBinnedImdim = Imdim/(2^(7+PreBinning));
    end
    
    set(handles.editSizeX, 'String', num2str(PostBinnedImdim));
    set(handles.editSizeY, 'String', num2str(PostBinnedImdim));
    set(handles.editSizeZ, 'String', num2str(PostBinnedImdim));

end

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------










% -------------------------------------------------------------------------
% Check Distributed Computing Toolbox
% -------------------------------------------------------------------------
function [YesOrNo] = CheckDistributedComputingToolbox()

% CheckDCT = ver;
% YesOrNo = 'no';
% for k=1:size(CheckDCT,2)
%     if strcmp(CheckDCT(1,k).Name, 'Distributed Computing Toolbox') == 1
%         YesOrNo = 'yes';
%     end
% end
YesOrNo = 'yes';
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% handles2Rec3d
% -------------------------------------------------------------------------
function [Rec3dControl, Rec3dProject, Rec3dReconstruction] = ...
handles2Rec3d(handles)

% handles -> Rec3d
if nargout == 1
    Rec3dControl = handles.Rec3dControl;
elseif nargout == 2
    Rec3dControl = handles.Rec3dControl;
    Rec3dProject = handles.Rec3dProject;
elseif nargout == 3
    Rec3dControl = handles.Rec3dControl;
    Rec3dProject = handles.Rec3dProject;
    Rec3dReconstruction = handles.Rec3dReconstruction;
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% getRec3dGUI
% -------------------------------------------------------------------------
function [Rec3dControl, Rec3dProject] = getRec3dGUI...
(handles, Rec3dControl, Rec3dProject)

if nargout == 1
    
    % Get Rec3dControl
    % ProtocolMode
    ProtocolMode = get(handles.ProtocolMode, 'Checked');
    Rec3dControl.ProtocolMode = ProtocolMode;
    % DemoMode
    DemoMode = get(handles.DemoMode, 'Checked');
    Rec3dControl.DemoMode = DemoMode;
    % DelMode
    DelMode = get(handles.DelMode, 'Checked');
    Rec3dControl.DelMode = DelMode;
    % ParallelMode
    ParallelMode = get(handles.ParallelMode, 'Checked');
    Rec3dControl.ParallelMode = ParallelMode;

end

if nargout == 2
    
    % Get Rec3dControl
    % ProtocolMode
    ProtocolMode = get(handles.ProtocolMode, 'Checked');
    Rec3dControl.ProtocolMode = ProtocolMode;
    % DemoMode
    DemoMode = get(handles.DemoMode, 'Checked');
    Rec3dControl.DemoMode = DemoMode;
    % DelMode
    DelMode = get(handles.DelMode, 'Checked');
    Rec3dControl.DelMode = DelMode;
    % ParallelMode
    ParallelMode = get(handles.ParallelMode, 'Checked');
    Rec3dControl.ParallelMode = ParallelMode;
    
    % Get HEAD
    ProjectName = get(handles.editProjectName, 'String');
    ProjectStatus = get(handles.dispProjectStatus, 'String');

    % Get ALIGNMENT
    % AlignmentMethod
    str = get(handles.popupAlignmentMethod, 'String');
    val = get(handles.popupAlignmentMethod, 'Value');
    AlignmentMethod = str{val};
    % ReferenceMarker
    str = get(handles.popupReferenceMarker, 'String');
    val = get(handles.popupReferenceMarker, 'Value');
    ReferenceMarker = str2num(str(val,:));
    
    % Get RECONSTRUCTION / Volume
    SizeX = str2num(get(handles.editSizeX, 'String'));
    SizeY = str2num(get(handles.editSizeY, 'String'));
    SizeZ = str2num(get(handles.editSizeZ, 'String'));
    PreBinning = get(handles.popupPreBinning, 'Value')-1;
    PostBinning = get(handles.popupPostBinning, 'Value')-1;

    % getRec3dGUI -> Rec3dProject
    % HEAD
    Rec3dProject.ProjectName   = ProjectName;
    Rec3dProject.ProjectStatus = ProjectStatus;
    % ALIGNMENT
    Rec3dProject.ALIGNMENT.AlignmentMethod = AlignmentMethod;
    Rec3dProject.ALIGNMENT.ReferenceMarker = ReferenceMarker;
    % RECONSTRUCTION / Volume
    Rec3dProject.VOLUME.SizeX = SizeX;
    Rec3dProject.VOLUME.SizeY = SizeY;
    Rec3dProject.VOLUME.SizeZ = SizeZ;
    Rec3dProject.VOLUME.PreBinning  = PreBinning;
    Rec3dProject.VOLUME.PostBinning = PostBinning;
    
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% setRec3dGUI
% -------------------------------------------------------------------------
function setRec3dGUI(handles, Rec3dControl, Rec3dProject)

% Return if Rec3dProject is empty
% Use GUI initialisation from guide
if isempty(Rec3dProject) == 1
    return;
end

% Rec3dControl -> setRec3dGUI
WaitbarMode = Rec3dControl.WaitbarMode;
ProtocolMode = Rec3dControl.ProtocolMode;
DemoMode = Rec3dControl.DemoMode;
DelMode = Rec3dControl.DelMode;
ParallelMode = Rec3dControl.ParallelMode;
DistributedComputingToolbox = Rec3dControl.DistributedComputingToolbox;
RecInMemory = Rec3dControl.RecInMemory;

% Rec3dProject -> setRec3dGUI
% HEAD
ProjectName = Rec3dProject.ProjectName;
ProjectStatus = Rec3dProject.ProjectStatus;
% PROJECT
TiltingGeometry = Rec3dProject.PROJECT.TiltingGeometry;
ProjectionsName1 = Rec3dProject.PROJECT.ProjectionsName1;
ProjectionsPath1 = Rec3dProject.PROJECT.ProjectionsPath1;
MarkerfileName1 = Rec3dProject.PROJECT.MarkerfileName1;
MarkerfilePath1 = Rec3dProject.PROJECT.MarkerfilePath1;
NumOfProj1 = Rec3dProject.PROJECT.NumOfProj1;
RefProj1 = Rec3dProject.PROJECT.RefProj1;
ProjectionsName2 = Rec3dProject.PROJECT.ProjectionsName2;
ProjectionsPath2 = Rec3dProject.PROJECT.ProjectionsPath2;
MarkerfileName2 = Rec3dProject.PROJECT.MarkerfileName2;
MarkerfilePath2 = Rec3dProject.PROJECT.MarkerfilePath2;
NumOfProj2 = Rec3dProject.PROJECT.NumOfProj2;
RefProj2 = Rec3dProject.PROJECT.RefProj2;
NumOfMarkers = Rec3dProject.PROJECT.NumOfMarkers;
Imdim = Rec3dProject.PROJECT.Imdim;
% ALIGNMENT
AlignmentMethod = Rec3dProject.ALIGNMENT.AlignmentMethod;
ReferenceMarker = Rec3dProject.ALIGNMENT.ReferenceMarker;
% RECONSTRUCTION / Volume
SizeX = Rec3dProject.VOLUME.SizeX;
SizeY = Rec3dProject.VOLUME.SizeY;
SizeZ = Rec3dProject.VOLUME.SizeZ;
PreBinning = Rec3dProject.VOLUME.PreBinning;
PostBinning = Rec3dProject.VOLUME.PostBinning;
% RECONSTRUCTION / Method
METHOD = Rec3dProject.METHOD;
% RECONSTRUCTION / Filter
FILTER = Rec3dProject.FILTER;
% RECONSTRUCTION / Parallel
PARALLEL = Rec3dProject.PARALLEL;

% Set Menue

%%%%%% Aktualize menue activation here

% Menue handles H
H = [handles.LoadProject;...
     handles.GenerateProject;...
     handles.ModifyProject;...
     handles.SaveProject;...
     handles.WaitbarMode;...
     handles.ProtocolMode;...
     handles.DemoMode;...
     handles.DelMode;...
     handles.ParallelMode;...
     handles.ShowProjections;...
     handles.ShowProject;...
     handles.CopyProject;...
     handles.Quit;...
     handles.DoAlignment;...
     handles.ClearAlignment;...
     handles.ShowAlignment;...
     handles.DoWarpAlignment;...
     handles.SetMethod;...
     handles.SetFilter;...
     handles.SetDetail;...
     handles.SetParallel;...
     handles.DoReconstruction;...
     handles.ClearReconstruction;...
     handles.LoadReconstruction;...
     handles.SaveReconstruction;...
     handles.ShowReconstruction;...
     handles.CreateTempFiles];
 
% Menue definition A 
%           ( Empty loaded aligned recYES recNO) 
A = {     'Enable'  'on'    'on'    'on'    'on'    'on';... %LoadProject                           
             'Enable'  'on'    'on'    'on'    'on'    'on';... %GenerateProject
             'Enable'  'on'    'on'    'off'   'off'   'off';...%ModifyProject
             'Enable'  'on'    'on'    'on'    'on'    'on';... %SaveProject
             'Enable'  'on'    'on'    'on'    'on'    'on';... %WaitbarMode
             'Enable'  'on'    'on'    'on'    'on'    'on';... %ProtocolMode
             'Enable'  'on'    'on'    'on'    'on'    'on';... %DemoMode
             'Enable'  'on'    'on'    'on'    'on'    'on';... %DelMode
             'Enable'  'on'    'on'    'on'    'on'    'on';... %ParallelMode
             'Enable'  'on'    'on'    'on'    'on'    'on';... %ShowProjections
             'Enable'  'on'    'on'    'on'    'on'    'on';... %ShowProject
             'Enable'  'on'    'on'    'on'    'on'    'on';... %CopyProject
             'Enable'  'on'    'on'    'on'    'on'    'on';... %Quit
             'Enable'  'on'    'on'    'off'   'off'   'off';...%DoAlignment
             'Enable'  'on'    'off'   'on'    'off'   'on';... %ClearAlignment
             'Enable'  'on'    'off'   'on'    'on'    'on';... %ShowAlignment
             'Enable'  'on'    'off'   'on'   'off'   'off';... %DoWarpAlignment
             'Enable'  'on'    'off'   'on'    'off'   'off';... %SetMethod
             'Enable'  'on'    'off'   'on'    'off'   'off';... %SetFilter
             'Enable'  'on'    'off'   'on'    'on'    'off';... %SetDetail
             'Enable'  'on'    'off'   'on'    'off'   'off';... %SetParallel
             'Enable'  'on'    'off'   'on'    'off'   'off';... %DoReconstruction
             'Enable'  'on'    'off'   'off'   'on'    'on';... %ClearReconstruction
             'Enable'  'on'    'off'   'off'   'off'   'on';... %LoadReconstruction
             'Enable'  'on'    'off'   'off'   'on'    'off';... %SaveReconstruction
             'Enable'  'on'    'off'   'off'   'on'    'off';... %ShowReconstruction
             'Enable'  'on'    'off'   'on'    'on'    'on'};   %CreateTempFiles

% Set menue activation
if strcmp(ProjectStatus, 'loaded') == 1
    for k=1:size(H,1)
         set(H(k,1), A(k,1), A(k,3));
    end
elseif strcmp(ProjectStatus, 'aligned') == 1
    for k=1:size(H,1)
         set(H(k,1), A(k,1), A(k,4));
    end
elseif strcmp(ProjectStatus, 'reconstructed') == 1 && ...
       strcmp(RecInMemory, 'yes') == 1
    for k=1:size(H,1)
         set(H(k,1), A(k,1), A(k,5));
    end
elseif strcmp(ProjectStatus, 'reconstructed') == 1 && ...
       strcmp(RecInMemory, 'no') == 1
    for k=1:size(H,1)
         set(H(k,1), A(k,1), A(k,6));
    end
end

% Set menue activation (special case)
if strcmp(ParallelMode, 'off') == 1
    set(handles.SetParallel, 'Enable', 'off');
end

% Set menue checkboxes
% WaitbarMode
if strcmp(WaitbarMode, 'on')
    set(handles.WaitbarMode, 'Checked', 'on');
else
    set(handles.WaitbarMode, 'Checked', 'off');
end
% ProtocolMode
if strcmp(ProtocolMode, 'on')
    set(handles.ProtocolMode, 'Checked', 'on');
else
    set(handles.ProtocolMode, 'Checked', 'off');
end
% DemoMode
if strcmp(DemoMode, 'on')
    set(handles.DemoMode, 'Checked', 'on');
else
    set(handles.DemoMode, 'Checked', 'off');
end
% DelMode
if strcmp(DelMode, 'on')
    set(handles.DelMode, 'Checked', 'on');
else
    set(handles.DelMode, 'Checked', 'off');
end
% ParallelMode
if strcmp(ParallelMode, 'on')
    set(handles.ParallelMode, 'Checked', 'on');
else
    set(handles.ParallelMode, 'Checked', 'off');
end

% Set menue checkboxes (special case)
if strcmp(DistributedComputingToolbox, 'no') == 1
    set(handles.ParallelMode, 'Checked', 'off');
    set(handles.ParallelMode, 'Enable', 'off');
end

% Reset HEAD
set(handles.editProjectName, 'String', '');
set(handles.dispProjectStatus, 'String', '');

% Reset PROJECT
set(handles.dispTiltingGeometry, 'String', '');
set(handles.dispProjectionsName1, 'String', '');
set(handles.dispProjectionsPath1, 'String', '');
set(handles.dispMarkerfileName1, 'String', '');
set(handles.dispMarkerfilePath1, 'String', '');
set(handles.dispNumOfProj1, 'String', '');
set(handles.dispRefProj1, 'String', '');
set(handles.dispProjectionsName2, 'String', '');
set(handles.dispProjectionsPath2, 'String', '');
set(handles.dispMarkerfileName2, 'String', '');
set(handles.dispMarkerfilePath2, 'String', '');
set(handles.dispNumOfProj2, 'String', '');
set(handles.dispRefProj2, 'String', '');
set(handles.dispNumOfMarkers, 'String', '');
set(handles.dispImdim, 'String', '');

% Set HEAD
set(handles.editProjectName, 'String', ProjectName);
set(handles.dispProjectStatus, 'String', ProjectStatus);

% Set PROJECT
set(handles.dispTiltingGeometry, 'String', TiltingGeometry);
set(handles.dispProjectionsName1, 'String', ProjectionsName1);
set(handles.dispProjectionsPath1, 'String', ProjectionsPath1);
set(handles.dispMarkerfileName1, 'String', MarkerfileName1);
set(handles.dispMarkerfilePath1, 'String', MarkerfilePath1);
set(handles.dispNumOfProj1, 'String', num2str(NumOfProj1));
set(handles.dispRefProj1, 'String', num2str(RefProj1));

if strcmp(TiltingGeometry, 'dualaxis') == 1
    
    set(handles.dispProjectionsName2, 'String', ProjectionsName2);
    set(handles.dispProjectionsPath2, 'String', ProjectionsPath2);
    set(handles.dispMarkerfileName2, 'String', MarkerfileName2);
    set(handles.dispMarkerfilePath2, 'String', MarkerfilePath2);
    set(handles.dispNumOfProj2, 'String', num2str(NumOfProj2));
    set(handles.dispRefProj2, 'String', num2str(RefProj2));

end

set(handles.dispNumOfMarkers, 'String', num2str(NumOfMarkers));
set(handles.dispImdim, 'String', num2str(Imdim));

% Set ALIGNMENT
% AlignmentMethod
if strcmp(AlignmentMethod, 'rigidbody') == 1
    set(handles.popupAlignmentMethod, 'Value', 1);
elseif strcmp(AlignmentMethod, 'freetilt') == 1
    set(handles.popupAlignmentMethod, 'Value', 2);
end
% ReferenceMarker
for k=1:(NumOfMarkers)
    MarkerIndex(k) = k;
end
set(handles.popupReferenceMarker, 'String', MarkerIndex, ...
                                  'Value', ReferenceMarker);

% Set RECONSTRUCTION / Volume
% Size
set(handles.editSizeX, 'String', num2str(SizeX));
set(handles.editSizeY, 'String', num2str(SizeY));
set(handles.editSizeZ, 'String', num2str(SizeZ));
% PreBinning
if PreBinning == 0
    set(handles.popupPreBinning, 'Value', 1);
elseif PreBinning == 1
    set(handles.popupPreBinning, 'Value', 2);
elseif PreBinning == 2
    set(handles.popupPreBinning, 'Value', 3);   
elseif PreBinning == 3
    set(handles.popupPreBinning, 'Value', 4);
elseif PreBinning == 4
    set(handles.popupPreBinning, 'Value', 5);    
elseif PreBinning == 5
    set(handles.popupPreBinning, 'Value', 6);
elseif PreBinning == 6
    set(handles.popupPreBinning, 'Value', 7);    
elseif PreBinning == 7
    set(handles.popupPreBinning, 'Value', 8);
end
% PostBinning
if PostBinning == 0
    set(handles.popupPostBinning, 'Value', 1);
elseif PostBinning == 1
    set(handles.popupPostBinning, 'Value', 2);
elseif PostBinning == 2
    set(handles.popupPostBinning, 'Value', 3);   
elseif PostBinning == 3
    set(handles.popupPostBinning, 'Value', 4);
elseif PostBinning == 4
    set(handles.popupPostBinning, 'Value', 5);    
elseif PostBinning == 5
    set(handles.popupPostBinning, 'Value', 6);
elseif PostBinning == 6
    set(handles.popupPostBinning, 'Value', 7);    
elseif PostBinning == 7
    set(handles.popupPostBinning, 'Value', 8);
end

% Set RECONSTRUCTION / Method
dispMETHODpre = METHOD;

dispMETHOD.Method = dispMETHODpre.ReconstructionMethod;
dispMETHOD.Normalization = dispMETHODpre.Normalization;
dispMETHOD.SmoothBorders = dispMETHODpre.SmoothBorders;
dispMETHOD.AlignTiltaxis = dispMETHODpre.AlignTiltaxis;
dispMETHOD.Handedness = dispMETHODpre.Handedness;
dispMETHOD.ApplyWeighting = dispMETHODpre.ApplyWeighting;
dispMETHOD.WeightingMethod = dispMETHODpre.WeightingMethod;
dispMETHOD.ObjectThickness = dispMETHODpre.ObjectThickness;
dispMETHOD.Taper = dispMETHODpre.Taper;
dispMETHOD.Iterations = dispMETHODpre.Iterations;
dispMETHOD.Relaxation = dispMETHODpre.Relaxation;
dispMETHOD.Pathlength = dispMETHODpre.Pathlength;

MethodFieldsAndValues = tom_Rec3dStructure2Cell(dispMETHOD);
set(handles.listboxMethod, 'String', MethodFieldsAndValues);

% Set RECONSTRUCTION / Filter
dispFILTERpre = FILTER;

dispFILTER.Type = dispFILTERpre.Type;
dispFILTER.Apply = dispFILTERpre.Apply;
dispFILTER.Times = dispFILTERpre.Value.times;
dispFILTER.Low = dispFILTERpre.Value.low;
dispFILTER.High = dispFILTERpre.Value.high;
dispFILTER.Smooth = dispFILTERpre.Value.smooth;
dispFILTER.Space = dispFILTERpre.Value.space;
dispFILTER.Method = dispFILTERpre.Value.method;
dispFILTER.Radius = dispFILTERpre.Value.radius;

FilterFieldsAndValues = tom_Rec3dStructure2Cell(dispFILTER);
set(handles.listboxFilter, 'String', FilterFieldsAndValues);

% Clear RECONSTRUCTION / Parallel
set(handles.listboxParallel, 'String', '');

% Set RECONSTRUCTION / Parallel
if strcmp(ParallelMode, 'on')
    
    ParallelFieldsAndValues = tom_Rec3dStructure2Cell(PARALLEL);
    set(handles.listboxParallel, 'String', ParallelFieldsAndValues);

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Rec3d2handles
% -------------------------------------------------------------------------
function [handles] = Rec3d2handles...
(handles, Rec3dControl, Rec3dProject, Rec3dReconstruction)

% Rec3d -> handles
if nargin == 2
    handles.Rec3dControl = Rec3dControl;
elseif nargin == 3
    handles.Rec3dControl = Rec3dControl;
    handles.Rec3dProject = Rec3dProject;
elseif nargin == 4
    handles.Rec3dControl = Rec3dControl;
    handles.Rec3dProject = Rec3dProject;
    handles.Rec3dReconstruction = Rec3dReconstruction;
end
% -------------------------------------------------------------------------










% -------------------------------------------------------------------------
% Rec3dGenerateProject
% -------------------------------------------------------------------------
function [Rec3dProject] = Rec3dGenerateProject()

% Which TiltingGeometry
TiltingGeometry = questdlg('Please choose Tilting Geometry of your project:', ...
                           'Generate Project', 'singleaxis', 'dualaxis', 'Cancel', 'singleaxis');

if strcmp(TiltingGeometry, 'Cancel') == 1 || strcmp(TiltingGeometry, '') == 1
    Rec3dProject = [];
    return;
end

% Locate projections and markerfile of your 'singleaxis' project
% and generate this project
if strcmp(TiltingGeometry, 'singleaxis') == 1

    % Get ProjectionsName1long and ProjectionsPath1
    [ProjectionsName1long, ProjectionsPath1] = uigetfile...
        ('*.em', 'Select arbitrary projection of your singleaxis tiltseries');
    
    if isequal(ProjectionsName1long, 0) || isequal(ProjectionsPath1, 0)
        Rec3dProject = [];
        return;
    end

    % Get ProjectionsExt1
    [pathstr, name, ProjectionsExt1] = fileparts(ProjectionsName1long);
    
    % Find underline index and seperate ProjectionsName1
    UnderlineIndex = findstr(ProjectionsName1long, '_');
    if isempty(UnderlineIndex) == 1
        % ERROR
        msgbox('Projections Name is not valid!', ...
                      'Generate Project', 'error');
        Rec3dProject = [];
        return;
    end
    WhichUnderline = size(UnderlineIndex, 2);
    UnderlineIndex = UnderlineIndex(1,WhichUnderline);
    
    for k = 1:(UnderlineIndex)
        ProjectionsName1(k) = ProjectionsName1long(k);
    end
    
    % Get MarkerfileName1 and MarkerfilePath1
    [MarkerfileName1, MarkerfilePath1] = uigetfile...
        ('*.em', 'Select markerfile of your singleaxis tiltseries');

    if isequal(MarkerfileName1, 0) || isequal(MarkerfilePath1, 0)
        Rec3dProject = [];
        return;
    end

    % Get MarkerfileExt1
    [pathstr, name, MarkerfileExt1] = fileparts(MarkerfileName1);

    % Which 'singleaxis' ProjectName
    ProjectName = inputdlg('Please choose Project Name of your singleaxis project:', ...
                                            'Generate Project');

    if isempty(ProjectName) == 1
        Rec3dProject = [];
        return;
    end

    ProjectName = ProjectName{1};

    % GenerateSingleaxisProject
    Rec3dProject = GenerateProject(ProjectName, TiltingGeometry, ...
        ProjectionsName1, ProjectionsPath1, ProjectionsExt1, ...
        MarkerfileName1, MarkerfilePath1, MarkerfileExt1, ...
        '', '', '', '', '', '');

end

% Locate projections and markerfiles of your 'dualaxis' project
% and generate this project
if strcmp(TiltingGeometry, 'dualaxis') == 1
    
    % Get ProjectionsName1long and ProjectionsPath1
    [ProjectionsName1long, ProjectionsPath1] = uigetfile...
        ('*.em', 'Select arbitrary projection of your FIRST dualaxis tiltseries');

    if isequal(ProjectionsName1long, 0) || isequal(ProjectionsPath1, 0)
        Rec3dProject = [];
        return;
    end

    % Get ProjectionsExt1
    [pathstr, name, ProjectionsExt1] = fileparts(ProjectionsName1long);

    % Find underline index and seperate ProjectionsName1
    UnderlineIndex = findstr(ProjectionsName1long, '_');
    if isempty(UnderlineIndex) == 1
        % ERROR
        msgbox('Projections Name is not valid!', ...
                      'Generate Project', 'error');
        Rec3dProject = [];
        return;
    end
    WhichUnderline = size(UnderlineIndex, 2);
    UnderlineIndex = UnderlineIndex(1,WhichUnderline);
    
    for k = 1:(UnderlineIndex)
        ProjectionsName1(k) = ProjectionsName1long(k);
    end
    
    % Get MarkerfileName1 and MarkerfilePath1
    [MarkerfileName1, MarkerfilePath1] = uigetfile...
        ('*.em', 'Select markerfile of your FIRST dualaxis tiltseries');

    if isequal(MarkerfileName1, 0) || isequal(MarkerfilePath1, 0)
        Rec3dProject = [];
        return;
    end

    % Get MarkerfileExt1
    [pathstr, name, MarkerfileExt1] = fileparts(MarkerfileName1);

    % Get ProjectionsName2long and ProjectionsPath2
    [ProjectionsName2long, ProjectionsPath2] = uigetfile...
        ('*.em', 'Select arbitrary projection of your SECOND dualaxis tiltseries');

    if isequal(ProjectionsName2long, 0) || isequal(ProjectionsPath2, 0)
        Rec3dProject = [];
        return;
    end

    % Get ProjectionsExt2
    [pathstr, name, ProjectionsExt2] = fileparts(ProjectionsName2long);

    % Find underline index and seperate ProjectionsName2
    UnderlineIndex = findstr(ProjectionsName2long, '_');
    if isempty(UnderlineIndex) == 1
        % ERROR
        msgbox('Projections Name is not valid!', ...
                      'Generate Project', 'error');
        Rec3dProject = [];
        return;
    end
    WhichUnderline = size(UnderlineIndex, 2);
    UnderlineIndex = UnderlineIndex(1,WhichUnderline);
    
    for k = 1:(UnderlineIndex)
        ProjectionsName2(k) = ProjectionsName2long(k);
    end
    
    % Get MarkerfileName2 and MarkerfilePath2
    [MarkerfileName2, MarkerfilePath2] = uigetfile...
        ('*.em', 'Select markerfile of your SECOND dualaxis tiltseries');

    if isequal(MarkerfileName2, 0) || isequal(MarkerfilePath2, 0)
        Rec3dProject = [];
        return;
    end

    % Get MarkerfileExt2
    [pathstr, name, MarkerfileExt2] = fileparts(MarkerfileName2);

    % Which 'dualaxis' ProjectName
    ProjectName = inputdlg('Please choose Project Name of your dualaxis project:', ...
                                            'Generate Project');

    if isempty(ProjectName) == 1
        Rec3dProject = [];
        return;
    end

    ProjectName = ProjectName{1};

    % GenerateDualaxisProject
    Rec3dProject = GenerateProject(ProjectName, TiltingGeometry, ...
        ProjectionsName1, ProjectionsPath1, ProjectionsExt1, ...
        MarkerfileName1, MarkerfilePath1, MarkerfileExt1, ...
        ProjectionsName2, ProjectionsPath2, ProjectionsExt2, ...
        MarkerfileName2, MarkerfilePath2, MarkerfileExt2);

end

function [Rec3dProject] = GenerateProject(ProjectName, TiltingGeometry, ...
ProjectionsName1, ProjectionsPath1, ProjectionsExt1, ...
MarkerfileName1, MarkerfilePath1, MarkerfileExt1, ...
ProjectionsName2, ProjectionsPath2, ProjectionsExt2, ...
MarkerfileName2, MarkerfilePath2, MarkerfileExt2)

% Generate 'singleaxis' project
if strcmp(TiltingGeometry, 'singleaxis') == 1

% Load Projection1_1
workingdir = pwd;
cd (ProjectionsPath1);

%     % Check is EM-file
%     CheckProjection1_1 = tom_isemfile...
%     ([ProjectionsName1 '1' ProjectionsExt1]);
% 
%     if CheckProjection1_1 == 0 || CheckProjection1_1 == -1
%     % ERROR
%     msgbox('Projections are not in EM-format!', ...
%                    'Generate Project', 'error');
%     Rec3dProject = [];
%     return;
%     end

Projection1 = tom_emread([ProjectionsName1 '1' ProjectionsExt1]);
cd (workingdir);

% % Check is projection
% if size(Projection1.Value, 3) ~= 1
%     % ERROR
%     msgbox('Projections are not projections!', ...
%                   'Generate Project', 'error');
%     Rec3dProject = [];
%     return;
% end

% Load Markerfile1
workingdir = pwd;
cd (MarkerfilePath1);

    % Check is EM-file
    CheckMarkerfile1 = tom_isemfile(MarkerfileName1);

    if CheckMarkerfile1 == 0 || CheckMarkerfile1 == -1
    % ERROR
    msgbox('Markerfile is not in EM-format!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
    end

Markerfile1 = tom_emread(MarkerfileName1);
cd (workingdir);

% Check is markerfile
if size(Markerfile1.Value, 3) == 1
    % ERROR
    msgbox('Markerfile is not a markerfile!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
end

% Number of projections of tiltseries1
NumOfProj1 = size(Markerfile1.Value, 2);

% Check number of projections
if NumOfProj1 > 999
    % ERROR
    msgbox('Number of projections is limited to max 999!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
end

% ProjectionsName and ProjectionsPath of Tiltseries1
for k=1:NumOfProj1
    if k<10
    ProjectionsName{k} = [ProjectionsName1 '00' num2str(k) ProjectionsExt1];
    ProjectionsPath{k} = ProjectionsPath1;
    end
    if k>= 10 && k<100
    ProjectionsName{k} = [ProjectionsName1 '0' num2str(k) ProjectionsExt1];
    ProjectionsPath{k} = ProjectionsPath1;
    end
    if k>= 100 && k<1000
    ProjectionsName{k} = [ProjectionsName1 num2str(k) ProjectionsExt1];
    ProjectionsPath{k} = ProjectionsPath1;
    end
end

% Number of markers
NumOfMarkers = size(Markerfile1.Value, 3);

% Measured tiltangles of tiltseries1
Tiltangles1(1:NumOfProj1) = Markerfile1.Value(1, 1:NumOfProj1, 1);
% Reference projection (minimum tilt) of tiltseries1
[notneeded, RefProj1] = min(abs(Tiltangles1));

% Projection dimension
Imdim = Projection1.Header.Size(1);

% Check projections
if Imdim ~= Projection1.Header.Size(2) 
    % ERROR
    msgbox('Projections are not quadratic!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
end

% Measured x coordinates of tiltseries1
x1(1:NumOfProj1, 1:NumOfMarkers) = ...
Markerfile1.Value(2, 1:NumOfProj1, 1:NumOfMarkers);
x1 = x1';
% Measured y coordinates of tiltseries1
y1(1:NumOfProj1, 1:NumOfMarkers) = ...
Markerfile1.Value(3, 1:NumOfProj1, 1:NumOfMarkers);
y1 = y1';
MarkersOnProj1(:,:,1) = x1;
MarkersOnProj1(:,:,2) = y1;

% Generate 'singleaxis' project structure
Rec3dProject = GenerateProjectStructure...
(ProjectName, 'loaded', TiltingGeometry, ...
 ProjectionsName1, ProjectionsPath1, ProjectionsExt1, ...
 MarkerfileName1, MarkerfilePath1, MarkerfileExt1, ...
 NumOfProj1, RefProj1, Tiltangles1, MarkersOnProj1, Markerfile1, ...
 '', '', '', ...
 '', '', '', ...
 0, 0, 0, 0, 0, ...
 NumOfMarkers, Imdim, ProjectionsName, ProjectionsPath);

end

% Generate 'dualaxis' project
if strcmp(TiltingGeometry, 'dualaxis') == 1

% Load Projection1_1
workingdir = pwd;
cd (ProjectionsPath1);

    % Check is EM-file
    CheckProjection1_1 = tom_isemfile...
    ([ProjectionsName1 '1' ProjectionsExt1]);

    if CheckProjection1_1 == 0 || CheckProjection1_1 == -1
    % ERROR
    msgbox('Projections are not in EM-format!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
    end

Projection1 = tom_emread([ProjectionsName1 '1' ProjectionsExt1]);
cd (workingdir);

% Check is projection
if size(Projection1.Value, 3) ~= 1
    % ERROR
    msgbox('Projections are not projections!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
end

% Load Markerfile1
workingdir = pwd;
cd (MarkerfilePath1);

    % Check is EM-file
    CheckMarkerfile1 = tom_isemfile(MarkerfileName1);

    if CheckMarkerfile1 == 0 || CheckMarkerfile1 == -1
    % ERROR
    msgbox('Markerfile is not in EM-format!', ...
                   'Generate Project', 'error');
    Rec3dProject = [];
    return;
    end

Markerfile1 = tom_emread(MarkerfileName1);
cd (workingdir);

% Check is markerfile
if size(Markerfile1.Value, 3) == 1
    % ERROR
    msgbox('Markerfile is not a markerfile!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
end

% Load Projection2_1
workingdir = pwd;
cd (ProjectionsPath2);

    % Check is EM-file
    CheckProjection2_1 = tom_isemfile...
    ([ProjectionsName2 '1' ProjectionsExt2]);

    if CheckProjection2_1 == 0 || CheckProjection2_1 == -1
    % ERROR
    msgbox('Projections are not in EM-format!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
    end

Projection2 = tom_emread([ProjectionsName2 '1' ProjectionsExt2]);
cd (workingdir);

% Check is projection
if size(Projection2.Value, 3) ~= 1
    % ERROR
    msgbox('Projections are not projections!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
end

% Load Markerfile2
workingdir = pwd;
cd (MarkerfilePath2);

    % Check is EM-file
    CheckMarkerfile2 = tom_isemfile(MarkerfileName2);

    if CheckMarkerfile2 == 0 || CheckMarkerfile2 == -1
    % ERROR
    msgbox('Markerfile is not in EM-format!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
    end

Markerfile2 = tom_emread(MarkerfileName2);
cd (workingdir);

% Check is markerfile
if size(Markerfile2.Value, 3) == 1
    % ERROR
    msgbox('Markerfile is not a markerfile!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
end

% Number of projections of tiltseries1
NumOfProj1 = size(Markerfile1.Value, 2);

% Check number of projections
if NumOfProj1 > 999
    % ERROR
    msgbox('Number of projections is limited to max 999!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
end

% ProjectionsName and ProjectionsPath of Tiltseries1
for k=1:NumOfProj1
    if k<10
    ProjectionsName{k} = [ProjectionsName1 '00' num2str(k) ProjectionsExt1];
    ProjectionsPath{k} = ProjectionsPath1;
    end
    if k>= 10 && k<100
    ProjectionsName{k} = [ProjectionsName1 '0' num2str(k) ProjectionsExt1];
    ProjectionsPath{k} = ProjectionsPath1;
    end
    if k>= 100 && k<1000
    ProjectionsName{k} = [ProjectionsName1 num2str(k) ProjectionsExt1];
    ProjectionsPath{k} = ProjectionsPath1;
    end
end

% Number of projections of tiltseries2
NumOfProj2 = size(Markerfile2.Value, 2);

% Check number of projections
if NumOfProj2 > 999
    % ERROR
    msgbox('Number of projections is limited to max 999!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
end

% ProjectionsName and ProjectionsPath of Tiltseries2
for k=1:NumOfProj2
    if k<10
    ProjectionsName{k+NumOfProj1} = [ProjectionsName2 '00' num2str(k) ProjectionsExt2];
    ProjectionsPath{k+NumOfProj1} = ProjectionsPath2;
    end
    if k>= 10 && k<100
    ProjectionsName{k+NumOfProj1} = [ProjectionsName2 '0' num2str(k) ProjectionsExt2];
    ProjectionsPath{k+NumOfProj1} = ProjectionsPath2;
    end
    if k>= 100 && k<1000
    ProjectionsName{k+NumOfProj1} = [ProjectionsName2 num2str(k) ProjectionsExt2];
    ProjectionsPath{k+NumOfProj1} = ProjectionsPath2;
    end
end

% Number of markers
NumOfMarkers1 = size(Markerfile1.Value, 3);
NumOfMarkers2 = size(Markerfile2.Value, 3);

% Check equal number of markers
if NumOfMarkers1 == NumOfMarkers2
    NumOfMarkers = NumOfMarkers1;
else
    % ERROR
    msgbox('Number of markers must be equal!', ...
                  'Generate Project', 'error');
    Rec3dProject = [];
    return;
end

% Measured tiltangles of tiltseries1
Tiltangles1(1:NumOfProj1) = Markerfile1.Value(1, 1:NumOfProj1, 1);
% Measured tiltangles of tiltseries2
Tiltangles2(1:NumOfProj2) = Markerfile2.Value(1, 1:NumOfProj2, 1);
% Reference projection (minimum tilt) of tiltseries1
[notneeded, RefProj1] = min(abs(Tiltangles1));
% Reference projection (minimum tilt) of tiltseries2
[notneeded, RefProj2] = min(abs(Tiltangles2));
% Projection dimension
Imdim1 = Projection1.Header.Size(1);
Imdim2 = Projection2.Header.Size(1);

% Check projections
if (Imdim1 == Projection1.Header.Size(2)) && ...
   (Imdim2 == Projection2.Header.Size(2))

    if Imdim1 == Imdim2
         Imdim = Imdim1;
    else
         % ERROR
         msgbox('Projections must have equal size!', ...
                       'Generate Project', 'error');
         Rec3dProject = [];
         return;
    end

else
    % ERROR
    msgbox('Projections are not quadratic!', ...
                   'Generate Project', 'error');
    Rec3dProject = [];
    return;
end

% Measured x coordinates of tiltseries1
x1(1:NumOfProj1, 1:NumOfMarkers) = ...
Markerfile1.Value(2, 1:NumOfProj1, 1:NumOfMarkers);
x1 = x1';
% Measured y coordinates of tiltseries1
y1(1:NumOfProj1, 1:NumOfMarkers) = ...
Markerfile1.Value(3, 1:NumOfProj1, 1:NumOfMarkers);
y1 = y1';
MarkersOnProj1(:,:,1) = x1;
MarkersOnProj1(:,:,2) = y1;

% Measured x coordinates of tiltseries2
x2(1:NumOfProj2, 1:NumOfMarkers) = ...
Markerfile2.Value(2, 1:NumOfProj2, 1:NumOfMarkers);
x2 = x2';
% Measured y coordinates of tiltseries2
y2(1:NumOfProj2, 1:NumOfMarkers) = ...
Markerfile2.Value(3, 1:NumOfProj2, 1:NumOfMarkers);
y2 = y2';
MarkersOnProj2(:,:,1) = x2;
MarkersOnProj2(:,:,2) = y2;

% Generate 'dualaxis' project structure
Rec3dProject = GenerateProjectStructure...
(ProjectName, 'loaded', TiltingGeometry, ...
 ProjectionsName1, ProjectionsPath1, ProjectionsExt1, ...
 MarkerfileName1, MarkerfilePath1, MarkerfileExt1, ...
 NumOfProj1, RefProj1, Tiltangles1, MarkersOnProj1, Markerfile1, ...
 ProjectionsName2, ProjectionsPath2, ProjectionsExt2, ...
 MarkerfileName2, MarkerfilePath2, MarkerfileExt2, ...
 NumOfProj2, RefProj2, Tiltangles2, MarkersOnProj2, Markerfile2, ...
 NumOfMarkers, Imdim, ProjectionsName, ProjectionsPath);

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% GenerateProjectStructure
% -------------------------------------------------------------------------
function [Rec3dProject] = GenerateProjectStructure...
(ProjectName, ProjectStatus, TiltingGeometry, ...
 ProjectionsName1, ProjectionsPath1, ProjectionsExt1, ...
 MarkerfileName1, MarkerfilePath1, MarkerfileExt1, ...
 NumOfProj1, RefProj1, Tiltangles1, MarkersOnProj1, Markerfile1, ...
 ProjectionsName2, ProjectionsPath2, ProjectionsExt2, ...
 MarkerfileName2, MarkerfilePath2, MarkerfileExt2, ...
 NumOfProj2, RefProj2, Tiltangles2, MarkersOnProj2, Markerfile2, ...
 NumOfMarkers, Imdim, ProjectionsName, ProjectionsPath)

%%%%%% Aktualize Rec3dProject here

% HEAD
Rec3dProject.ProjectName = ProjectName;
Rec3dProject.ProjectStatus = ProjectStatus;
Rec3dProject.ProjectVersion = '1.0';
% PROJECT
Rec3dProject.PROJECT.TiltingGeometry = TiltingGeometry;
Rec3dProject.PROJECT.ProjectionsName1 = ProjectionsName1;
Rec3dProject.PROJECT.ProjectionsPath1 = ProjectionsPath1;
Rec3dProject.PROJECT.ProjectionsExt1 = ProjectionsExt1;
Rec3dProject.PROJECT.MarkerfileName1 = MarkerfileName1;
Rec3dProject.PROJECT.MarkerfilePath1 = MarkerfilePath1;
Rec3dProject.PROJECT.MarkerfileExt1 = MarkerfileExt1;
Rec3dProject.PROJECT.NumOfProj1 = NumOfProj1;
Rec3dProject.PROJECT.NumOfProjOriginal1 = NumOfProj1;
Rec3dProject.PROJECT.RefProj1 = RefProj1;
Rec3dProject.PROJECT.Tiltangles1 = Tiltangles1;
Rec3dProject.PROJECT.MarkersOnProj1 = MarkersOnProj1;
Rec3dProject.PROJECT.MarkersOnProjOriginal1 = MarkersOnProj1;
Rec3dProject.PROJECT.Markerfile1 = Markerfile1;
Rec3dProject.PROJECT.MarkerfileOriginal1 = Markerfile1;
Rec3dProject.PROJECT.ProjectionsName2 = ProjectionsName2;
Rec3dProject.PROJECT.ProjectionsPath2 = ProjectionsPath2;
Rec3dProject.PROJECT.ProjectionsExt2 = ProjectionsExt2;
Rec3dProject.PROJECT.MarkerfileName2 = MarkerfileName2;
Rec3dProject.PROJECT.MarkerfilePath2 = MarkerfilePath2;
Rec3dProject.PROJECT.MarkerfileExt2 = MarkerfileExt2;
Rec3dProject.PROJECT.NumOfProj2 = NumOfProj2;
Rec3dProject.PROJECT.NumOfProjOriginal2 = NumOfProj2;
Rec3dProject.PROJECT.RefProj2 = RefProj2;
Rec3dProject.PROJECT.Tiltangles2 = Tiltangles2;
Rec3dProject.PROJECT.MarkersOnProj2 = MarkersOnProj2;
Rec3dProject.PROJECT.MarkersOnProjOriginal2 = MarkersOnProj2;
Rec3dProject.PROJECT.Markerfile2 = Markerfile2;
Rec3dProject.PROJECT.MarkerfileOriginal2 = Markerfile2;
Rec3dProject.PROJECT.NumOfProj = (NumOfProj1+NumOfProj2);
Rec3dProject.PROJECT.NumOfMarkers = NumOfMarkers;
Rec3dProject.PROJECT.Imdim = Imdim;
% ALIGNMENT
Rec3dProject.ALIGNMENT.AlignmentMethod = 'rigidbody';
Rec3dProject.ALIGNMENT.ReferenceMarker = 1;
Rec3dProject.ALIGNMENT.Origin1 = [0 0 0];
Rec3dProject.ALIGNMENT.m3d1 = 0;
Rec3dProject.ALIGNMENT.Tiltaxis1 = 0;
Rec3dProject.ALIGNMENT.tx1 = 0;
Rec3dProject.ALIGNMENT.ty1 = 0;
Rec3dProject.ALIGNMENT.isoscale1 = 1;
Rec3dProject.ALIGNMENT.WarpDone1 = 'no';
Rec3dProject.ALIGNMENT.WarpAlignment1 = 0;
Rec3dProject.ALIGNMENT.Origin2 = [0 0 0];
Rec3dProject.ALIGNMENT.m3d2 = 0;
Rec3dProject.ALIGNMENT.Tiltaxis2 = 0;
Rec3dProject.ALIGNMENT.tx2 = 0;
Rec3dProject.ALIGNMENT.ty2 = 0;
Rec3dProject.ALIGNMENT.isoscale2 = 1;
Rec3dProject.ALIGNMENT.WarpDone2 = 'no';
Rec3dProject.ALIGNMENT.WarpAlignment2 = 0;
Rec3dProject.ALIGNMENT.RotMatrix = 0;
Rec3dProject.ALIGNMENT.Psi = 0;
Rec3dProject.ALIGNMENT.Theta = 0;
Rec3dProject.ALIGNMENT.Phi = 0;
% ALIGNMENT / Residuals
Rec3dProject.ALGRESIDUALS.ResidualMatrix1 = 0;
Rec3dProject.ALGRESIDUALS.Sigma1 = 0;
Rec3dProject.ALGRESIDUALS.AveragePerProjection1 = 0;
Rec3dProject.ALGRESIDUALS.AveragePerMarker1 = 0;
Rec3dProject.ALGRESIDUALS.ResidualMatrix2 = 0;
Rec3dProject.ALGRESIDUALS.Sigma2 = 0;
Rec3dProject.ALGRESIDUALS.AveragePerProjection2 = 0;
Rec3dProject.ALGRESIDUALS.AveragePerMarker2 = 0;
Rec3dProject.ALGRESIDUALS.EulerAnglesResidual = 0;
Rec3dProject.ALGRESIDUALS.MaximumResidual = 0;
Rec3dProject.ALGRESIDUALS.ResidualSpheres = 0;
Rec3dProject.ALGRESIDUALS.AverageResidualSphere = 0;
% RECONSTRUCTION
Rec3dProject.NAMEPATHEXT.ReconstructionName = '';
Rec3dProject.NAMEPATHEXT.ReconstructionPath = '';
Rec3dProject.NAMEPATHEXT.ReconstructionExt = '';
Rec3dProject.NAMEPATHEXT.TempFilesName = '';
Rec3dProject.NAMEPATHEXT.TempFilesPath = '';
Rec3dProject.NAMEPATHEXT.TempFilesExt = '';
% RECONSTRUCTION / Parameter
Rec3dProject.PARAMETER.ProjectionsName = ProjectionsName;
Rec3dProject.PARAMETER.ProjectionsNameOriginal = ProjectionsName;
Rec3dProject.PARAMETER.ProjectionsPath = ProjectionsPath;
Rec3dProject.PARAMETER.ProjectionsPathOriginal = ProjectionsPath;
Rec3dProject.PARAMETER.Tiltangles = 0;
Rec3dProject.PARAMETER.Tiltaxis = 0;
Rec3dProject.PARAMETER.ProjDir = 0;
Rec3dProject.PARAMETER.tx = 0;
Rec3dProject.PARAMETER.ty = 0;
Rec3dProject.PARAMETER.isoscale = 1;
% RECONSTRUCTION / Volume
Rec3dProject.VOLUME.SizeX = Imdim;
Rec3dProject.VOLUME.SizeY = Imdim;
Rec3dProject.VOLUME.SizeZ = Imdim;
Rec3dProject.VOLUME.PreBinning = 0;
Rec3dProject.VOLUME.PostBinning = 0;
% RECONSTRUCTION / Method
Rec3dProject.METHOD.ReconstructionMethod = 'WBP';
Rec3dProject.METHOD.Normalization = 'phase';
Rec3dProject.METHOD.SmoothBorders = 0;
Rec3dProject.METHOD.AlignTiltaxis = 'ProjDir';
Rec3dProject.METHOD.Handedness = 0;
Rec3dProject.METHOD.ApplyWeighting = 'on';
Rec3dProject.METHOD.WeightingMethod = 'exact';
Rec3dProject.METHOD.ObjectThickness = Imdim;
Rec3dProject.METHOD.Taper = 'off';
Rec3dProject.METHOD.Iterations = 12;
Rec3dProject.METHOD.Relaxation = 0.05;
Rec3dProject.METHOD.Pathlength = 'ones';
% RECONSTRUCTION / Filter
Rec3dProject.FILTER.Apply = 1;
Rec3dProject.FILTER.Value.times = 1;
Rec3dProject.FILTER.Value.low = 0;
Rec3dProject.FILTER.Value.high = Imdim;
Rec3dProject.FILTER.Value.smooth = 0;
Rec3dProject.FILTER.Value.space = 'real';
Rec3dProject.FILTER.Value.method = 'quadr';
Rec3dProject.FILTER.Value.radius = 0;
Rec3dProject.FILTER.Type = 'bandpass';
% RECONSTRUCTION / Detail
Rec3dProject.DETAIL.DetailMode = 'off';
Rec3dProject.DETAIL.OverviewSizeZ = 0;
Rec3dProject.DETAIL.OverviewPreBinning = 0;
Rec3dProject.DETAIL.OverviewPostBinning = 0;
Rec3dProject.DETAIL.NumberOfDetails = 0;
Rec3dProject.DETAIL.DetailCoordinates = [];
Rec3dProject.DETAIL.av3_alignstruct = [];
% RECONSTRUCTION / Parallel
Rec3dProject.PARALLEL.jobmanager = 'default_jobmanager';
Rec3dProject.PARALLEL.packageloss = 0;
Rec3dProject.PARALLEL.number_of_tasks = 1;
Rec3dProject.PARALLEL.workers.min =  1;
Rec3dProject.PARALLEL.workers.max =  16;
Rec3dProject.PARALLEL.timeout = 3600;
Rec3dProject.PARALLEL.restart_workers = 0;
% RECONSTRUCTION / Sirtresiduals
Rec3dProject.SIRTRESIDUALS.EuclidianDistancePerProjection = 0;
Rec3dProject.SIRTRESIDUALS.EuclidianDistancePerIteration = 0;
Rec3dProject.SIRTRESIDUALS.dProjVal = 0;
Rec3dProject.SIRTRESIDUALS.DifVolVal = 0;
Rec3dProject.SIRTRESIDUALS.RecVolVal = 0;
% -------------------------------------------------------------------------

