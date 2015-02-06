function varargout = tom_Rec3dModifyProjectGUI(varargin)
%TOM_REC3DMODIFYPROJECTGUI is a module of TOM_REC3DGUI.
%
%   Rec3dProject = tom_Rec3dModifyProjectGUI(Rec3dProject);
%
%   It allows to exclude projections from the tiltseries and
%   from reconstruction.
%
%PARAMETERS
%
%  INPUT
%   Rec3dProject     main-structure of Rec3d
%  
%  OUTPUT
%   Rec3dProject     aktualized main-structure of Rec3d
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
                   'gui_OpeningFcn', @tom_Rec3dModifyProjectGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_Rec3dModifyProjectGUI_OutputFcn, ...
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
function tom_Rec3dModifyProjectGUI_OpeningFcn(hObject, eventdata, handles, varargin)

if (size(varargin, 2) < 1)
    
    % Input REC3DPROJECT
    handles.inputREC3DPROJECT = [];
    
else
    
    % Input REC3DPROJECT
    handles.inputREC3DPROJECT = varargin{1};
    
end

% handles.inputREC3DPROJECT -> handles.outputREC3DPROJECT 
handles.outputREC3DPROJECT = handles.inputREC3DPROJECT;

% Display REC3DPROJECT
setRec3dModifyProjectGUIValues(handles);

% Update handles structure
guidata(hObject, handles);
uiwait;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% OutputFunction
% -------------------------------------------------------------------------
function varargout = tom_Rec3dModifyProjectGUI_OutputFcn(hObject, eventdata, handles)

% Output METHOD
varargout{1} = handles.outputREC3DPROJECT;

% Delete Rec3dModifyProjectGUI
delete(handles.Rec3dModifyProjectGUI);
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

% Aktualize outputREC3DPROJECT
handles = getRec3dModifyProjectGUIValues(handles);

% Update handles structure
guidata(hObject, handles);
uiresume;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% setRec3dModifyProjectGUIValues
% -------------------------------------------------------------------------
function setRec3dModifyProjectGUIValues(handles);

if isempty(handles.inputREC3DPROJECT) == 1
    set(handles.pushbuttonT1toS1, 'Enable', 'off');
    set(handles.pushbuttonDeleteS1, 'Enable', 'off');
    set(handles.pushbuttonT2toS2, 'Enable', 'off');
    set(handles.pushbuttonDeleteS2, 'Enable', 'off');
    return;
end

% Simplify
TiltingGeometry = handles.inputREC3DPROJECT.PROJECT.TiltingGeometry;
NumOfProj1 = handles.inputREC3DPROJECT.PROJECT.NumOfProj1;
NumOfProjOriginal1 = handles.inputREC3DPROJECT.PROJECT.NumOfProjOriginal1;
NumOfProj2 = handles.inputREC3DPROJECT.PROJECT.NumOfProj2;
NumOfProjOriginal2 = handles.inputREC3DPROJECT.PROJECT.NumOfProjOriginal2;
ProjectionsName = handles.inputREC3DPROJECT.PARAMETER.ProjectionsName;
ProjectionsNameOriginal = handles.inputREC3DPROJECT.PARAMETER.ProjectionsNameOriginal;

% Background color listbox
set(handles.listboxTiltseries1, 'BackgroundColor', 'White');
set(handles.listboxSelection1, 'BackgroundColor', 'White');
set(handles.listboxTiltseries2, 'BackgroundColor', 'White');
set(handles.listboxSelection2, 'BackgroundColor', 'White');

% Populate listbox tiltseries1
for k=1:NumOfProjOriginal1
    ProjectionsNameOriginal1{k} = ProjectionsNameOriginal{k};
end
    set(handles.listboxTiltseries1, 'String', ProjectionsNameOriginal1);
% Populate listbox tiltseries2
if strcmp(TiltingGeometry, 'dualaxis')
for k=1:NumOfProjOriginal2
    ProjectionsNameOriginal2{k} = ProjectionsNameOriginal{k+NumOfProjOriginal1};
end
    set(handles.listboxTiltseries2, 'String', ProjectionsNameOriginal2);
else
    set(handles.listboxTiltseries2, 'Enable', 'off');
    set(handles.pushbuttonT2toS2, 'Enable', 'off');
end

% Populate listbox selection1
for k=1:NumOfProj1
    ProjectionsName1{k} = ProjectionsName{k};
end
    set(handles.listboxSelection1, 'String', ProjectionsName1);
% Populate listbox selection2
if strcmp(TiltingGeometry, 'dualaxis')
for k=1:NumOfProj2
    ProjectionsName2{k} = ProjectionsName{k+NumOfProj1};
end
    set(handles.listboxSelection2, 'String', ProjectionsName2);
else
    set(handles.listboxSelection2, 'Enable', 'off');
    set(handles.pushbuttonDeleteS2, 'Enable', 'off');
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% getRec3dModifyProjectGUIValues
% -------------------------------------------------------------------------
function [handles] = getRec3dModifyProjectGUIValues(handles)

% Simplify
TiltingGeometry = handles.inputREC3DPROJECT.PROJECT.TiltingGeometry;
ProjectionsPath1 = handles.inputREC3DPROJECT.PROJECT.ProjectionsPath1;
ProjectionsExt1 = handles.inputREC3DPROJECT.PROJECT.ProjectionsExt1;
MarkerfileOriginal1 = handles.inputREC3DPROJECT.PROJECT.MarkerfileOriginal1;
ProjectionsPath2 = handles.inputREC3DPROJECT.PROJECT.ProjectionsPath2;
ProjectionsExt2 = handles.inputREC3DPROJECT.PROJECT.ProjectionsExt2;
MarkerfileOriginal2 = handles.inputREC3DPROJECT.PROJECT.MarkerfileOriginal2;
NumOfMarkers = handles.inputREC3DPROJECT.PROJECT.NumOfMarkers;

% ELICIT NEW MARKERFILE 1

% ProjectionsName1
ProjectionsName1 = get(handles.listboxSelection1, 'String');
ProjectionsName1 = ProjectionsName1';

% NumOfProj1
NumOfProj1 = size(ProjectionsName1, 2);

% Find underline index and seperate index
for k=1:NumOfProj1
    
    UnderlineIndex1 = findstr(ProjectionsName1{k}, '_');
    WhichUnderline = size(UnderlineIndex1, 2);
    UnderlineIndex1 = UnderlineIndex1(1,WhichUnderline);
    
    for j = (UnderlineIndex1+1):(UnderlineIndex1+3)
         ProjectionsIndexCell1{k}(j) = ProjectionsName1{k}(j);
    end
end

for k=1:NumOfProj1
    ProjectionsIndex1(k) = str2num(ProjectionsIndexCell1{k});
end

% Copy indexed columns and header
NumOfLines1 = size(MarkerfileOriginal1.Value, 1);

for i=1:NumOfLines1
    for j=1:NumOfProj1
         for k=1:NumOfMarkers
            
Markerfile1.Value(i, j, k) = ...
         MarkerfileOriginal1.Value(i, ProjectionsIndex1(j), k);
         
         end
    end
end
Markerfile1.Header = MarkerfileOriginal1.Header;

% Tiltangles1
Tiltangles1(1:NumOfProj1) = Markerfile1.Value(1, 1:NumOfProj1, 1);

% RefProj1
[notneeded, RefProj1] = min(abs(Tiltangles1));

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

% ELICIT NEW MARKERFILE 2

% ProjectionsName2
if strcmp(TiltingGeometry, 'dualaxis')
    ProjectionsName2 = get(handles.listboxSelection2, 'String');
    ProjectionsName2 = ProjectionsName2';
end

% NumOfProj2
if strcmp(TiltingGeometry, 'dualaxis')
    NumOfProj2 = size(ProjectionsName2, 2);
else
    NumOfProj2 = 0;
end

% Find underline index and seperate index
if strcmp(TiltingGeometry, 'dualaxis')

for k=1:NumOfProj2
    
    UnderlineIndex2 = findstr(ProjectionsName2{k}, '_');
    WhichUnderline = size(UnderlineIndex2, 2);
    UnderlineIndex2 = UnderlineIndex2(1,WhichUnderline);
    
    for j = (UnderlineIndex2+1):(UnderlineIndex2+3)
         ProjectionsIndexCell2{k}(j) = ProjectionsName2{k}(j);
    end
end

for k=1:NumOfProj2
    ProjectionsIndex2(k) = str2num(ProjectionsIndexCell2{k});
end

% Copy indexed columns and header
NumOfLines2 = size(MarkerfileOriginal2.Value, 1);

for i=1:NumOfLines2
    for j=1:NumOfProj2
         for k=1:NumOfMarkers
            
Markerfile2.Value(i, j, k) = ...
         MarkerfileOriginal2.Value(i, ProjectionsIndex2(j), k);
         
         end
    end
end
Markerfile2.Header = MarkerfileOriginal2.Header;

% Tiltangles2
Tiltangles2(1:NumOfProj2) = Markerfile2.Value(1, 1:NumOfProj2, 1);

% RefProj2
[notneeded, RefProj2] = min(abs(Tiltangles2));

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

end

% Combine ProjectionsName1 and ProjectionsName2
for k=1:NumOfProj1
    ProjectionsName{k} = ProjectionsName1{k};
    ProjectionsPath{k} = [ProjectionsPath1];
end
if strcmp(TiltingGeometry, 'dualaxis')
    for k=1:NumOfProj2
        ProjectionsName{k+NumOfProj1} = ProjectionsName2{k};
        ProjectionsPath{k+NumOfProj1} = [ProjectionsPath2];
    end
end

% NumOfProj
NumOfProj = NumOfProj1 + NumOfProj2;

% Aktualize outputREC3DPROJECT
% PROJECT
handles.outputREC3DPROJECT.PROJECT.NumOfProj1 = NumOfProj1;
handles.outputREC3DPROJECT.PROJECT.RefProj1 = RefProj1;
handles.outputREC3DPROJECT.PROJECT.Tiltangles1 = Tiltangles1;
handles.outputREC3DPROJECT.PROJECT.MarkersOnProj1 = MarkersOnProj1;
handles.outputREC3DPROJECT.PROJECT.Markerfile1 = Markerfile1;
if strcmp(TiltingGeometry, 'dualaxis')
handles.outputREC3DPROJECT.PROJECT.NumOfProj2 = NumOfProj2;
handles.outputREC3DPROJECT.PROJECT.RefProj2 = RefProj2;
handles.outputREC3DPROJECT.PROJECT.Tiltangles2 = Tiltangles2;
handles.outputREC3DPROJECT.PROJECT.MarkersOnProj2 = MarkersOnProj2;
handles.outputREC3DPROJECT.PROJECT.Markerfile2 = Markerfile2;
end
handles.outputREC3DPROJECT.PROJECT.NumOfProj = NumOfProj;
% RECONSTRUCTION / Parameter
handles.outputREC3DPROJECT.PARAMETER.ProjectionsName = ProjectionsName;
handles.outputREC3DPROJECT.PARAMETER.ProjectionsPath = ProjectionsPath;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% pushbuttonT1toS1
% -------------------------------------------------------------------------
function pushbuttonT1toS1_Callback(hObject, eventdata, handles)

% Get selection tiltseries1
T1str = get(handles.listboxTiltseries1, 'String');
T1val = get(handles.listboxTiltseries1, 'Value');

% Get selection1
S1str = get(handles.listboxSelection1, 'String');

% Get marked cells
zaehler = 1;
for k=1:size(T1val, 2)
    index_marked_cell = T1val(k);
    marked_cells{zaehler} = T1str{index_marked_cell};
    zaehler = zaehler + 1;
end
marked_cells = marked_cells';

% Put marked cells to selection listbox
sizeS1str = size(S1str, 1);
for k=1:sizeS1str
    newS1str{k} = S1str{k};
end

for k=1:size(marked_cells, 1)
    newS1str{sizeS1str+k} = marked_cells{k};
end
newS1str = newS1str';

% Find double selection
rds_index = [];
zaehler = 1;
for k=1:size(newS1str, 1)
    
    for j=1:size(newS1str, 1)
    
        if (k<j) && (k~=j) && strcmp(newS1str{k}, newS1str{j}) == 1
       
        rds_index(zaehler) = j;    
        
        zaehler = zaehler + 1;
        
        end
    end
end
     
% Remove double selection
if isempty(rds_index) == 1
    newnewS1str = newS1str;
else

zaehler = 1;    
for k=1:size(newS1str, 1)
    
    DoThis = 'Hold';
    
    for j=1:size(rds_index, 2)
        
        if rds_index(j) == k
            DoThis = 'Delete';
        end
        
    end
      
    if strcmp(DoThis, 'Hold')
        
        newnewS1str{zaehler} = newS1str{k};
        
        zaehler = zaehler + 1;
    
    end

end
newnewS1str = newnewS1str';

end

% Sort entries
newsort=sort(newnewS1str);

% Set selection listbox
set(handles.listboxSelection1, 'String', newsort);

% Background color listbox
set(handles.listboxTiltseries1, 'BackgroundColor', 'White');
set(handles.listboxSelection1, 'BackgroundColor', 'White');

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------
    
% -------------------------------------------------------------------------
% pushbuttonDeleteS1
% -------------------------------------------------------------------------
function pushbuttonDeleteS1_Callback(hObject, eventdata, handles)

list_entries = get(handles.listboxSelection1, 'String');
index_selected = get(handles.listboxSelection1, 'Value');

if size(list_entries,1) == size(index_selected,2)
    return;
end

zaehler = 1;    
for k=1:size(list_entries, 1)
    
    WhatToDo = 'Hold';
    
    for j=1:size(index_selected, 2)
        
        if index_selected(j) == k
            
            WhatToDo = 'Del';
        
        end
        
    end
      
    if strcmp(WhatToDo, 'Hold')
        
        new_list_entries{zaehler} = list_entries{k};
        
        zaehler = zaehler + 1;
    
    end

end

new_list_entries = new_list_entries';

set(handles.listboxSelection1, 'Value', 1);
set(handles.listboxSelection1, 'String', new_list_entries);

% Background color listbox
set(handles.listboxTiltseries1, 'BackgroundColor', 'White');
set(handles.listboxSelection1, 'BackgroundColor', 'White');

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% pushbuttonT2toS2
% -------------------------------------------------------------------------
function pushbuttonT2toS2_Callback(hObject, eventdata, handles)

% Get selection tiltseries1
T2str = get(handles.listboxTiltseries2, 'String');
T2val = get(handles.listboxTiltseries2, 'Value');

% Get selection1
S2str = get(handles.listboxSelection2, 'String');

% Get marked cells
zaehler = 1;
for k=1:size(T2val, 2)
    index_marked_cell = T2val(k);
    marked_cells{zaehler} = T2str{index_marked_cell};
    zaehler = zaehler + 1;
end
marked_cells = marked_cells';

% Put marked cells to selection listbox
sizeS2str = size(S2str, 1);
for k=1:sizeS2str
    newS2str{k} = S2str{k};
end

for k=1:size(marked_cells, 1)
    newS2str{sizeS2str+k} = marked_cells{k};
end
newS2str = newS2str';

% Find double selection
rds_index = [];
zaehler = 1;
for k=1:size(newS2str, 1)
    
    for j=1:size(newS2str, 1)
    
        if (k<j) && (k~=j) && strcmp(newS2str{k}, newS2str{j}) == 1
       
        rds_index(zaehler) = j;    
        
        zaehler = zaehler + 1;
        
        end
    end
end
     
% Remove double selection
if isempty(rds_index) == 1
    newnewS2str = newS2str;
else

zaehler = 1;    
for k=1:size(newS2str, 1)
    
    DoThis = 'Hold';
    
    for j=1:size(rds_index, 2)
        
        if rds_index(j) == k
            DoThis = 'Delete';
        end
        
    end
      
    if strcmp(DoThis, 'Hold')
        
        newnewS2str{zaehler} = newS2str{k};
        
        zaehler = zaehler + 1;
    
    end

end
newnewS2str = newnewS2str';

end

% Sort entries
newsort=sort(newnewS2str);

% Set selection listbox
set(handles.listboxSelection2, 'String', newsort);

% Background color listbox
set(handles.listboxTiltseries2, 'BackgroundColor', 'White');
set(handles.listboxSelection2, 'BackgroundColor', 'White');

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% pushbuttonDelete2
% -------------------------------------------------------------------------
function pushbuttonDeleteS2_Callback(hObject, eventdata, handles)

list_entries = get(handles.listboxSelection2, 'String');
index_selected = get(handles.listboxSelection2, 'Value');

if size(list_entries,1) == size(index_selected,2)
    return;
end

zaehler = 1;    
for k=1:size(list_entries, 1)
    
    WhatToDo = 'Hold';
    
    for j=1:size(index_selected, 2)
        
        if index_selected(j) == k
            
            WhatToDo = 'Del';
        
        end
        
    end
      
    if strcmp(WhatToDo, 'Hold')
        
        new_list_entries{zaehler} = list_entries{k};
        
        zaehler = zaehler + 1;
    
    end

end

new_list_entries = new_list_entries';

set(handles.listboxSelection2, 'Value', 1);
set(handles.listboxSelection2, 'String', new_list_entries);

% Background color listbox
set(handles.listboxTiltseries2, 'BackgroundColor', 'White');
set(handles.listboxSelection2, 'BackgroundColor', 'White');

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------

