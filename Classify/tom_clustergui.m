function varargout = tom_clustergui(varargin)
%TOM_CLUSTERGUI creates ...
%
%   tom_clustergui(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout   		...
%
%EXAMPLE
%   ... = tom_clustergui(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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

%error(nargchk(0, 1, nargin, 'struct'))

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @tom_clustergui_OpeningFcn, ...
    'gui_OutputFcn',  @tom_clustergui_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Opening Function                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tom_clustergui_OpeningFcn(hObject, eventdata, handles, varargin)

if nargin == 4
    in_struct = varargin{1};
    eigstring = '';
    for i=1:length(in_struct.no_eigenvectors)
        eigstring = [eigstring ',' num2str(in_struct.no_eigenvectors(i))];
    end
    eigstring = eigstring(2:end);
    set(handles.edit_cluster_eigenvectors,'String',eigstring);
    set(handles.edit_cluster_numclasses,'String',num2str(in_struct.no_classes));
    handles.scores = in_struct.scores;
    handles.maxeigs = size(in_struct.scores,1);
else
    error('No input given!');
end

% Choose default command line output for tom_clustergui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_clustergui wait for user response (see UIRESUME)
uiwait(handles.tom_clustergui);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_clustergui_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;

close(handles.tom_clustergui);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Algorithm select box                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dropdown_cluster_algorithm_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns dropdown_cluster_algorithm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dropdown_cluster_algorithm

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calc optimal number of classes                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_cluster_calcoptcluster_Callback(hObject, eventdata, handles)

eigs = str2num(get(handles.edit_cluster_eigenvectors,'String'));
inp = handles.scores(eigs,:)';

[dists opt_num_of_classes]=tom_find_opt_num_of_classes(inp,str2num(get(handles.edit_cluster_optclasses,'String')),'k-means',[],1);
figure; bar(dists(:,2),dists(:,1));

set(handles.edit_cluster_numclasses,'String',num2str(opt_num_of_classes));

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Start Clustering                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_cluster_start_Callback(hObject, eventdata, handles)

%select eigenvectors

eigs = str2num(get(handles.edit_cluster_eigenvectors,'String'));
numclasses = str2num(get(handles.edit_cluster_numclasses,'String'));
%prepare input from scores
inp = handles.scores(eigs,:)';

[classes centroid sum_dist dists] = kmeans(inp,numclasses);

refpoint = 1;
newclasses = 1;
centroid2 = centroid;
centroid2(1,:) = zeros(length(eigs),1)-9e100;
%sort classes
for i=1:numclasses-1
    refpoint = tom_nearestpoint(centroid(refpoint,:), centroid2);
    centroid2(refpoint,:) = zeros(length(eigs),1)-9e100;
    newclasses = [newclasses refpoint];
end

for i=1:numclasses
    newcent(i,:) = centroid(newclasses(i),:);
end


classes2 = classes;
for i=1:numclasses
    idx = find(classes==newclasses(i));
    classes2(idx) = i;
end

handles.output = struct();
handles.output.centroid = newcent;
handles.output.sum_dist = sum_dist;
handles.output.dists = dists;
handles.output.classes = classes2';
handles.output.eigenvectors = eigs;
handles.output.noclasses = str2num(get(handles.edit_cluster_numclasses,'String'));
uiresume(handles.tom_clustergui);
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Cancel                                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_cluster_cancel_Callback(hObject, eventdata, handles)

handles.output = [];
guidata(hObject, handles);
uiresume(handles.tom_clustergui);
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_cluster_numclasses_Callback(hObject, eventdata, handles)
function edit_cluster_numclasses_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cluster_eigenvectors_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cluster_eigenvectors_Callback(hObject, eventdata, handles)
function dropdown_cluster_algorithm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cluster_optclasses_Callback(hObject, eventdata, handles)
function edit_cluster_optclasses_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




