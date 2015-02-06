function varargout = tom_markAlignRigidBody(varargin)
% TOM_MARKALIGNRIGIDBODY M-file for tom_markAlignRigidBody.fig
%      TOM_MARKALIGNRIGIDBODY, by itself, creates a new TOM_MARKALIGNRIGIDBODY or raises the existing
%      singleton*.
%
%      H = TOM_MARKALIGNRIGIDBODY returns the handle to a new TOM_MARKALIGNRIGIDBODY or the handle to
%      the existing singleton*.
%
%      TOM_MARKALIGNRIGIDBODY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_MARKALIGNRIGIDBODY.M with the given input arguments.
%
%      TOM_MARKALIGNRIGIDBODY('Property','Value',...) creates a new TOM_MARKALIGNRIGIDBODY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_markAlignRigidBody_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_markAlignRigidBody_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OpeningFcn, ...
                   'gui_OutputFcn',  @Outputfcn, ...
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
%% --- Executes just before tom_markAlignRigidBody is made visible.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_markAlignRigidBody (see VARARGIN)

markGui_hfigure = findall(0, 'Type', 'figure', 'Tag', 'figure_markGui');
if (length(markGui_hfigure) ~= 1)
    %errordlg('Call tom_markGui_Alignment3dItr from tom_markGui', 'markGui_Alignment3dItr', 'modal');
    warning('Call tom_markGui_Alignment3dItr from tom_markGui');
    tom_markGui;
    delete(handles.figure_markGui_Alignment3dItr);
    return;
end;
markGui_data = guidata(markGui_hfigure);


handles.data2 = markGui_data.data;
handles.fcn_tom_markGui = markGui_data.fcn_local;
handles.data.markGui_hfigure = markGui_hfigure;

handles.data.guiloadad = false;

if (isempty(handles.data2.markersets) || ...
    handles.data2.sel_markerset < 1 || ...
    handles.data2.sel_markerset > length(handles.data2.markersets) || ...
    isempty(handles.data2.filenames))
    errordlg('Load the markerset in tom_markGui first', 'markGui_Alignment3dItr', 'modal');
    delete(handles.figure_markGui_Alignment3dItr);
end;


markerset = handles.data2.markersets(handles.data2.sel_markerset);


% Construct the alignment configuration....
if (isfield(markerset, 'alignment3ditr') && isstruct(markerset.alignment3ditr) && numel(markerset.alignment3ditr)==1)
    aligncfg = markerset.alignment3ditr;
else
    aligncfg = struct();
end;

ni = size(markerset.markerset, 2);
nm = size(markerset.markerset, 3);

[emheader_tiltangles, emheader_imsize] = handles.fcn_tom_markGui.getTiltangles(handles.data2.emheader);

if (~isfield(aligncfg, 'tiltangles') || numel(aligncfg.tiltangles)~=ni)
    aligncfg.tiltangles = emheader_tiltangles;
    aligncfg.done = false;
end;
if (~isfield(aligncfg, 'imsize') || numel(aligncfg.imsize)~=1 || ~isnumeric(aligncfg.imsize) ||~(aligncfg.imsize>0))
    aligncfg.imsize = mean(emheader_imsize(:));
    aligncfg.done = false;
end;    
if (~isfield(aligncfg, 'psi') || numel(aligncfg.psi)~=1 || ~isnumeric(aligncfg.psi))
    aligncfg.psi = 0;
    aligncfg.done = false;
end;
if (~isfield(aligncfg, 'sigma') || numel(aligncfg.sigma)~=1 || ~isnumeric(aligncfg.sigma))
    aligncfg.sigma = 0;
    aligncfg.done = false;
end;
if (~isfield(aligncfg, 'tx') || ~isnumeric(aligncfg.tx) || numel(aligncfg.tx)~=ni)
    aligncfg.tx = zeros(1,ni);
    aligncfg.done = false;
end;
if (~isfield(aligncfg, 'ty') || ~isnumeric(aligncfg.ty) || numel(aligncfg.ty)~=ni)
    aligncfg.ty = zeros(1,ni);
    aligncfg.done = false;
end;
if (~isfield(aligncfg, 'imintilt') || numel(aligncfg.imintilt)~=1 || round(aligncfg.imintilt)~=aligncfg.imintilt || aligncfg.imintilt<1 || aligncfg.imintilt>ni)
    aligncfg.imintilt = returnLastImTilt;
    if (isempty(aligncfg.imintilt))
        aligncfg.imintilt = returnLastImTilt(find(min(abs(aligncfg.tiltangles)) == abs(aligncfg.tiltangles), 1, 'first'));
    end;
    aligncfg.done = false;
end;

if (~isfield(aligncfg, 'done'))
    aligncfg.done = false;
elseif (~islogical(aligncfg.done))
    aligncfg.done = false;
end;


handles.data.aligncfg = aligncfg;
handles.data.aligndata = [];
handles.data.markerset = markerset;

handles = loadGui(handles);

handles.status = 'cancel';
guidata(handles.figure_markGui_Alignment3dItr, handles);


%set(handles.figure_markGui_Alignment3dItr, 'WindowStyle', 'modal');

if (0)    
    warning('No uiwait for debugging purpose!');
    Callback_pushbutton_align(hObject, eventdata, handles)
else
    uiwait(handles.figure_markGui_Alignment3dItr);
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- Outputs from this function are returned to the command line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = Outputfcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = 'cancel';
varargout{2} = [];
varargout{3} = [];

if (~ishandle(hObject))
    return;
end;


switch (handles.status)
    case {'ok', 'okimport'}
        varargout{1} = handles.status;
        varargout{2} = handles.data.aligncfg;
        varargout{3} = handles.data.aligndata;        
    case 'cancel'
    otherwise
        warning('Reached unexpected point in program run');
end;
delete(handles.figure_markGui_Alignment3dItr); 









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function i = returnLastImTilt(setval)

persistent imintilt;
if (exist('setval', 'var'))
    imintilt = setval;
end;
i = imintilt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = loadGui(handles)

markerset = handles.data.markerset;
aligncfg = handles.data.aligncfg;


msize = [1, 1, 1];
[msize(1), msize(2), msize(3)] = size(markerset.markerset);

set(handles.figure_markGui_Alignment3dItr, 'Name', ['Rigid Body Alignment with iteration for "' markerset.name '"'])
set(handles.text_nprojections, 'String', num2str(msize(2)));
setNumericEditField(handles.edit_imsize, handles.data.aligncfg.imsize);


ss = cell(1,msize(2));
for (i=1:msize(2))
    ss{i} = [num2str(i) ': ' handles.data2.fileinfo(i).filename ' (' sprintf('%7.3f', handles.data.aligncfg.tiltangles(i)) '??)'];
end;
set(handles.popupmenu_refim, 'String', ss, 'Value', handles.data.aligncfg.imintilt);


otheralignments = false;
for (i=1:length(handles.data2.markersets))
    if (i ~= handles.data2.sel_markerset && ...
        isfield(handles.data2.markersets(i), 'alignment3ditr') && ...
        isstruct(handles.data2.markersets(i).alignment3ditr) && ...
        numel(handles.data2.markersets(i).alignment3ditr)==1 && ...
        isfield(handles.data2.markersets(i).alignment3ditr, 'done') && ...
        handles.data2.markersets(i).alignment3ditr.done)
        otheralignments = true;
        break;
    end;
    if (isfield(handles.data2.markersets(i), 'alignment3d') && ...
        isstruct(handles.data2.markersets(i).alignment3d) && ...
        numel(handles.data2.markersets(i).alignment3d)==1 && ...
        isfield(handles.data2.markersets(i).alignment3d, 'done') && ...
        handles.data2.markersets(i).alignment3d.done)
        otheralignments = true;
        break;
    end;
end;

if (otheralignments)
    set(handles.radiobutton_init_other, 'Enable', 'on', 'Value', 0);
else
    set(handles.radiobutton_init_other, 'Enable', 'off', 'Value', 0);
end;

if (aligncfg.done)
    set(handles.togglebutton_align_view, 'Visible', 'on');
    set(handles.radiobutton_init_last, 'Enable', 'on', 'Value', 1);
    %handles = Callback_togglebuttons(handles.togglebutton_align_view, [], handles);
    handles = Callback_togglebuttons(handles.togglebutton_align_compute, [], handles);
else
    set(handles.togglebutton_align_view, 'Visible', 'off');
    set(handles.radiobutton_init_last, 'Enable', 'off', 'Value', 0);
    set(handles.radiobutton_init_other, 'Value', otheralignments);
    handles = Callback_togglebuttons(handles.togglebutton_align_compute, [], handles);
end;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setRefMarkXCoords(handles)

pos = handles.data.markerset(1:2, handles.data.aligncfg.imintilt, handles.data.aligncfg.refmark)';

if (~all(pos > 0) || ~all(pos <= handles.data.aligncfg.imsize))
    set(handles.radiobutton_refmarkx_manually, 'Value', true);
    set(handles.radiobutton_refmarkx_auto, 'Value', false, 'Enable', 'off');
    set(handles.edit_refmarkx_x_auto, 'String', 'NONE');
    set(handles.edit_refmarkx_y_auto, 'String', 'NONE');
else
    set(handles.radiobutton_refmarkx_auto, 'Enable', 'on');
    set(handles.edit_refmarkx_x_auto, 'String', num2str(pos(1)));
    set(handles.edit_refmarkx_y_auto, 'String', num2str(pos(2)));
end;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setNumericEditField(hObject, value)

userdata = get(hObject, 'Userdata');

userdata.numvalue = value;
set(hObject, 'String', num2str(value), 'Userdata', userdata);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeNumericEditField(hObject);


userdata = get(hObject, 'Userdata');

val = str2double(get(hObject, 'String'));
if (isnan(val))
    warndlg(['The entered number is not a number as expected.'], 'Warning', 'modal');
else
    setNumericEditField(hObject, val);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_cancel(hObject, eventdata, handles)


handles.status = 'cancel';
guidata(handles.figure_markGui_Alignment3dItr, handles); 

uiresume();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_align(hObject, eventdata, handles)




if (get(handles.togglebutton_align_compute, 'Value'))

    markerset = handles.data.markerset.markerset;

    aligncfg = handles.data.aligncfg;

    userdata = get(handles.edit_imsize, 'UserData');
    aligncfg.imsize = userdata.numvalue;

    aligncfg.imintilt = get(handles.popupmenu_refim, 'Value');

    initshift = [];
    if (get(handles.radiobutton_init_last, 'Value'))
        [P, x, X_trian] = tom_mark_cvaf_alignment3d_reproj(aligncfg.tiltangles, aligncfg.psi, aligncfg.tx, aligncfg.ty, aligncfg.imsize, [], markerset, true);
        X_trian = mean(X_trian(:,all(isfinite(X_trian), 1)), 2);
        if (isempty(X_trian))
            initshift = [aligncfg.tx; aligncfg.ty];
        else
            X_trian(4) = 1;
            initshift = nan(3,size(P,3));
            for (i=1:size(P,3))
                initshift(1:3,i) = P(1:3,1:4,i) * X_trian;
            end;
            initshift = initshift(1:2,:);
        end;
    elseif (get(handles.radiobutton_init_other, 'Value'))
        [answw, initshift] = getOtherInitialisation(handles);
    end;
    
    
    
    delete(findall(0, 'type', 'figure', 'Tag', 'tom_markGui_Alignment3dItr'));
    hwaitbar = waitbar(0, 'please wait... doing alignment', 'CreateCancelBtn', 'set(gcf, ''Userdata'', true);', 'Userdata', false, 'Name', 'alignment3d iterative', 'Tag', 'tom_markGui_Alignment3dItr');
    set(findobj(hwaitbar, 'Tag', 'TMWWaitbarCancelButton', 'Type', 'uicontrol', 'Style', 'pushbutton'), 'String', 'Break');
    %set(waitbarh, 'WindowStyle', 'modal');
    set(hwaitbar, 'WindowStyle', 'normal');

    [P, X, aligncfg] = tom_mark_cvaf_alignment3dItr(markerset, aligncfg, 50000, 100000, initshift, hwaitbar);

    close(hwaitbar);
    delete(hwaitbar);


    x = ones(3, size(markerset,2), size(markerset,3));
    for (i=1:size(markerset,2))
        x(1:3,i,:) = P(1:3,1:4,i) * X;
    end;
    handles.data.aligndata.P = P;
    handles.data.aligndata.x = x(1:2,:,:);
    handles.data.aligndata.X = X(1:3,:,:);
    handles.data.aligncfg = aligncfg;
    handles.data.aligncfg.done = true;
    
    handles.data.guiloadad = false;
    
    
    set(handles.radiobutton_init_last, 'Value', 1, 'Enable', 'on');
    set(handles.radiobutton_init_other, 'Value', 0);
    set(handles.togglebutton_align_view, 'Visible', 'on');    
    set(handles.pushbutton_import_reproj, 'Visible', 'on');
    
    guidata(handles.figure_markGui_Alignment3dItr, handles);
    
    handles = Callback_togglebuttons(handles.togglebutton_align_view, [], handles);

elseif (get(handles.togglebutton_align_view, 'Value'))
    handles.status = 'ok';
    guidata(handles.figure_markGui_Alignment3dItr, handles);
    uiresume();
else
    warning('Reached unexpected point in program run');
end;
        
        
         
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [answw, initshift] = getOtherInitialisation(handles)

answw = false;
initshift = [];


nother = cell(1, length(handles.data2.markersets));
ss = {cell(1, length(handles.data2.markersets))};
useidx = [];
for (i=1:length(handles.data2.markersets))
    notheri = {};
    if (isfield(handles.data2.markersets(i), 'alignment3d') && ...
        isstruct(handles.data2.markersets(i).alignment3d) && ...
        numel(handles.data2.markersets(i).alignment3d)==1 && ...
        isfield(handles.data2.markersets(i).alignment3d, 'done') && ...
        handles.data2.markersets(i).alignment3d.done)
        otheralignments = true;
        notheri = { notheri{:}, 'alignment3d'};
    end;
    if (i ~= handles.data2.sel_markerset && ...
        isfield(handles.data2.markersets(i), 'alignment3ditr') && ...
        isstruct(handles.data2.markersets(i).alignment3ditr) && ...
        numel(handles.data2.markersets(i).alignment3ditr)==1 && ...
        isfield(handles.data2.markersets(i).alignment3ditr, 'done') && ...
        handles.data2.markersets(i).alignment3ditr.done)
        notheri = { notheri{:}, 'alignment3ditr'};
    end;
    nother{i} = notheri;
    if (~isempty(notheri))
        useidx(end+1) = i;
    end;
    ss{i} = [num2str(i) ': ' handles.data2.markersets(i).name ' (' num2str(size(handles.data2.markersets(i).markerset,3)) ' marker)'];
end;

notherverbose = cell(size(nother));
for (i=1:length(nother))
    notheri = nother{i};
    notherverbosei = {};
    for (j=1:length(notheri))
        switch (notheri{j})
            case 'alignment3d'
                notherverbosei{j} = 'Rigid Body Alignmentd3d';
            case 'alignment3ditr'
                notherverbosei{j} = 'Rigid Body Alignmentd3d (iterative)';
            otherwise
                notherverbosei{j} = notheri{j};
        end;        
    end; 
    notherverbose{i} = notherverbosei;
end;

if (isempty(useidx))
    warning('This should not happen!!!! check reason');
    return;
end

[idx_markerset, answ] = listdlg('ListString', ss(useidx), 'SelectionMode', 'single', 'Listsize', [350,200], 'InitialValue', 1, 'Name', 'Iterative alignment', 'Promptstring', 'Take the initalisation from which markerset?');

if (~answ || numel(idx_markerset)~=1)
    return;
end;
idx_markerset = useidx(idx_markerset);

notherverbosei = notherverbose{idx_markerset};
[idx_type, answ] = listdlg('ListString', notherverbosei, 'SelectionMode', 'single', 'Listsize', [350,200], 'InitialValue', 1, 'Name', 'Iterative alignment', 'Promptstring', 'Which alignment-type?');

if (~answ || numel(idx_type)~=1)
    return;
end;

notheri = nother{idx_markerset};
switch (notheri{idx_type})
    case {'alignment3d', 'alignment3ditr'}
        aligncfg2 = handles.data2.markersets(idx_markerset).(notheri{idx_type});
        [P, xdummy, X_trian] = tom_mark_cvaf_alignment3d_reproj(aligncfg2.tiltangles, aligncfg2.psi, aligncfg2.tx, aligncfg2.ty, aligncfg2.imsize, [], handles.data.markerset.markerset, true);
        X_trian = mean(X_trian(:,all(isfinite(X_trian), 1)), 2);
        if (isempty(X_trian))
            initshift = [aligncfg2.tx; aligncfg2.ty];
        else
            X_trian(4) = 1;
            initshift = nan(3,size(P,3));
            for (i=1:size(P,3))
                initshift(1:3,i) = P(1:3,1:4,i) * X_trian;
            end;
            initshift = initshift(1:2,:);
        end;
    otherwise
        warning('This should not happen!!!! check reason');
        return;
end;

answw = true;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = Callback_togglebuttons(hObject, eventdata, handles)
 

tag = get(hObject, 'Tag');
switch (tag)
    case 'togglebutton_align_compute'
        set(handles.togglebutton_align_compute, 'Value', 1);
        set(handles.togglebutton_align_view, 'Value', 0);
        set(handles.uipanel_align_compute, 'Visible', 'on');
        set(handles.uipanel_align_view, 'Visible', 'off');
        set(handles.pushbutton_align, 'String', 'do alignment');
        set(handles.pushbutton_import_reproj, 'Visible', 'off');
    case 'togglebutton_align_view'
        set(handles.togglebutton_align_compute, 'Value', 0);
        set(handles.togglebutton_align_view, 'Value', 1);
        set(handles.uipanel_align_compute, 'Visible', 'off');
        set(handles.uipanel_align_view, 'Visible', 'on');
        set(handles.pushbutton_align, 'String', 'Ok');
        set(handles.pushbutton_import_reproj, 'Visible', 'on');
        handles = loadGuiView(handles);
    otherwise
        warning('Reached unexpected point in program run');
end;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
function handles = loadGuiView(handles)        

if (~handles.data.guiloadad)
    
    aligncfg = handles.data.aligncfg;
    
    if (isempty(handles.data.aligndata) || ~isstruct(handles.data.aligndata) || ~isfield(handles.data.aligndata, 'x') || any(size(handles.data.markerset.markerset)~=size(handles.data.aligndata.x)))
        [handles.data.aligndata.P, dummyx, handles.data.aligndata.X, handles.data.aligndata.x] = ...
            tom_mark_cvaf_alignment3d_reproj(aligncfg.tiltangles, aligncfg.psi, aligncfg.tx, aligncfg.ty, aligncfg.imsize, [], handles.data.markerset.markerset, true);
        handles.data.aligndata.x = handles.data.aligndata.x(1:2,:,:);
    end;

    distance = sqrt(sum((handles.data.markerset.markerset - handles.data.aligndata.x) .^ 2, 1));
    
    set(handles.text_view_tiltaxis, 'String', [num2str(aligncfg.psi * 180/pi()) '??']);
    set(handles.text_view_distance, 'String', [sprintf('%7.4f', min(distance(:))) ' ~ ' sprintf('%7.4f', mean(distance(isfinite(distance)))) ' ~ ' sprintf('%7.4f', max(distance(:)))]);
    
    guidata(handles.figure_markGui_Alignment3dItr, handles);

end;
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_import_reproj(hObject, eventdata, handles)
        
        
handles.status = 'okimport';
guidata(handles.figure_markGui_Alignment3dItr, handles);
uiresume();
    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_view_plot3(hObject, eventdata, handles)
 


ftag = ['tom_markAlignRigidBody_plot3'];
h = findobj('Type', 'figure', 'Tag', ftag);

if (length(h) ~= 1)
    delete(h);
    h = figure();
    set(h, 'Tag', ftag);
end;

set(0, 'CurrentFigure', h);


set(handles.figure_markGui_Alignment3dItr, 'WindowStyle', 'normal');
X = handles.data.aligndata.X;    

cla;
plot3(X(1,:), X(2,:), X(3,:), 'b.');
axis equal;
        
        
        
        

