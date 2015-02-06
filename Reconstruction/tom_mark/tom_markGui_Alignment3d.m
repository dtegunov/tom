function varargout = tom_markGui_Alignment3d(varargin)
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
                   'gui_OpeningFcn', @tom_markAlignRigidBody_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_markAlignRigidBody_OutputFcn, ...
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
function tom_markAlignRigidBody_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_markAlignRigidBody (see VARARGIN)

if (length(varargin) ~= 4)
    error([mfilename() ' were called with the wrong number of arguments']);
end;

handles.data.markerset = varargin{1};
handles.data.msname = varargin{2};
handles.data.imname = varargin{3};
handles.data.aligncfg = varargin{4};



handles.data.aligncfg_old = handles.data.aligncfg;
handles.data.aligncfg = parseAlignCfg(handles);

handles.data.markerset_inbound_mask = all(isfinite(handles.data.markerset), 1);

handles.data.aligndata = [];

handles = loadGui(handles);

handles.status = 'cancel';
guidata(handles.figure_markGui_Alignment3d, handles);


%set(handles.figure_markGui_Alignment3d, 'WindowStyle', 'modal');

uiwait(handles.figure_markGui_Alignment3d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function aligncfg = parseAlignCfg(handles)

aligncfg = handles.data.aligncfg;
ni = size(handles.data.markerset,2);
nm = size(handles.data.markerset,3);

if (isnumeric(aligncfg) && numel(aligncfg)==ni)
    aligncfg = struct('tiltangles', {aligncfg});
elseif (isstruct(aligncfg) && numel(aligncfg)==1)
else
    error([mfilename ' needs the tiltangles of the images']);
end;


if (~isfield(aligncfg, 'tiltangles') || ~isnumeric(aligncfg.tiltangles) || numel(aligncfg.tiltangles)~=ni)
    error([mfilename ' needs the tiltangles of the images']);
end;
if (~isfield(aligncfg, 'imsize') || numel(aligncfg.imsize)~=1 || ~isnumeric(aligncfg.imsize) ||~(aligncfg.imsize>0))
    aligncfg.imsize = 1024;
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
inbound = all(handles.data.markerset>0 & handles.data.markerset<aligncfg.imsize, 1);    
if (~isfield(aligncfg, 'refmark') || numel(aligncfg.refmark)~=1 || round(aligncfg.refmark)~=aligncfg.refmark || aligncfg.refmark<1 || aligncfg.refmark>nm)
    [i, aligncfg.refmark] = max(squeeze(sum(inbound, 2) .* inbound(1,aligncfg.imintilt,:)));
    aligncfg.done = false;
end;    
if (~isfield(aligncfg, 'refmarkX') || numel(aligncfg.refmarkX)~=3 || ~all(isfinite(aligncfg.refmarkX)))
    if (inbound(1,aligncfg.imintilt, aligncfg.refmark))
        aligncfg.refmarkX = [handles.data.markerset(1:2, aligncfg.imintilt, aligncfg.refmark)', aligncfg.imsize/2];
    else
        aligncfg.refmarkX = ones(1,3) * aligncfg.imsize/2;
    end;
    aligncfg.done = false;
end;

if (~isfield(aligncfg, 'done'))
    aligncfg.done = false;
elseif (~islogical(aligncfg.done))
    aligncfg.done = false;
end;





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
[msize(1), msize(2), msize(3)] = size(markerset);

set(handles.figure_markGui_Alignment3d, 'Name', ['Rigid Body Alignment for "' handles.data.msname '"'])
set(handles.text_nprojections, 'String', num2str(msize(2)));
setNumericEditField(handles.edit_imsize, handles.data.aligncfg.imsize);
SetRefMarkList(handles);

ss = cell(1,msize(2));
for (i=1:msize(2))
    ss{i} = [num2str(i) ': ' handles.data.imname{i} ' (' sprintf('%7.3f', handles.data.aligncfg.tiltangles(i)) ' deg)'];
end;
set(handles.popupmenu_refim, 'String', ss, 'Value', handles.data.aligncfg.imintilt);

set(handles.radiobutton_refmarkx_manually, 'Value', true);
set(handles.radiobutton_refmarkx_auto, 'Value', false);

setNumericEditField(handles.edit_refmarkx_x, handles.data.aligncfg.refmarkX(1));
setNumericEditField(handles.edit_refmarkx_y, handles.data.aligncfg.refmarkX(2));
setNumericEditField(handles.edit_refmarkx_z, handles.data.aligncfg.refmarkX(3));
 
Callback_popupmenu_refmark(handles.popupmenu_refmark, [], handles);


if (aligncfg.done)
    set(handles.togglebutton_align_view, 'Visible', 'on');
    set(handles.pushbutton_import_reproj, 'Visible', 'on');
    handles = Callback_togglebuttons(handles.togglebutton_align_view, [], handles);
else
    set(handles.togglebutton_align_view, 'Visible', 'off');
    set(handles.pushbutton_import_reproj, 'Visible', 'off');
    handles = Callback_togglebuttons(handles.togglebutton_align_compute, [], handles);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- Outputs from this function are returned to the command line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_markAlignRigidBody_OutputFcn(hObject, eventdata, handles) 
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
delete(handles.figure_markGui_Alignment3d); 






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
%% Sets the RefMarker-Popuplist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetRefMarkList(handles)

ni = size(handles.data.markerset, 2);
nm = size(handles.data.markerset, 3);
imintilt = handles.data.aligncfg.imintilt;

ss = cell(1,nm); 
for (i=1:nm)
    if (handles.data.markerset_inbound_mask(1,imintilt,i))
        ss{i} = ['#' num2str(i) ' (in ' num2str(sum(handles.data.markerset_inbound_mask(1,:,i))) ' images)*'];
    else
        ss{i} = ['#' num2str(i) ' (in ' num2str(sum(handles.data.markerset_inbound_mask(1,:,i))) ' images)'];
    end;
end;


set(handles.popupmenu_refmark, 'String', ss, 'Value', handles.data.aligncfg.refmark);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = Callback_popupmenu_refmark(hObject, eventdata, handles)

refmark = get(handles.popupmenu_refmark, 'Value');

handles.data.aligncfg.refmark = refmark;

setRefMarkXCoords(handles);

guidata(handles.figure_markGui_Alignment3d, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = Callback_popupmenu_refim(hObject, eventdata, handles)

handles.data.aligncfg.imintilt = get(handles.popupmenu_refim, 'Value');

SetRefMarkList(handles);
setRefMarkXCoords(handles);

guidata(handles.figure_markGui_Alignment3d, handles);



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
guidata(handles.figure_markGui_Alignment3d, handles); 

uiresume();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_align(hObject, eventdata, handles)




if (get(handles.togglebutton_align_compute, 'Value'))

    markerset = handles.data.markerset;

    aligncfg = handles.data.aligncfg;

    userdata = get(handles.edit_imsize, 'UserData');
    aligncfg.imsize = userdata.numvalue;

    inbound_mask = all(markerset>0 & markerset<aligncfg.imsize, 1);
    
    
    if (get(handles.radiobutton_refmarkx_auto, 'Value') && inbound_mask(1,aligncfg.imintilt,aligncfg.refmark))
        aligncfg.refmarkX(1:2) = markerset(1:2, aligncfg.imintilt,aligncfg.refmark);
    else
        userdata = get(handles.edit_refmarkx_x, 'UserData');
        aligncfg.refmarkX(1) = userdata.numvalue;
        userdata = get(handles.edit_refmarkx_y, 'UserData');
        aligncfg.refmarkX(2) = userdata.numvalue;
    end;
    userdata = get(handles.edit_refmarkx_z, 'UserData');
    aligncfg.refmarkX(3) = userdata.numvalue;

    [aligncfg.psi, aligncfg.tx, aligncfg.ty, aligncfg.sigma] = tom_mark_cvaf_alignment3d(markerset, aligncfg.tiltangles, aligncfg.refmark, aligncfg.refmarkX, aligncfg.imsize);
    aligncfg.done = true;
    
    handles.data.aligndata = [];
    
    handles.data.aligncfg = aligncfg;
    
    
    set(handles.togglebutton_align_view, 'Visible', 'on');    
    set(handles.pushbutton_import_reproj, 'Visible', 'on');
    
    guidata(handles.figure_markGui_Alignment3d, handles);
    
    handles = Callback_togglebuttons(handles.togglebutton_align_view, [], handles);


elseif (get(handles.togglebutton_align_view, 'Value'))
    handles.status = 'ok';
    guidata(handles.figure_markGui_Alignment3d, handles);
    uiresume();
else
    warning('Reached unexpected point in program run');
end;
        
        
         
        


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
    case 'togglebutton_align_view'
        set(handles.togglebutton_align_compute, 'Value', 0);
        set(handles.togglebutton_align_view, 'Value', 1);
        set(handles.uipanel_align_compute, 'Visible', 'off');
        set(handles.uipanel_align_view, 'Visible', 'on');
        set(handles.pushbutton_align, 'String', 'Ok');
        handles = loadGuiView(handles);
    otherwise
        warning('Reached unexpected point in program run');
end;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
function handles = loadGuiView(handles)        

if (isempty(handles.data.aligndata))
    
    aligncfg = handles.data.aligncfg;
    markerset = handles.data.markerset;

    [P, xdummy, X, x] = tom_mark_cvaf_alignment3d_reproj(aligncfg.tiltangles, aligncfg.psi, aligncfg.tx, aligncfg.ty, aligncfg.imsize, [], markerset, false);
    
    x = x(1:2,:,:);
    X = X(1:3,:);

    handles.data.aligndata.P = P;
    handles.data.aligndata.x = x;
    handles.data.aligndata.X = X;
    
    
    distance = sqrt(sum((markerset - x) .^ 2, 1));
    
 
    
    set(handles.text_view_tiltaxis, 'String', [num2str(aligncfg.psi * 180/pi()) ' deg']);
    set(handles.text_view_rms, 'String', [num2str(aligncfg.sigma) ' pixel']);
    set(handles.text_view_distance, 'String', [sprintf('%7.4f', min(distance(:))) ' ~ ' sprintf('%7.4f', mean(distance(isfinite(distance)))) ' ~ ' sprintf('%7.4f', max(distance(:)))]);
    
    guidata(handles.figure_markGui_Alignment3d, handles);

end;
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_import_reproj(hObject, eventdata, handles)
        
        
handles.status = 'okimport';
guidata(handles.figure_markGui_Alignment3d, handles);
uiresume();
    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_view_plot3(hObject, eventdata, handles)
 


ftag = ['tom_markAlignRigidBody_plot3'];
h = findobj('Type', 'figure', 'Tag', ftag);

if (isempty(h))
    h = figure();
    set(h, 'Tag', ftag);
end;

figure(h);

set(handles.figure_markGui_Alignment3d, 'WindowStyle', 'normal');
X = handles.data.aligndata.X;    

cla;
axis equal;
plot3(X(1,:), X(2,:), X(3,:), 'b.');
        
        
        
        

