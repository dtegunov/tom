function varargout = tom_stattool(varargin)
% TOM_STATTOOL M-file for tom_stattool.fig
%      TOM_STATTOOL, by itself, creates a new TOM_STATTOOL or raises the existing
%      singleton*.
%
%      H = TOM_STATTOOL returns the handle to a new TOM_STATTOOL or the handle to
%      the existing singleton*.
%
%      TOM_STATTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_STATTOOL.M with the given input arguments.
%
%      TOM_STATTOOL('Property','Value',...) creates a new TOM_STATTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_stattool_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_stattool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_stattool

% Last Modified by GUIDE v2.5 23-Nov-2006 13:51:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_stattool_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_stattool_OutputFcn, ...
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
function tom_stattool_OpeningFcn(hObject, eventdata, handles, varargin)


handles.data = varargin{1};
handles.datatype = varargin{2};


if isequal(handles.datatype,'2dstack')
    im = tom_emreadc(handles.data);
    handles.data = im.Value;
    handles.data = reshape(handles.data,[size(handles.data,2).*size(handles.data,1),size(handles.data,3)]);
elseif isequal(handles.datatype,'statmatrix')
    handles.data = varargin{1}';
end

if nargin < 3
    handles.mask_ind=find(ones(size(handles.data,1),1));
    
else
    handles.mask_ind= find(reshape(varargin{3},size(varargin{3},2).*size(varargin{3},1),1));
end;

handles = calc_vals(handles,handles.mask_ind);

% Choose default command line output for tom_stattool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_stattool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_stattool_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  checkbox mean                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_stattool_mean_Callback(hObject, eventdata, handles)

handles = replot(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  checkbox min                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_stattool_min_Callback(hObject, eventdata, handles)

handles = replot(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  checkbox max                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_stattool_max_Callback(hObject, eventdata, handles)

handles = replot(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  checkbox std                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_stattool_std_Callback(hObject, eventdata, handles)

handles = replot(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  checkbox var                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_stattool_var_Callback(hObject, eventdata, handles)

handles = replot(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display image                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_stattool_showimage_Callback(hObject, eventdata, handles)

binning = str2num(get(handles.edit_stattool_binning,'String'));

im = zeros(size(handles.data,1)./2^binning,size(handles.data,2));
for i=1:size(handles.data,2)
    im(:,i) = tom_bin(handles.data(:,i),binning);
end

axes(handles.stattool_display);
imagesc(im);colormap gray;

set(handles.checkbox_stattool_mean,'Value',0);
set(handles.checkbox_stattool_min,'Value',0);
set(handles.checkbox_stattool_max,'Value',0);
set(handles.checkbox_stattool_std,'Value',0);
set(handles.checkbox_stattool_var,'Value',0);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = calc_vals(handles,mask_ind)

numpoints = size(handles.data,2);

handles.mean = zeros(numpoints,1);
handles.min = zeros(numpoints,1);
handles.max = zeros(numpoints,1);
handles.std = zeros(numpoints,1);
handles.var = zeros(numpoints,1);

for i=1:numpoints   
    tmp=handles.data(:,i);
    [handles.mean(i), handles.max(i), handles.min(i), handles.std(i), handles.variance(i)] = tom_dev(tmp(mask_ind),'noinfo');
    %[handles.mean(i), handles.max(i), handles.min(i), handles.std(i), handles.variance(i)] = tom_dev(tmp,'noinfo');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  update plot                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = replot(handles)

axes(handles.stattool_display);
cla;reset(handles.stattool_display);
hold on;
if get(handles.checkbox_stattool_mean,'Value') == 1
    plot(handles.mean,'Color',[1 0 0]);
end

if get(handles.checkbox_stattool_min,'Value') == 1
    plot(handles.min,'Color',[0 1 0]);
end

if get(handles.checkbox_stattool_max,'Value') == 1
    plot(handles.max,'Color',[0 0 1]);
end

if get(handles.checkbox_stattool_std,'Value') == 1
    plot(handles.std,'Color',[0 1 1]);
end

if get(handles.checkbox_stattool_var,'Value') == 1
    plot(handles.variance,'Color',[1 1 0]);
end


hold off;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  unused callbacks / create functions                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_stattool_binning_Callback(hObject, eventdata, handles)
function edit_stattool_binning_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




