function varargout = tom_xmipp_psd_enhance_gui(varargin)


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_xmipp_psd_enhance_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_xmipp_psd_enhance_gui_OutputFcn, ...
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

% -------------------------------------------------------------------------
% Opening function
% -------------------------------------------------------------------------
function tom_xmipp_psd_enhance_gui_OpeningFcn(hObject, eventdata, handles, varargin)

set(handles.tom_xmipp_psd_enhance_gui,'Currentaxes',handles.psdaxes);
axis off;

set(handles.tom_xmipp_psd_enhance_gui,'Currentaxes',handles.imageaxes);
axis off;



handles.clim = [];

handles.output = hObject;
guidata(hObject, handles);


% -------------------------------------------------------------------------
% Output function
% -------------------------------------------------------------------------
function varargout = tom_xmipp_psd_enhance_gui_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% -------------------------------------------------------------------------
% edit load image
% -------------------------------------------------------------------------
function edit_loadimage_Callback(hObject, eventdata, handles)

fn = get(handles.edit_loadimage,'String');
handles = show_thumb(handles,fn);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% button browse
% -------------------------------------------------------------------------
function button_browse_Callback(hObject, eventdata, handles)

[FileName,PathName] = uigetfile({'*.em';'*.spi';'*.*'},'Select file to open');

fn = [ PathName '/' FileName];

handles = show_thumb(handles,fn);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% Load image
% -------------------------------------------------------------------------
function button_loadimage_Callback(hObject, eventdata, handles)

fn = get(handles.edit_loadimage,'String');

if tom_isemfile(fn)
    handles.image = tom_emreadc(fn);
elseif tom_isspiderfile(fn)
    handles.image = tom_spiderread(fn);
else
    errordlg('Only em and spider files supported.');
    return;
end

handles.image.Value = tom_calc_periodogram(handles.image.Value,128);

handles = enhance_psd(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% Show image
% -------------------------------------------------------------------------
function button_show_Callback(hObject, eventdata, handles)

fn = get(handles.edit_loadimage,'String');

if tom_isemfile(fn)
    image = tom_emreadc(fn);
elseif tom_isspiderfile(fn)
    image = tom_spiderread(fn);
else
    errordlg('Only em and spider files supported.');
    return;
end

imtool(image.Value);

guidata(hObject, handles);

% -------------------------------------------------------------------------
% edit parameter  filter w1
% -------------------------------------------------------------------------
function edit_param_w1_Callback(hObject, eventdata, handles)

handles = enhance_psd(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% edit parameter  filter w2
% -------------------------------------------------------------------------
function edit_param_w2_Callback(hObject, eventdata, handles)

handles = enhance_psd(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% edit parameter  mask m1
% -------------------------------------------------------------------------
function edit_param_m1_Callback(hObject, eventdata, handles)

handles = enhance_psd(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% edit parameter  mask m2
% -------------------------------------------------------------------------
function edit_param_m2_Callback(hObject, eventdata, handles)

handles = enhance_psd(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% edit parameter  mask m2
% -------------------------------------------------------------------------
function edit_param_decay_Callback(hObject, eventdata, handles)

handles = enhance_psd(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% button contrast
% -------------------------------------------------------------------------
function button_contrast_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    handles.contrasttool = imcontrast(handles.psdaxes);
else
    delete(handles.contrasttool);
end

guidata(hObject, handles);


% -------------------------------------------------------------------------
% render thumbnail
% -------------------------------------------------------------------------
function handles = show_thumb(handles,fn)

if tom_isemfile(fn)
    imagethumb = tom_emreadc(fn,'resample',[2 2 0]);
elseif tom_isspiderfile(fn)
    imagethumb = tom_binc(tom_spiderread(fn),2);
else
    errordlg('Only em and spider files supported.');
    return;
end

set(handles.edit_loadimage,'String',fn);

set(handles.tom_xmipp_psd_enhance_gui,'Currentaxes',handles.imageaxes);
imagesc(imagethumb.Value');axis ij;colormap gray;axis off;


% -------------------------------------------------------------------------
% enhance psd
% -------------------------------------------------------------------------
function handles = enhance_psd(handles)

f_w1 = str2double(get(handles.edit_param_w1,'String'));
f_w2 = str2double(get(handles.edit_param_w2,'String'));
decay = str2double(get(handles.edit_param_decay,'String'));
m1 = str2double(get(handles.edit_param_m1,'String'));
m2 = str2double(get(handles.edit_param_m2,'String'));

psd = tom_xmipp_psd_enhance(handles.image.Value,true,true,f_w1,f_w2,decay,m1,m2);

set(handles.tom_xmipp_psd_enhance_gui,'Currentaxes',handles.psdaxes);

try
    handles.clim = get(handles.psdaxes,'Clim');
end

imagesc(psd');axis ij;colormap gray;axis off;

if ~isempty(handles.clim)
    set(handles.psdaxes,'Clim',handles.clim);
end

if get(handles.button_contrast,'Value') == 1
    handles.contrasttool = imcontrast(handles.psdaxes);
end





% -------------------------------------------------------------------------
% Create Functions
% -------------------------------------------------------------------------
function edit_loadimage_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_param_m1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_param_decay_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_param_m2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_param_w2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_param_w1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





