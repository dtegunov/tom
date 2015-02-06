function varargout = tom_isosurface(varargin)
% TOM_ISOSURFACE M-file for tom_isosurface.fig
%      TOM_ISOSURFACE, by itself, creates a new TOM_ISOSURFACE or raises the existing
%      singleton*.
%
%      H = TOM_ISOSURFACE returns the handle to a new TOM_ISOSURFACE or the handle to
%      the existing singleton*.
%
%      TOM_ISOSURFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_ISOSURFACE.M with the given input arguments.
%
%      TOM_ISOSURFACE('Property','Value',...) creates a new TOM_ISOSURFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_isosurface_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_isosurface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_isosurface

% Last Modified by GUIDE v2.5 09-Jan-2003 14:02:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_isosurface_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_isosurface_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before tom_isosurface is made visible.
function tom_isosurface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_isosurface (see VARARGIN)

% Choose default command line output for tom_isosurface
handles.output = hObject;
% Update handles structure
if size(varargin,1) <1 
[filename, pathname] = uigetfile({'*.vol';'*.em';'*.*'}, 'Pick an EM-file');
    if isequal(filename,0) | isequal(pathname,0) error('No data loaded.');
        return; 
    end;
    t=tom_emread([pathname,filename]);
    handles.volume=t.Value;
else
handles.volume=cell2mat(varargin(1));
end;

handles.Resmap='';

if (nargin > 4)
    handles.Resmap=varargin{2};
end;
clear varargin;
handles.Histogram=tom_hist3d(handles.volume);
[handles.Histogram_mean handles.Histogram_max handles.Histogram_min handles.Histogram_std]=tom_dev(handles.Histogram);
[handles.mean handles.max handles.min handles.std]=tom_dev(handles.volume);
handles.molmass=0;
handles.pixelsize=5.0;
% Update handles structure
tmp_obj=findobj('Tag','Surface');
axes(tmp_obj);cla;axis vis3d;axis off;
set(gcf,'Renderer','OpenGL');
view(3); camlight; light;
lighting phong;


guidata(hObject, handles);

% UIWAIT makes tom_isosurface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_isosurface_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp_obj=findobj('Tag','thresh');

set(tmp_obj,'Min',handles.min);
set(tmp_obj,'Max',handles.max);
set(tmp_obj,'SliderStep',[1./1000 1./100]);
tmp_obj=findobj('Tag','Histogram');
axes(tmp_obj);
plot(handles.Histogram);
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function thresh_Callback(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.thresh=(get(hObject,'Value'));
tmp_obj=findobj('Tag','thresh_txt');
set(tmp_obj,'String',num2str(handles.thresh));
update_surface(hObject, eventdata, handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function thresh_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function thresh_txt_Callback(hObject, eventdata, handles)
% hObject    handle to thresh_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh_txt as text
%        str2double(get(hObject,'String')) returns contents of thresh_txt as a double
tmp_obj=findobj('Tag','thresh_txt');
get(tmp_obj,'Value');
slicexy=(eval(get(tmp_obj,'String')));
if slicexy<handles.min slicexy=handles.min ; set(tmp_obj,'String',slicexy); end;
if slicexy>handles.max slicexy=handles.max; set(tmp_obj,'String',slicexy); end;
tmp_obj=findobj('Tag','thresh');
set(tmp_obj,'Value',slicexy);
handles.thresh=slicexy;
update_surface(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function mass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function mass_Callback(hObject, eventdata, handles)
% hObject    handle to mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mass as text
%        str2double(get(hObject,'String')) returns contents of mass as a double

tmp_obj=findobj('Tag','mass');
get(tmp_obj,'Value');
handles.molmass=(eval(get(tmp_obj,'String')));
nrmass=sum(sum(sum(handles.volume>=handles.max)));
rho=1.3; % g/cm^3 for protein
pixelsize_cm3=((handles.pixelsize.*1.0e-10).^3)./(0.01.^3);

mass_g=nrmass.*rho.*pixelsize_cm3;
mass_kg=mass_g./(10.^3);
mass=mass_kg./(1.66e-27)./1000;
lauf=handles.max;

% coarse search
while mass<handles.molmass
    lauf=lauf-.1;
    nrmass=sum(sum(sum(handles.volume>=lauf)));
    handles.mass=nrmass;
    rho=1.3; % g/cm^3 for protein
    pixelsize_cm3=(handles.pixelsize.*1.0e-10).^3./(0.01.^3);
    nr=nrmass;
    mass_g=nr.*rho.*pixelsize_cm3;
    mass_kg=mass_g./(10.^3);
    mass=mass_kg./(1.66e-27)./1000;
end;
handles.thresh=lauf+.1;
guidata(hObject, handles);
update_txt(hObject, eventdata, handles);
update_surface(hObject, eventdata, handles);
%guidata(hObject, handles);

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function pixelsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function pixelsize_Callback(hObject, eventdata, handles)
% hObject    handle to pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixelsize as text
%        str2double(get(hObject,'String')) returns contents of pixelsize as a double

tmp_obj=findobj('Tag','pixelsize');
get(tmp_obj,'Value');
handles.pixelsize=(eval(get(tmp_obj,'String')));
guidata(hObject, handles);

function update_surface(hObject, eventdata, handles)
tmp_obj=findobj('Tag','Surface');
axes(tmp_obj);cla;axis vis3d;axis off;

if (isempty(handles.Resmap))
    p = patch(isosurface(handles.volume, handles.thresh),'FaceColor', [230./255 213./255 17./255], 'EdgeColor', 'none');
    isonormals(handles.volume,p);
    set(gcf,'Renderer','OpenGL'); 
else
    %p = patch(isosurface(handles.volume, handles.thresh,handles.Resmap),'FaceColor', [230./255 213./255 17./255], 'EdgeColor', 'none');
    isosurface(handles.volume, handles.thresh,handles.Resmap);
end;


view(3); camlight; light; material shiny; camlight (270, 270,'infinite');
lighting phong; rotate3d;
handles.mass=sum(sum(sum(handles.volume>=handles.thresh)));
guidata(hObject, handles);
voxel_in_molmass(hObject, eventdata, handles);
guidata(hObject, handles);

function update_txt(hObject, eventdata, handles)
tmp_obj=findobj('Tag','thresh_txt');
set(tmp_obj,'String',num2str(handles.thresh));
tmp_obj=findobj('Tag','thresh');
set(tmp_obj,'Value',handles.thresh);
tmp_obj=findobj('Tag','mass');
set(tmp_obj,'Value',handles.molmass);

function voxel_in_molmass(hObject, eventdata, handles)
tmp_obj=findobj('Tag','pixelsize');
get(tmp_obj,'Value');
handles.pixelsize=(eval(get(tmp_obj,'String')));
rho=1.3; % g/cm^3 for protein
pixelsize_cm3=(handles.pixelsize.*1.0e-10).^3./(0.01.^3);
nr=handles.mass;
mass_g=nr.*rho.*pixelsize_cm3;
mass_kg=mass_g./(10.^3);
handles.molmass=mass_kg./(1.66e-27)./1000;
tmp_obj=findobj('Tag','mass');
set(tmp_obj,'String',num2str(handles.molmass));
guidata(hObject, handles);
