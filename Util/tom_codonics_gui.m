function varargout = tom_codonics_gui(varargin)
%TOM_CODONICS_GUI M-file for tom_codonics_gui.fig
%
%   tom_amira_createisosurface(filelabel,threshold, label, color, transformmatrix, iconposition, host)
%
%      TOM_CODONICS_GUI, by itself, creates a new TOM_CODONICS_GUI or raises the existing
%      singleton*.
%
%      H = TOM_CODONICS_GUI returns the handle to a new TOM_CODONICS_GUI or the handle to
%      the existing singleton*.
%
%      TOM_CODONICS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_CODONICS_GUI.M with the given input arguments.
%
%      TOM_CODONICS_GUI('Property','Value',...) creates a new TOM_CODONICS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_codonics_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_codonics_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout           ...
%
%EXAMPLE
%   ... = tom_codonics_gui(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   GUIDE, GUIDATA, GUIHANDLES
%
% Edit the above text to modify the response to help tom_codonics_gui
%
% Last Modified by GUIDE v2.5 17-Jun-2004 12:16:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_codonics_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_codonics_gui_OutputFcn, ...
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


% --- Executes just before tom_codonics_gui is made visible.
function tom_codonics_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_codonics_gui (see VARARGIN)

% Choose default command line output for tom_codonics_gui
handles.output = hObject;
set(handles.overview,'Value',0);

set(handles.size_x,'Enable','on');
set(handles.size_y,'Enable','on');
set(handles.contrast_value,'Enable','on');
set(handles.gamma_value,'Enable','on');
set(handles.tcr_value,'Enable','on');
set(handles.mcm_value,'Enable','on');
set(handles.contrast_min,'Enable','off');
set(handles.contrast_max,'Enable','off');
set(handles.gamma_min,'Enable','off');
set(handles.gamma_max,'Enable','off');
set(handles.num_x,'Enable','off');
set(handles.num_y,'Enable','off');
set(handles.tcr_min,'Enable','off');
set(handles.tcr_max,'Enable','off');
set(handles.mcm_min,'Enable','off');
set(handles.mcm_max,'Enable','off');
set(handles.use_contrast,'Enable','off');
set(handles.use_gamma,'Enable','off');
set(handles.use_tcr,'Enable','off');
set(handles.use_mcr,'Enable','off');

set(handles.size_default,'Enable','on');
set(handles.contrast_default,'Enable','on');
set(handles.gamma_default,'Enable','on');
set(handles.tcr_default,'Enable','on');
set(handles.mcm_default,'Enable','on');

    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_codonics_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_codonics_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in overview.
function overview_Callback(hObject, eventdata, handles)
% hObject    handle to overview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of overview

state=get(handles.overview,'Value');

if (state==1)
    set(handles.size_x,'Enable','off');
    set(handles.size_y,'Enable','off');
    set(handles.contrast_value,'Enable','off');
    set(handles.gamma_value,'Enable','off');
    set(handles.tcr_value,'Enable','off');
    set(handles.mcm_value,'Enable','off');
    set(handles.contrast_min,'Enable','on');
    set(handles.contrast_max,'Enable','on');
    set(handles.gamma_min,'Enable','on');
    set(handles.gamma_max,'Enable','on');
    set(handles.num_x,'Enable','on');
    set(handles.num_y,'Enable','on');
    set(handles.tcr_min,'Enable','on');
    set(handles.tcr_max,'Enable','on');
    set(handles.mcm_min,'Enable','on');
    set(handles.mcm_max,'Enable','on');
    set(handles.use_contrast,'Enable','on');
    set(handles.use_gamma,'Enable','on');
    set(handles.use_tcr,'Enable','on');
    set(handles.use_mcr,'Enable','on');
    
    set(handles.size_default,'Enable','off');
    set(handles.contrast_default,'Enable','off');
    set(handles.gamma_default,'Enable','off');
    set(handles.tcr_default,'Enable','off');
    set(handles.mcm_default,'Enable','off');
   
    
    
else
    set(handles.size_x,'Enable','on');
    set(handles.size_y,'Enable','on');
    set(handles.contrast_value,'Enable','on');
    set(handles.gamma_value,'Enable','on');
    set(handles.tcr_value,'Enable','on');
    set(handles.mcm_value,'Enable','on');
    set(handles.contrast_min,'Enable','off');
    set(handles.contrast_max,'Enable','off');
    set(handles.gamma_min,'Enable','off');
    set(handles.gamma_max,'Enable','off');
    set(handles.num_x,'Enable','off');
    set(handles.num_y,'Enable','off');
    set(handles.tcr_min,'Enable','off');
    set(handles.tcr_max,'Enable','off');
    set(handles.mcm_min,'Enable','off');
    set(handles.mcm_max,'Enable','off');
    set(handles.use_contrast,'Enable','off');
    set(handles.use_gamma,'Enable','off');
    set(handles.use_tcr,'Enable','off');
    set(handles.use_mcr,'Enable','off');
    
    set(handles.size_default,'Enable','on');
    set(handles.contrast_default,'Enable','on');
    set(handles.gamma_default,'Enable','on');
    set(handles.tcr_default,'Enable','on');
    set(handles.mcm_default,'Enable','on');
    
end;


% --- Executes on button press in go.
function go_Callback(hObject, eventdata, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
state=get(handles.overview,'Value');

path=get(handles.browse_txt,'String');

if (state==1)
    
    use=[ get(handles.use_contrast,'Value') get(handles.use_gamma,'Value') ...
            get(handles.use_tcr,'Value') get(handles.use_mcr,'Value') ];
    
    kk=1;
    if (use(1)==1)
        range(kk,1)=str2num(get(handles.contrast_min,'String'));
        range(kk,2)=str2num(get(handles.contrast_max,'String'));
        kk=kk+1;
    end;
    
    if (use(2)==1)
        range(kk,1)=str2num(get(handles.gamma_min,'String'));
        range(kk,2)=str2num(get(handles.gamma_max,'String'));
        kk=kk+1;
    end;
    
    if (use(3)==1)
        range(kk,1)=str2num(get(handles.tcr_min,'String'));
        range(kk,2)=str2num(get(handles.tcr_max,'String'));
        kk=kk+1;
    end;
    
    if (use(4)==1)
        range(kk,1)=str2num(get(handles.mcm_min,'String'));
        range(kk,2)=str2num(get(handles.mcm_max,'String'));
        kk=kk+1;
    end;
    num(1)=str2num(get(handles.num_x,'String'));
    num(2)=str2num(get(handles.num_y,'String'));
    m=tom_codonics_overview(path,range(1,:),range(2,:),num,use);
    m
else
    default=[get(handles.size_default,'Value') get(handles.contrast_default,'Value') ...
            get(handles.gamma_default,'Value') get(handles.tcr_default,'Value') get(handles.mcm_default,'Value')];
    size(1)=str2num(get(handles.size_x,'String'));
    size(2)=str2num(get(handles.size_y,'String'));
    contrast=str2num(get(handles.contrast_value,'String'));
    gamma=str2num(get(handles.gamma_value,'String'));
    tcr=str2num(get(handles.tcr_value,'String'));
    mcm=str2num(get(handles.mcm_value,'String'));
    
    tom_codonics(path,gamma,contrast,tcr,mcm,size,default);
end;



% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[myfile mypath] = uigetfile('*.tif');
mypath=[mypath myfile];
set(handles.browse_txt,'String',mypath);





% --- Executes during object creation, after setting all properties.
function contrast_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contrast_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





% --- Executes during object creation, after setting all properties.
function gamma_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





% --- Executes during object creation, after setting all properties.
function contrast_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contrast_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes during object creation, after setting all properties.
function gamma_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





% --- Executes during object creation, after setting all properties.
function contrast_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contrast_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





% --- Executes during object creation, after setting all properties.
function gamma_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes during object creation, after setting all properties.
function browse_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to browse_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





% --- Executes during object creation, after setting all properties.
function size_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to size_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





% --- Executes during object creation, after setting all properties.
function size_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to size_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function num_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function num_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function tcr_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tcr_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function tcr_min_Callback(hObject, eventdata, handles)
% hObject    handle to tcr_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tcr_min as text
%        str2double(get(hObject,'String')) returns contents of tcr_min as a double


% --- Executes during object creation, after setting all properties.
function mcm_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mcm_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function mcm_min_Callback(hObject, eventdata, handles)
% hObject    handle to mcm_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mcm_min as text
%        str2double(get(hObject,'String')) returns contents of mcm_min as a double


% --- Executes during object creation, after setting all properties.
function tcr_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tcr_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function tcr_max_Callback(hObject, eventdata, handles)
% hObject    handle to tcr_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tcr_max as text
%        str2double(get(hObject,'String')) returns contents of tcr_max as a double


% --- Executes during object creation, after setting all properties.
function mcm_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mcm_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function mcm_max_Callback(hObject, eventdata, handles)
% hObject    handle to mcm_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mcm_max as text
%        str2double(get(hObject,'String')) returns contents of mcm_max as a double


% --- Executes during object creation, after setting all properties.
function tcr_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tcr_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function tcr_value_Callback(hObject, eventdata, handles)
% hObject    handle to tcr_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tcr_value as text
%        str2double(get(hObject,'String')) returns contents of tcr_value as a double


% --- Executes during object creation, after setting all properties.
function mcm_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mcm_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function mcm_value_Callback(hObject, eventdata, handles)
% hObject    handle to mcm_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mcm_value as text
%        str2double(get(hObject,'String')) returns contents of mcm_value as a double


% --- Executes on button press in use_contrast.
function use_contrast_Callback(hObject, eventdata, handles)
% hObject    handle to use_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_contrast


% --- Executes on button press in use_gamma.
function use_gamma_Callback(hObject, eventdata, handles)
% hObject    handle to use_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_gamma


% --- Executes on button press in use_tcr.
function use_tcr_Callback(hObject, eventdata, handles)
% hObject    handle to use_tcr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_tcr


% --- Executes on button press in use_mcr.
function use_mcr_Callback(hObject, eventdata, handles)
% hObject    handle to use_mcr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_mcr


% --- Executes on button press in contrast_default.
function contrast_default_Callback(hObject, eventdata, handles)
% hObject    handle to contrast_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of contrast_default


% --- Executes on button press in gamma_default.
function gamma_default_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gamma_default


% --- Executes on button press in tcr_default.
function tcr_default_Callback(hObject, eventdata, handles)
% hObject    handle to tcr_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tcr_default


% --- Executes on button press in mcm_default.
function mcm_default_Callback(hObject, eventdata, handles)
% hObject    handle to mcm_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mcm_default


% --- Executes on button press in size_default.
function size_default_Callback(hObject, eventdata, handles)
% hObject    handle to size_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of size_default


