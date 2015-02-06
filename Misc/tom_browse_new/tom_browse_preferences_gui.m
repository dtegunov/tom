function varargout = tom_browse_preferences_gui(varargin)
% TOM_BROWSE_PREFERENCES_GUI M-file for tom_browse_preferences_gui.fig
%      TOM_BROWSE_PREFERENCES_GUI, by itself, creates a new TOM_BROWSE_PREFERENCES_GUI or raises the existing
%      singleton*.
%
%      H = TOM_BROWSE_PREFERENCES_GUI returns the handle to a new TOM_BROWSE_PREFERENCES_GUI or the handle to
%      the existing singleton*.
%
%      TOM_BROWSE_PREFERENCES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_BROWSE_PREFERENCES_GUI.M with the given input arguments.
%
%      TOM_BROWSE_PREFERENCES_GUI('Property','Value',...) creates a new TOM_BROWSE_PREFERENCES_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_browse_preferences_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_browse_preferences_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_browse_preferences_gui

% Last Modified by GUIDE v2.5 27-Aug-2010 16:40:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_browse_preferences_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_browse_preferences_gui_OutputFcn, ...
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


% --- Executes just before tom_browse_preferences_gui is made visible.
function tom_browse_preferences_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_browse_preferences_gui (see VARARGIN)

% Choose default command line output for tom_browse_preferences_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_browse_preferences_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_browse_preferences_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browse_for_chimera_path.
function browse_for_chimera_path_Callback(hObject, eventdata, handles)
% hObject    handle to browse_for_chimera_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ispref('tom_toolbox','tom_browse_chimera_path')
    p=getpref('tom_toolbox','tom_browse_chimera_path');
    disp(['Current path for chimera: ' p]);
else
    disp(['Path preferences don''t exist.']);
end;
[filename, pathname] = uigetfile('*.*', 'Pick the chimera Exe-file');
if isequal(filename,0) || isequal(pathname,0)
    disp('No file selected. No changes.');return;
else
    disp(['Selected file: ', fullfile(pathname, filename)])
end
set(handles.path_for_chimera_gui,'String',fullfile(pathname, filename));
if ispref('tom_toolbox','tom_browse_chimera_path')
    setpref('tom_toolbox','tom_browse_chimera_path',fullfile(pathname, filename))
else
    addpref('tom_toolbox','tom_browse_chimera_path',fullfile(pathname, filename))
end;


% --- Executes on button press in chimera_path_which.
function chimera_path_which_Callback(hObject, eventdata, handles)
% hObject    handle to chimera_path_which (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isunix
    [s,p] = unix('which chimera');
    if s==0
        if ispref('tom_toolbox','tom_browse_chimera_path')
            setpref('tom_toolbox','tom_browse_chimera_path',fullfile(p))
        else
            addpref('tom_toolbox','tom_browse_chimera_path',fullfile(p))
        end;
        set(handles.path_for_chimera_gui,'String',fullfile(p(1:end-1)));
        disp(['chimera path: ', fullfile(p(1:end-1))])
    else
        disp('Cannot find chimera');
    end;
else
    disp('For Linux only.');
end;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
all_fontsizes = cellstr(get(hObject,'String'));
box_fontsize=str2num(all_fontsizes{get(hObject,'Value')});
if ispref('tom_toolbox','tom_browse_box_fontsize')
    setpref('tom_toolbox','tom_browse_box_fontsize',box_fontsize);
else
    addpref('tom_toolbox','tom_browse_box_fontsize',box_fontsize);
end;
disp(['Box fontsize: ', num2str(box_fontsize) ' pixel'])
set(handles.test_text,'Fontsize',box_fontsize)

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

all_schemes = cellstr(get(hObject,'String'));
color_scheme=all_schemes{get(hObject,'Value')};
if ispref('tom_toolbox','tom_browse_color_scheme')
    setpref('tom_toolbox','tom_browse_color_scheme',color_scheme);
else
    addpref('tom_toolbox','tom_browse_color_scheme',color_scheme);
end;
disp(['Color scheme: ', color_scheme])
if strcmp(color_scheme,'white on black')
    set(handles.test_text,'ForegroundColor','white')
    set(handles.test_text,'BackgroundColor','black')
end;
if strcmp(color_scheme,'black on white')
    set(handles.test_text,'ForegroundColor','black')
    set(handles.test_text,'BackgroundColor','white')
end;
if strcmp(color_scheme,'white on green')
    set(handles.test_text,'ForegroundColor','white')
    set(handles.test_text,'BackgroundColor',[.3 .6 .3])
end;
% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
