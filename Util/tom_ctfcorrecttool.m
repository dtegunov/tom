function varargout = tom_ctfcorrecttool(varargin)
% TOM_CTFCORRECTTOOL M-file for tom_ctfcorrecttool.fig
%      TOM_CTFCORRECTTOOL, by itself, creates a new TOM_CTFCORRECTTOOL or raises the existing
%      singleton*.
%
%      H = TOM_CTFCORRECTTOOL returns the handle to a new TOM_CTFCORRECTTOOL or the handle to
%      the existing singleton*.
%
%      TOM_CTFCORRECTTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_CTFCORRECTTOOL.M with the given input arguments.
%
%      TOM_CTFCORRECTTOOL('Property','Value',...) creates a new TOM_CTFCORRECTTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_ctfcorrecttool_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_ctfcorrecttool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_ctfcorrecttool

% Last Modified by GUIDE v2.5 10-Jan-2007 09:53:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_ctfcorrecttool_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_ctfcorrecttool_OutputFcn, ...
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


% --- Executes just before tom_ctfcorrecttool is made visible.
function tom_ctfcorrecttool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_ctfcorrecttool (see VARARGIN)

% Choose default command line output for tom_ctfcorrecttool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_ctfcorrecttool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_ctfcorrecttool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function input_ctfdeconvtool_file_Callback(hObject, eventdata, handles)
% hObject    handle to input_ctfdeconvtool_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_ctfdeconvtool_file as text
%        str2double(get(hObject,'String')) returns contents of input_ctfdeconvtool_file as a double


% --- Executes during object creation, after setting all properties.
function input_ctfdeconvtool_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_ctfdeconvtool_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_ctfdeconvtool_loadfile.
function button_ctfdeconvtool_loadfile_Callback(hObject, eventdata, handles)
% hObject    handle to button_ctfdeconvtool_loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function param_defocus_Callback(hObject, eventdata, handles)
% hObject    handle to param_defocus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of param_defocus as text
%        str2double(get(hObject,'String')) returns contents of param_defocus as a double


% --- Executes during object creation, after setting all properties.
function param_defocus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param_defocus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


