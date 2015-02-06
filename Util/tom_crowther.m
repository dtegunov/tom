function varargout = tom_crowther(varargin)
% TOM_CROWTHER M-file for tom_crowther.fig
%      TOM_CROWTHER, by itself, creates a new TOM_CROWTHER or raises the existing
%      singleton*.
%
%      H = TOM_CROWTHER returns the handle to a new TOM_CROWTHER or the handle to
%      the existing singleton*.
%
%      TOM_CROWTHER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_CROWTHER.M with the given input arguments.
%
%      TOM_CROWTHER('Property','Value',...) creates a new TOM_CROWTHER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_crowther_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_crowther_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_crowther

% Last Modified by GUIDE v2.5 11-Aug-2003 15:58:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_crowther_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_crowther_OutputFcn, ...
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


% --- Executes just before tom_crowther is made visible.
function tom_crowther_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_crowther (see VARARGIN)

% Choose default command line output for tom_crowther
handles.output = hObject;


% UIWAIT makes tom_crowther wait for user response (see UIRESUME)
% uiwait(handles.figure1);
handles.objectthickness=200; % Object thickness in nm
handles.angularrange=140; % in deg
handles.numberofprojections=75; % number of projections
handles.finalresolution=2; % resolution in nm
handles.finalresolution=handles.objectthickness.*pi./(handles.numberofprojections.*180./handles.angularrange);
handles.increment=140./handles.numberofprojections; % in deg

set(handles.resolution,'String',num2str(handles.finalresolution));

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = tom_crowther_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function thickness_Callback(hObject, eventdata, handles)
% hObject    handle to thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thickness as text
%        str2double(get(hObject,'String')) returns contents of thickness as a double

handles.objectthickness=str2num(get(handles.thickness,'String'));
handles.finalresolution=handles.objectthickness.*pi./(handles.numberofprojections.*180./handles.angularrange);

set(handles.resolution,'String',num2str(handles.finalresolution));

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tiltincrement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tiltincrement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function tiltincrement_Callback(hObject, eventdata, handles)
% hObject    handle to tiltincrement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tiltincrement as text
%        str2double(get(hObject,'String')) returns contents of tiltincrement as a double
handles.increment=str2num(get(handles.tiltincrement,'String'));
handles.numberofprojections=180./handles.increment;
handles.finalresolution=handles.objectthickness.*pi./(handles.numberofprojections.*180./handles.angularrange);

set(handles.number,'String',num2str(handles.numberofprojections));
set(handles.tiltincrement,'String',num2str(handles.increment));
set(handles.resolution,'String',num2str(handles.finalresolution));

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function number_Callback(hObject, eventdata, handles)
% hObject    handle to number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number as text
%        str2double(get(hObject,'String')) returns contents of number as a double

handles.numberofprojections=str2num(get(handles.number,'String'));
handles.increment=180./handles.numberofprojections;
handles.finalresolution=handles.objectthickness.*pi./(handles.numberofprojections.*180./handles.angularrange);

set(handles.number,'String',num2str(handles.numberofprojections));
set(handles.tiltincrement,'String',num2str(handles.increment));
set(handles.resolution,'String',num2str(handles.finalresolution));

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function resolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function resolution_Callback(hObject, eventdata, handles)
% hObject    handle to resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resolution as text
%        str2double(get(hObject,'String')) returns contents of resolution as a double


handles.finalresolution=str2num(get(handles.resolution,'String'));
handles.numberofprojections=handles.objectthickness.*pi./handles.finalresolution;

handles.increment=180./(handles.numberofprojections.*180./handles.angularrange);


set(handles.number,'String',num2str(handles.numberofprojections));
set(handles.tiltincrement,'String',num2str(handles.increment));
set(handles.resolution,'String',num2str(handles.finalresolution));

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function range_Callback(hObject, eventdata, handles)
% hObject    handle to range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of range as text
%        str2double(get(hObject,'String')) returns contents of range as a double


