function varargout = tom_makemoviegui(varargin)
%TOM_MAKEMOVIEGUI by itself, creates a new TOM_MAKEMOVIEGUI or raises the existing
%      singleton*.
%
%   varargout = tom_makemoviegui(varargin)
%
%   H = TOM_MAKEMOVIEGUI returns the handle to a new TOM_MAKEMOVIEGUI or the handle to
%   the existing singleton*.
%
%   TOM_MAKEMOVIEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%   function named CALLBACK in TOM_MAKEMOVIEGUI.M with the given input arguments.
%
%   TOM_MAKEMOVIEGUI('Property','Value',...) creates a new TOM_MAKEMOVIEGUI or raises the
%   existing singleton*.  Starting from the left, property value pairs are
%   applied to the GUI before tom_makemoviegui_OpeningFunction gets called.  An
%   unrecognized property name or invalid value makes property application
%   stop.  All inputs are passed to tom_makemoviegui_OpeningFcn via varargin.
%
%   *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%   instance to run (singleton)".
%
%PARAMETERS
%
%  INPUT
%   filelabel           ...
%   threshold           ...
%   label               ...
%   color               ...
%   transformmatrix     ...
%   iconposition        ...
%   host                ...
%  
%  OUTPUT
%   data		...
%
%EXAMPLE
%   tom_makemoviegui(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   GUIDE, GUIDATA, GUIHANDLES
%
% Copyright 2002-2003 The MathWorks, Inc.
%
% Edit the above text to modify the response to help tom_makemoviegui
%
% Last Modified by GUIDE v2.5 28-Sep-2005 10:59:39


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_makemoviegui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_makemoviegui_OutputFcn, ...
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


% --- Executes just before tom_makemoviegui is made visible.
function tom_makemoviegui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_makemoviegui (see VARARGIN)

% Choose default command line output for tom_makemoviegui
handles.output = hObject;

handles.xyz = varargin{1};

tmpobj = findobj('Tag','frame_start');
set(tmpobj,'String','1');

tmpobj = findobj('Tag','frame_end');
set(tmpobj,'String',num2str(handles.xyz.dimensions.z));

tmpobj = findobj('Tag','direction');
set(tmpobj,'Value',3);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_makemoviegui wait for user response (see UIRESUME)
% uiwait(handles.moviemaker);


% --- Outputs from this function are returned to the command line.
function varargout = tom_makemoviegui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function size_x_Callback(hObject, eventdata, handles)
% hObject    handle to size_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of size_x as text
%        str2double(get(hObject,'String')) returns contents of size_x as a double


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



function size_y_Callback(hObject, eventdata, handles)
% hObject    handle to size_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of size_y as text
%        str2double(get(hObject,'String')) returns contents of size_y as a double


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


% --- Executes on selection change in direction.
function direction_Callback(hObject, eventdata, handles)
% hObject    handle to direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns direction contents as cell array
%        contents{get(hObject,'Value')} returns selected item from direction

tmpobj = findobj('Tag','frame_end');
if get(hObject,'Value') == 1
	set(tmpobj,'String',num2str(handles.xyz.dimensions.x));
elseif get(hObject,'Value') == 2
	set(tmpobj,'String',num2str(handles.xyz.dimensions.y));
else
	set(tmpobj,'String',num2str(handles.xyz.dimensions.z));
end

% --- Executes during object creation, after setting all properties.
function direction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function frame_start_Callback(hObject, eventdata, handles)
% hObject    handle to frame_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_start as text
%        str2double(get(hObject,'String')) returns contents of frame_start as a double


% --- Executes during object creation, after setting all properties.
function frame_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function frame_end_Callback(hObject, eventdata, handles)
% hObject    handle to frame_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_end as text
%        str2double(get(hObject,'String')) returns contents of frame_end as a double


% --- Executes during object creation, after setting all properties.
function frame_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function quality_Callback(hObject, eventdata, handles)
% hObject    handle to quality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of quality as text
%        str2double(get(hObject,'String')) returns contents of quality as a double


% --- Executes during object creation, after setting all properties.
function quality_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function fps_Callback(hObject, eventdata, handles)
% hObject    handle to fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fps as text
%        str2double(get(hObject,'String')) returns contents of fps as a double


% --- Executes during object creation, after setting all properties.
function fps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function filename_movie_Callback(hObject, eventdata, handles)
% hObject    handle to filename_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename_movie as text
%        str2double(get(hObject,'String')) returns contents of filename_movie as a double


% --- Executes during object creation, after setting all properties.
function filename_movie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[File Path] = uiputfile('*.avi');
tmpobj = findobj('Tag','filename_movie');
set(tmpobj,'String',strcat(Path,File));


% --- Executes on button press in button_makemovie.
function button_makemovie_Callback(hObject, eventdata, handles)
% hObject    handle to button_makemovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

in = handles.xyz.Volume.Value;
out = get(findobj('Tag','filename_movie'),'String');
fps = str2num(get(findobj('Tag','fps'),'String'));
quality = str2num(get(findobj('Tag','quality'),'String'));
contrast = handles.xyz.DataScale;
average.xy = handles.xyz.average.x;
average.xz = handles.xyz.average.y;
average.yz = handles.xyz.average.z;
zoom.xy = handles.xyz.actualaxis.xy;
zoom.xz = handles.xyz.actualaxis.xz;
zoom.yz = handles.xyz.actualaxis.yz;
frames = [str2num(get(findobj('Tag','frame_start'),'String')) str2num(get(findobj('Tag','frame_end'),'String'))];
direction = get(findobj('Tag','direction'),'Value');
if direction == 1
	direction = 'x';
elseif direction == 2
	direction = 'y';
else
	direction = 'z';
end 
position.xy = handles.xyz.position.xy;
position.xz = handles.xyz.position.xz;
position.yz = handles.xyz.position.yz;
moviesize = [str2num(get(findobj('Tag','size_x'),'String')) str2num(get(findobj('Tag','size_y'),'String'))];
bandpass.filter_low = handles.xyz.filter.low;
bandpass.filter_high = handles.xyz.filter.high;
bandpass.filter_xy = handles.xyz.filter.xy;
bandpass.filter_xz = handles.xyz.filter.xz;
bandpass.filter_yz = handles.xyz.filter.yz;
linewidth = [str2num(get(findobj('Tag','outerwidth'),'String')) str2num(get(findobj('Tag','innerwidth'),'String'))];
redlinewidth = str2num(get(findobj('Tag','redwidth'),'String'));
numbersflag = get(handles.slicenumbers,'Value');


%contrast = [-.1418,0.124891];
tom_makemovie3d(in,out,fps,quality,contrast,average,zoom,frames,direction,position,moviesize,bandpass,linewidth,redlinewidth,numbersflag);

close(findobj('Tag','moviemaker'));



function outerwidth_Callback(hObject, eventdata, handles)
% hObject    handle to outerwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outerwidth as text
%        str2double(get(hObject,'String')) returns contents of outerwidth as a double


% --- Executes during object creation, after setting all properties.
function outerwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outerwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function innerwidth_Callback(hObject, eventdata, handles)
% hObject    handle to innerwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of innerwidth as text
%        str2double(get(hObject,'String')) returns contents of innerwidth as a double


% --- Executes during object creation, after setting all properties.
function innerwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to innerwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in slicenumbers.
function slicenumbers_Callback(hObject, eventdata, handles)
% hObject    handle to slicenumbers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slicenumbers



function redwidth_Callback(hObject, eventdata, handles)
% hObject    handle to redwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of redwidth as text
%        str2double(get(hObject,'String')) returns contents of redwidth as a double


% --- Executes during object creation, after setting all properties.
function redwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to redwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


