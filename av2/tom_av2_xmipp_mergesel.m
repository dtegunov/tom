function varargout = tom_av2_xmipp_mergesel(varargin)
% TOM_AV2_XMIPP_MERGESEL M-file for tom_av2_xmipp_mergesel.fig
%      TOM_AV2_XMIPP_MERGESEL, by itself, creates a new TOM_AV2_XMIPP_MERGESEL or raises the existing
%      singleton*.
%
%      H = TOM_AV2_XMIPP_MERGESEL returns the handle to a new TOM_AV2_XMIPP_MERGESEL or the handle to
%      the existing singleton*.
%
%      TOM_AV2_XMIPP_MERGESEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_AV2_XMIPP_MERGESEL.M with the given input arguments.
%
%      TOM_AV2_XMIPP_MERGESEL('Property','Value',...) creates a new TOM_AV2_XMIPP_MERGESEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_av2_xmipp_mergesel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_av2_xmipp_mergesel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_av2_xmipp_mergesel

% Last Modified by GUIDE v2.5 13-Aug-2009 17:33:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_av2_xmipp_mergesel_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_xmipp_mergesel_OutputFcn, ...
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


% --- Executes just before tom_av2_xmipp_mergesel is made visible.
function tom_av2_xmipp_mergesel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_av2_xmipp_mergesel (see VARARGIN)

% Choose default command line output for tom_av2_xmipp_mergesel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_av2_xmipp_mergesel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_av2_xmipp_mergesel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browse_selfile.
function browse_selfile_Callback(hObject, eventdata, handles)
% hObject    handle to browse_selfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 str=get(handles.listbox_sel,'String');

[file path] = uigetfile('*.sel');

if (isempty(str))
    new=[path file];
else
    if (ischar(str))
        new{1}=str;
        new{2}=[path file];
    else
        new=cat(1,str,[path file]);
    end;
end;

set(handles.listbox_sel,'String',new);


% --- Executes on selection change in listbox_sel.
function listbox_sel_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_sel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_sel


% --- Executes during object creation, after setting all properties.
function listbox_sel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_selfile.
function save_selfile_Callback(hObject, eventdata, handles)
% hObject    handle to save_selfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file path] = uiputfile('*.sel');

set(handles.output_path,'String',[path file]);


function output_path_Callback(hObject, eventdata, handles)
% hObject    handle to output_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_path as text
%        str2double(get(hObject,'String')) returns contents of output_path as a double


% --- Executes during object creation, after setting all properties.
function output_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out=get(handles.output_path,'String');
str = get(handles.listbox_sel,'String');
new = 'cat ';

for i=1:length(str)
    
    new = [new ' ' str{i}];

end


new=[new ' > ' out];

unix(new);






