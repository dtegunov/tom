function varargout = tom_av2_compress_stack_gui(varargin)
% TOM_AV2_COMPRESS_STACK_GUI MATLAB code for tom_av2_compress_stack_gui.fig
%      TOM_AV2_COMPRESS_STACK_GUI, by itself, creates a new TOM_AV2_COMPRESS_STACK_GUI or raises the existing
%      singleton*.
%
%      H = TOM_AV2_COMPRESS_STACK_GUI returns the handle to a new TOM_AV2_COMPRESS_STACK_GUI or the handle to
%      the existing singleton*.
%
%      TOM_AV2_COMPRESS_STACK_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_AV2_COMPRESS_STACK_GUI.M with the given input arguments.
%
%      TOM_AV2_COMPRESS_STACK_GUI('Property','Value',...) creates a new TOM_AV2_COMPRESS_STACK_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_av2_compress_stack_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_av2_compress_stack_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_av2_compress_stack_gui

% Last Modified by GUIDE v2.5 12-May-2011 12:18:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_av2_compress_stack_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_compress_stack_gui_OutputFcn, ...
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


% --- Executes just before tom_av2_compress_stack_gui is made visible.
function tom_av2_compress_stack_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_av2_compress_stack_gui (see VARARGIN)

% Choose default command line output for tom_av2_compress_stack_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_av2_compress_stack_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%defaults
st.input_str='';
st.output_str='';
st.e_num_of_classes=1000;
st.e_equal='not_equal';
st.e_pre_alg_th=0.2;
st.e_eigs='1:15';
st.e_cluster_method='k-means';
st.e_verbose=0;
st.e_bin=1;



if (nargin > 0)
    in_st=varargin{1};
    
    if (isfield(in_st,'input_str'))
        st.input_str=in_st.input_str;
    end;
    
    if (isfield(in_st,'output_str'))
        st.output_str=in_st.output_str;
    end;
    
    if (isfield(in_st,'e_num_of_classes'))
        st.e_num_of_classes=in_st.e_num_of_classes;
    end;
    
    if (isfield(in_st,'e_pre_alg_th'))
        st.e_pre_alg_th=in_st.e_pre_alg_th;
    end;
    
    if (isfield(in_st,'e_eigs'))
        st.e_eigs=in_st.e_eigs;
    end;
    
    if (isfield(in_st,'e_bin'))
        st.e_bin=in_st.e_bin;
    end;
    
    if (isfield(in_st,'e_cluster_method'))
        st.e_cluster_method=in_st.e_cluster_method;
    end;
    
    if (isfield(in_st,'e_verbose'))
        st.e_verbose=in_st.e_verbose;
    end;
    
    set_gui_values(handles,st);
    
end;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_os3_align_images wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_av2_compress_stack_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;

close(handles.figure1);

% --- Executes on button press in br_input.
function br_input_Callback(hObject, eventdata, handles)
% hObject    handle to br_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[sfile spath]=uigetfile('*.em');
set(handles.input_str,'String',[spath '/' sfile]);


function input_str_Callback(hObject, eventdata, handles)
% hObject    handle to input_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_str as text
%        str2double(get(hObject,'String')) returns contents of input_str as a double


% --- Executes during object creation, after setting all properties.
function input_str_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function e_num_of_classes_Callback(hObject, eventdata, handles)
% hObject    handle to e_num_of_classes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e_num_of_classes as text
%        str2double(get(hObject,'String')) returns contents of e_num_of_classes as a double


% --- Executes during object creation, after setting all properties.
function e_num_of_classes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_num_of_classes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function e_equal_Callback(hObject, eventdata, handles)
% hObject    handle to e_equal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e_equal as text
%        str2double(get(hObject,'String')) returns contents of e_equal as a double


% --- Executes during object creation, after setting all properties.
function e_equal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_equal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function e_pre_alg_th_Callback(hObject, eventdata, handles)
% hObject    handle to e_pre_alg_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e_pre_alg_th as text
%        str2double(get(hObject,'String')) returns contents of e_pre_alg_th as a double


% --- Executes during object creation, after setting all properties.
function e_pre_alg_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_pre_alg_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function e_eigs_Callback(hObject, eventdata, handles)
% hObject    handle to e_eigs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e_eigs as text
%        str2double(get(hObject,'String')) returns contents of e_eigs as a double


% --- Executes during object creation, after setting all properties.
function e_eigs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_eigs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function e_cluster_method_Callback(hObject, eventdata, handles)
% hObject    handle to e_cluster_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e_cluster_method as text
%        str2double(get(hObject,'String')) returns contents of e_cluster_method as a double


% --- Executes during object creation, after setting all properties.
function e_cluster_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_cluster_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function e_verbose_Callback(hObject, eventdata, handles)
% hObject    handle to e_verbose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e_verbose as text
%        str2double(get(hObject,'String')) returns contents of e_verbose as a double


% --- Executes during object creation, after setting all properties.
function e_verbose_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_verbose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in br_output.
function br_output_Callback(hObject, eventdata, handles)
% hObject    handle to br_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[sfile spath]=uiputfile;
set(handles.output_str,'String',[spath '/' sfile]);


function output_str_Callback(hObject, eventdata, handles)
% hObject    handle to output_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_str as text
%        str2double(get(hObject,'String')) returns contents of output_str as a double


% --- Executes during object creation, after setting all properties.
function output_str_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in go.
function go_Callback(hObject, eventdata, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);


handles.output=st;

guidata(hObject, handles);


uiresume(handles.figure1);
function st=get_gui_values(handles)

st.input_str=get(handles.input_str,'String');
st.output_str=get(handles.output_str,'String');
st.e_num_of_classes=str2double(get(handles.e_num_of_classes,'String'));
st.e_equal=get(handles.e_equal,'String');
st.e_pre_alg_th=str2double(get(handles.e_pre_alg_th,'String'));
st.e_eigs=get(handles.e_eigs,'String');
st.e_cluster_method=get(handles.e_cluster_method,'String');
st.e_verbose=str2double(get(handles.e_verbose,'String'));
st.e_bin=str2double(get(handles.e_bin,'String'));

function set_gui_values(handles,st)

set(handles.input_str,'String',st.input_str);
set(handles.output_str,'String',st.output_str);
set(handles.e_num_of_classes,'String',num2str(st.e_num_of_classes));
set(handles.e_equal,'String',st.e_equal);
set(handles.e_pre_alg_th,'String',num2str(st.e_pre_alg_th));
set(handles.e_eigs,'String',st.e_eigs);
set(handles.e_cluster_method,'String',st.e_cluster_method);
set(handles.e_verbose,'String',num2str(st.e_verbose));
set(handles.e_bin,'String',num2str(st.e_bin));


function e_bin_Callback(hObject, eventdata, handles)
% hObject    handle to e_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e_bin as text
%        str2double(get(hObject,'String')) returns contents of e_bin as a double


% --- Executes during object creation, after setting all properties.
function e_bin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
