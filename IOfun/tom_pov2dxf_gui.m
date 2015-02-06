function varargout = tom_pov2dxf_gui(varargin)
% TOM_POV2DXF_GUI is the user interface for tom_pov2dxf which is a format
% converter for 3D graphics.
% tom_pov2dxf converts a .pov file created by chimera into a .dxf file
% which can be used by Bryce. TOM_POV2DXF can handle large files and allows
% surface reduction.
%
%   tom_pov2dxf(filename,dxf_filename,reduce_patch,matlab_surface_file)
%
% PARAMETERS
%   INPUT
%   filename        name of the .pov file
%   dxf_filename    name of the .dxf file 
%   reduce_patch    reduces the number of triangles (0<reduce_patch<=1).
%                   1.0 is 100 percent (no reduction).
%                   A number > 1 defines the number of faces which will be created.
%   matlab_surface_file Matlab .mat filename with surface object (optional).
%
%   OUTPUT
%   total_faces     total number of faces created
%
% EXAMPLE
%   tom_pov2dxf('lid_sym_atpase.pov','lid_sym_atpase',.5,'surface.mat');
%   load surface.mat
%   patch('Faces',surface.faces,'Vertices',surface.vertices,'FaceVertexCData'...
%   ,ones(size(surface.faces,1),3).*.8,'FaceColor','flat');
%   view(3);grid on;
%
%REFERENCES
%
% SEE ALSO
%   tom_pov2dxf, tom_surface2dxf, tom_surface2stl, tom_stl2surface, patch,
%   reducepatch
%
%   last change
%   07/13/08 SN
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2008
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_pov2dxf_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_pov2dxf_gui_OutputFcn, ...
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

% --- Executes just before tom_pov2dxf_gui is made visible.
function tom_pov2dxf_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_pov2dxf_gui (see VARARGIN)

% Choose default command line output for tom_pov2dxf_gui
handles.output = hObject;
handles.number_of_faces=0.5;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes tom_pov2dxf_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_pov2dxf_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function pov_filename_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pov_filename_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.pov_filename='';
set(hObject,'String',handles.pov_filename);
handles.pov_filename=get(hObject,'String');
guidata(hObject,handles);


function pov_filename_display_Callback(hObject, eventdata, handles)
% hObject    handle to pov_filename_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pov_filename_display as text
%        str2double(get(hObject,'String')) returns contents of pov_filename_display as a double
handles.pov_filename=get(hObject,'String');
w=what;
handles.pov_filename=fullfile(w.path,handles.pov_filename);
if exist(handles.pov_filename)
    set(hObject,'String',handles.pov_filename);
else
    handles.pov_filename='POV_Filename';
    set(hObject,'String',handles.pov_filename);
    disp('File doesn''t exist');
end;

% Save the new pov_filename_display value
guidata(hObject,handles);

% --- Executes on button press in convert.
function convert_Callback(hObject, eventdata, handles)
% hObject    handle to convert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename=handles.pov_filename;
dxf_filename=handles.dxf_filename;
matlab_surface_file=handles.mat_filename;
reduce_patch=handles.number_of_faces;
if get(handles.checkbox1,'Value')
    total_faces=tom_pov2dxf(filename,dxf_filename,reduce_patch,matlab_surface_file);
else
    total_faces=tom_pov2dxf(filename,dxf_filename,reduce_patch);
end;
msgbox([ num2str(total_faces) ' faces created'],'tom_pov2dxf','modal');

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close (gcf);

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the exit flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to exit the data.

% Update handles structure
guidata(handles.figure1, handles);


% --- Executes on button press in browse_pov_filename.
function browse_pov_filename_Callback(hObject, eventdata, handles)
% hObject    handle to browse_pov_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
[filename, pathname] = uigetfile('*.pov', 'Pick a POV-file');
    if isequal(filename,0) || isequal(pathname,0)
       disp('cancel')
    else
       handles.pov_filename=fullfile(pathname, filename);
    end
 
set(handles.pov_filename_display,'String',handles.pov_filename);
guidata(handles.figure1, handles);

% --- Executes on button press in browse_dxf_filename.
function browse_dxf_filename_Callback(hObject, eventdata, handles)
% hObject    handle to browse_dxf_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile('*.dxf', 'Pick a DXF-file');
    if isequal(filename,0) || isequal(pathname,0)
       disp('cancel')
    else
       handles.dxf_filename=fullfile(pathname, filename);
    end
 
set(handles.dxf_filename_display,'String',handles.dxf_filename);
guidata(handles.figure1, handles);



function dxf_filename_display_Callback(hObject, eventdata, handles)
% hObject    handle to dxf_filename_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dxf_filename_display as text
%        str2double(get(hObject,'String')) returns contents of dxf_filename_display as a double
handles.dxf_filename=get(hObject,'String');
w=what;
handles.dxf_filename=fullfile(w.path,handles.dxf_filename);
set(hObject,'String',handles.dxf_filename);
guidata(handles.figure1, handles);


% --- Executes during object creation, after setting all properties.
function dxf_filename_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dxf_filename_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.dxf_filename='';
set(hObject,'String',handles.dxf_filename);
handles.dxf_filename=get(hObject,'String');
guidata(hObject,handles);



function mat_filename_display_Callback(hObject, eventdata, handles)
% hObject    handle to mat_filename_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mat_filename_display as text
%        str2double(get(hObject,'String')) returns contents of mat_filename_display as a double
 
handles.mat_filename=get(hObject,'String');
w=what;
handles.mat_filename=fullfile(w.path,handles.mat_filename);
set(hObject,'String',handles.mat_filename);
guidata(handles.figure1, handles);


% --- Executes during object creation, after setting all properties.
function mat_filename_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mat_filename_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.mat_filename='';
set(hObject,'String',handles.mat_filename);
handles.mat_filename=get(hObject,'String');
guidata(hObject,handles);


% --- Executes on button press in browse_mat_filename.
function browse_mat_filename_Callback(hObject, eventdata, handles)
% hObject    handle to browse_mat_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile('*.mat', 'Pick a MAT-file');
    if isequal(filename,0) || isequal(pathname,0)
       disp('cancel')
    else
       handles.mat_filename=fullfile(pathname, filename);
    end
 
set(handles.mat_filename_display,'String',handles.mat_filename);
guidata(handles.figure1, handles);



function number_of_faces_display_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_faces_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_faces_display as text
%        str2double(get(hObject,'String')) returns contents of number_of_faces_display as a double
handles.number_of_faces=str2num(get(handles.number_of_faces_display,'String'));
handles.number_of_faces=handles.number_of_faces./100; % in per cent
guidata(handles.figure1, handles);

% --- Executes during object creation, after setting all properties.
function number_of_faces_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_faces_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
if get(hObject,'Value')
    set(handles.mat_filename_display,'Enable','on');
    set(handles.browse_mat_filename,'Enable','on');    
else
    set(handles.mat_filename_display,'Enable','off');
    set(handles.browse_mat_filename,'Enable','off');
end;
