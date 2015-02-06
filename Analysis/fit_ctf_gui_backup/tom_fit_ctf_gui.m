function varargout = tom_fit_ctf_gui(varargin)
% TOM_FIT_CTF_GUI M-file for tom_fit_ctf_gui.fig
%      TOM_FIT_CTF_GUI, by itself, creates a new TOM_FIT_CTF_GUI or raises the existing
%      singleton*.
%
%      H = TOM_FIT_CTF_GUI returns the handle to a new TOM_FIT_CTF_GUI or
%      the handle to
%      the existing singleton*.
%
%      TOM_FIT_CTF_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_FIT_CTF_GUI.M with the given input arguments.
%
%      TOM_FIT_CTF_GUI('Property','Value',...) creates a new
%      TOM_FIT_CTF_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before tom_fit_ctf_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_fit_ctf_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_fit_ctf_gui

% Last Modified by GUIDE v2.5 26-Feb-2010 11:18:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_fit_ctf_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_fit_ctf_gui_OutputFcn, ...
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


% --- Executes just before tom_fit_ctf_gui is made visible.
function tom_fit_ctf_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_fit_ctf_gui (see VARARGIN)

% Choose default command line output for tom_fit_ctf_gui
handles.output = hObject;

handles.size=[256 256];
handles.mask_inner_radius=9;
handles.mask_outer_radius=handles.size(1)./2-1;
handles.correctbackground_inner_radius=8;
handles.correctbackground_outer_radius=handles.size(1)./2;
handles.filter=0;
handles.correctbackground_method='calc_decay';
set(handles.gui_mask_inner_radius,'String',handles.mask_inner_radius);
set(handles.gui_mask_outer_radius,'String',handles.mask_outer_radius);
set(handles.gui_background_mask_inner_radius,'String',handles.correctbackground_inner_radius);
set(handles.gui_background_mask_outer_radius,'String',handles.correctbackground_outer_radius);
set(handles.gui_defocus_fixed_button,'Value',1);
set(handles.gui_defocus_delta_fixed_button,'Value',1);
set(handles.gui_Phi_0_fixed_button,'Value',1);
set(handles.gui_amplitude_contrast_button,'Value',1);
set(handles.gui_decay_button,'Value',1);
set(handles.gui_decay_part_coh_ill_button,'Value',1);
set(handles.gui_decay_energy_spread_button,'Value',1);
set(handles.show_minima_red_button,'Value',1);
set(handles.show_maxima_green_button,'Value',1);
set(handles.show_zeros_yellow_button,'Value',1);

set(handles.Defocus_Slider,'SliderStep',[.1./2000 .1./2000]);

handles.ctf_st=[];
handles.old_name=[];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_fit_ctf_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_fit_ctf_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function gui_filename_Callback(hObject, eventdata, handles)
% hObject    handle to gui_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_filename as text
%        str2double(get(hObject,'String')) returns contents of gui_filename as a double


% --- Executes during object creation, after setting all properties.
function gui_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in gui_Size.
function gui_Size_Callback(hObject, eventdata, handles)
% hObject    handle to gui_Size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns gui_Size contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gui_Size
contents = get(hObject,'String');
sz=contents{get(hObject,'Value')};
if ~isempty(findstr(sz,'User'))
    handles.size=str2num(sz(2:8));
    handles.region(1)=str2num(sz(12:15));
    handles.region(2)=handles.region(1);
else
    handles.size=str2num(sz);
    handles.region=handles.image_fullsize;
end;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function gui_Size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_Size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in gui_generator.
function gui_generator_Callback(hObject, eventdata, handles)
% hObject    handle to gui_generator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns gui_generator contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gui_generator


% --- Executes during object creation, after setting all properties.
function gui_generator_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_generator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {'*.em;', 'All TOM Files (*.em)'; ...
    '*.*',                   'All Files (*.*)'}, ...
    'Pick an image file');
 
if filename==0
    disp('image load cancelled.');
else
    handles=reload_it_man(handles,pathname,filename,0);
    set(handles.gui_defocus_search_start,'String',handles.intended_defocus-1);
    tmp=get(handles.gui_defocus_search_step,'String');
    
    if (strcmp(tmp,'Step')==1)
        set(handles.gui_defocus_search_step,'String',.1);
    end;
    
    set(handles.gui_defocus_search_stop,'String',handles.intended_defocus+1);
    set(handles.gui_defocus_fixed,'String',handles.defocus_fixed);
    %set(handles.gui_defocus_fixed_button,'Value',1);
    set(handles.Defocus_Slider,'Value',handles.defocus_fixed);
    set(handles.gui_intended_defocus,'String',handles.intended_defocus);
    
    tmp_fnames=dir([handles.pathname '/*' handles.file_ext ]);
    for i=1:length(tmp_fnames)
        new_filenames{i}=[handles.pathname '/' tmp_fnames(i).name];
    end;
    handles.folder_filenames=new_filenames;
    num_of_files=length(new_filenames);
    
    for i=1:length(handles.folder_filenames); 
        if isempty(strfind(handles.folder_filenames{i},filename))==0
            val=i;
            break;
        end; 
    end;
    set(handles.gui_file_slider,'Max',num_of_files);
    set(handles.gui_file_slider,'Min',1);
    set(handles.gui_file_slider,'Value',val);
    set(handles.gui_file_slider,'SliderStep',[(1./(num_of_files-1)) (1./(num_of_files-1))]);
end;
guidata(hObject, handles);


% --- Executes on button press in create.
function create_Callback(hObject, eventdata, handles)
% hObject    handle to create (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=create_ps(handles);
set(handles.gui_mask_outer_radius,'String',handles.mask_outer_radius);

guidata(hObject, handles);

function gui_magnification_Callback(hObject, eventdata, handles)
% hObject    handle to gui_magnification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_magnification as text
%        str2double(get(hObject,'String')) returns contents of gui_magnification as a double
handles.magnification=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function gui_magnification_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_magnification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_cc_Callback(hObject, eventdata, handles)
% hObject    handle to gui_cc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_cc as text
%        str2double(get(hObject,'String')) returns contents of gui_cc as a double
handles.cc=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function gui_cc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_cc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_pixelsize_Callback(hObject, eventdata, handles)
% hObject    handle to gui_pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_pixelsize as text
%        str2double(get(hObject,'String')) returns contents of gui_pixelsize as a double
handles.pixelsize=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function gui_pixelsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_voltage_Callback(hObject, eventdata, handles)
% hObject    handle to gui_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_voltage as text
%        str2double(get(hObject,'String')) returns contents of gui_voltage as a double
handles.voltage=str2double(get(hObject,'String')).*1000;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function gui_voltage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_cs_Callback(hObject, eventdata, handles)
% hObject    handle to gui_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_cs as text
%        str2double(get(hObject,'String')) returns contents of gui_cs as a double

handles.cs=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function gui_cs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in gui_filter.
function gui_filter_Callback(hObject, eventdata, handles)
% hObject    handle to gui_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns gui_filter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gui_filter
contents = get(hObject,'String');
if strcmp(contents{get(hObject,'Value')},'real space 3x3')==1
    handles.filter=3;
else
    handles.filter=0;
end;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function gui_filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gui_filter_create.
function gui_filter_create_Callback(hObject, eventdata, handles)
% hObject    handle to gui_filter_create (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filter~=0
    handles.ps=tom_filter(handles.ps,handles.filter,'quadr','real');
    axes(handles.main_axes);imagesc(handles.ps');axis image;colormap gray;
end;
guidata(hObject, handles);



function gui_mask_inner_radius_Callback(hObject, eventdata, handles)
% hObject    handle to gui_mask_inner_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_mask_inner_radius as text
%        str2double(get(hObject,'String')) returns contents of gui_mask_inner_radius as a double
handles.mask_inner_radius=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function gui_mask_inner_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_mask_inner_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_mask_outer_radius_Callback(hObject, eventdata, handles)
% hObject    handle to gui_mask_outer_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_mask_outer_radius as text
%        str2double(get(hObject,'String')) returns contents of gui_mask_outer_radius as a double
handles.mask_outer_radius=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function gui_mask_outer_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_mask_outer_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gui_apply_mask.
function gui_apply_mask_Callback(hObject, eventdata, handles)
% hObject    handle to gui_apply_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=apply_mask_bg(handles);
guidata(hObject, handles);


% --- Executes on selection change in gui_correctbackground.
function gui_correctbackground_Callback(hObject, eventdata, handles)
% hObject    handle to gui_correctbackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns gui_correctbackground contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gui_correctbackground
handles.correctbackground_method = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function gui_correctbackground_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_correctbackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gui_correctbackgound.
function gui_correctbackgound_Callback(hObject, eventdata, handles)
% hObject    handle to gui_correctbackgound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=correct_bg(handles);

guidata(hObject, handles);


% --- Executes on button press in gui_reload.
function gui_reload_Callback(hObject, eventdata, handles)
% hObject    handle to gui_reload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pathname=handles.pathname;
filename=handles.filename;

[header] = tom_reademheader([handles.pathname handles.filename]);
set(handles.gui_filename,'String',handles.filename);
%handles.magnification=header.Header.Magnification;
%handles.pixelsize=header.Header.Objectpixelsize;
%handles.intended_defocus=header.Header.Defocus./10000;
%handles.defocus_fixed=handles.intended_defocus;
%handles.cs=header.Header.Cs;
%handles.cc=2.2;
%handles.voltage=header.Header.Voltage;
in=tom_emreadc([pathname filename],'resample',[header.Header.Size(1)./512 header.Header.Size(2)./512 1]);
axes(handles.main_axes);imagesc(in.Value');axis image;colormap gray;
set(handles.gui_magnification,'String',handles.magnification);
set(handles.gui_pixelsize,'String',handles.pixelsize);
set(handles.gui_voltage,'String',handles.voltage./1000);
set(handles.gui_cs,'String',handles.cs);
set(handles.gui_cc,'String',handles.cc);
set(handles.gui_intended_defocus,'String',handles.intended_defocus);
set(handles.gui_defocus_fixed,'String',handles.defocus_fixed);
disp('reloaded');
guidata(hObject, handles);



function gui_background_mask_inner_radius_Callback(hObject, eventdata, handles)
% hObject    handle to gui_background_mask_inner_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_background_mask_inner_radius as text
%        str2double(get(hObject,'String')) returns contents of gui_background_mask_inner_radius as a double
handles.correctbackground_inner_radius=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function gui_background_mask_inner_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_background_mask_inner_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_background_mask_outer_radius_Callback(hObject, eventdata, handles)
% hObject    handle to gui_background_mask_outer_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_background_mask_outer_radius as text
%        str2double(get(hObject,'String')) returns contents of gui_background_mask_outer_radius as a double

handles.correctbackground_outer_radius=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function gui_background_mask_outer_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_background_mask_outer_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function gui_intended_defocus_Callback(hObject, eventdata, handles)
% hObject    handle to gui_intended_defocus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_intended_defocus as text
%        str2double(get(hObject,'String')) returns contents of gui_intended_defocus as a double
handles.intended_defocus=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function gui_intended_defocus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_intended_defocus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gui_defocus_fixed_button.
function gui_defocus_fixed_button_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_fixed_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_defocus_fixed_button


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in gui_defocus_search_button.
function gui_defocus_search_button_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_search_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_defocus_search_button



function gui_defocus_fixed_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_defocus_fixed as text
%        str2double(get(hObject,'String')) returns contents of gui_defocus_fixed as a double
defocus=str2double(get(hObject,'String'));
set(handles.Defocus_Slider,'Value',defocus);

% --- Executes during object creation, after setting all properties.
function gui_defocus_fixed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_defocus_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
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



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_do_fit_Callback(hObject, eventdata, handles)
% hObject    handle to gui_do_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_do_fit as text
%        str2double(get(hObject,'String')) returns contents of gui_do_fit as a double




% --- Executes during object creation, after setting all properties.
function gui_do_fit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_do_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function gui_defocus_search_start_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_defocus_search_start as text
%        str2double(get(hObject,'String')) returns contents of gui_defocus_search_start as a double


% --- Executes during object creation, after setting all properties.
function gui_defocus_search_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_defocus_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_defocus_search_step_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_defocus_search_step as text
%        str2double(get(hObject,'String')) returns contents of gui_defocus_search_step as a double
step=str2double(get(hObject,'String'));
set(handles.Defocus_Slider,'SliderStep',[step./2000 1]);

% --- Executes during object creation, after setting all properties.
function gui_defocus_search_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_defocus_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_defocus_search_stop_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_defocus_search_stop as text
%        str2double(get(hObject,'String')) returns contents of gui_defocus_search_stop as a double


% --- Executes during object creation, after setting all properties.
function gui_defocus_search_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_defocus_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gui_defocus_center.
function gui_defocus_center_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mid=str2num(get(handles.gui_defocus_fixed,'String'));
step=str2num(get(handles.gui_defocus_search_step,'String'));
mid_diff=mid./10;
set(handles.gui_defocus_search_start,'String',sprintf('%0.5g',round((mid+mid_diff.*2)./step).*step));
set(handles.gui_defocus_search_stop,'String',sprintf('%0.5g',round((mid-mid_diff.*2)./step).*step));


function [image EM Search]=get_gui_values(hObject, eventdata, handles)

% PS
image=handles.ps;

% load the EM structure with the parameters from the GUI:
EM.Objectpixelsize=str2double(get(handles.gui_pixelsize,'String')).*1e-10;
EM.Voltage=str2double(get(handles.gui_voltage,'String')).*1000.0;
EM.Cs=str2double(get(handles.gui_cs,'String')).*1e-3;
EM.Cc=str2double(get(handles.gui_cc,'String')).*1e-3;

% load the Search structure with the parameters from the GUI:
% Defocus
if get(handles.gui_defocus_fixed_button,'Value')
    Search.Dz_search=str2double(get(handles.gui_defocus_fixed,'String'));
else
    Search.Dz_search=[  str2double(get(handles.gui_defocus_search_start,'String')) : ...
                        str2double(get(handles.gui_defocus_search_step,'String')) : ...
                        str2double(get(handles.gui_defocus_search_stop,'String'))];
end;
Search.Dz_search=Search.Dz_search.*1e-6; % in m !!!
% Defocus delta, Astigmatism
if get(handles.gui_defocus_delta_fixed_button,'Value')
    Search.Dz_delta_search=str2double(get(handles.gui_defocus_delta_fixed,'String'));
else
    Search.Dz_delta_search=[  str2double(get(handles.gui_defocus_delta_search_start,'String')) : ...
                        str2double(get(handles.gui_defocus_delta_search_step,'String')) : ...
                        str2double(get(handles.gui_defocus_delta_search_stop,'String'))];
end;
Search.Dz_delta_search=Search.Dz_delta_search.*1e-6; % in m !!!


% Defocus delta angle, Astigmatism angle
if get(handles.gui_Phi_0_fixed_button,'Value')
    Search.Phi_0_search=str2double(get(handles.gui_Phi_0_fixed,'String'));
else
    Search.Phi_0_search=[  str2double(get(handles.gui_Phi_0_search_start,'String')) : ...
                        str2double(get(handles.gui_Phi_0_search_step,'String')) : ...
                        str2double(get(handles.gui_Phi_0_search_stop,'String'))];
end;

% Amplitude contrast, Ice is 0.07, negative stain is 0.15.
if get(handles.gui_amplitude_contrast_button,'Value')
    Search.amplitude_contrast_search=str2double(get(handles.gui_amplitude_contrast_fixed,'String'));
else
    Search.amplitude_contrast_search=[  str2double(get(handles.gui_ampltude_contrast_search_start,'String')) : ...
                        str2double(get(handles.gui_ampltude_contrast_search_step,'String')) : ...
                        str2double(get(handles.gui_ampltude_contrast_search_stop,'String'))];
end;

% Decay k^2
if get(handles.gui_decay_button,'Value')
    Search.decay_search=str2double(get(handles.gui_decay_fixed,'String'));
else
    Search.decay_search=[  str2double(get(handles.gui_decay_search_start,'String')) : ...
                        str2double(get(handles.gui_decay_search_step,'String')) : ...
                        str2double(get(handles.gui_decay_search_stop,'String'))];
end;

% Decay partially coherent illumination
if get(handles.gui_decay_part_coh_ill_button,'Value')
    Search.decay_part_coh_ill_search=str2double(get(handles.gui_decay_part_coh_ill_fixed,'String'));
else
    Search.decay_part_coh_ill_search=[  str2double(get(handles.gui_decay_part_coh_ill_search_start,'String')) : ...
                        str2double(get(handles.gui_decay_part_coh_ill_search_step,'String')) : ...
                        str2double(get(handles.gui_decay_part_coh_ill_search_stop,'String'))];
end;

% Decay due to energy spread
if get(handles.gui_decay_energy_spread_button,'Value')
    Search.decay_energy_spread_search=str2double(get(handles.gui_energy_spread_fixed,'String'));
else
    Search.decay_energy_spread_search=[  str2double(get(handles.gui_decay_energy_spread_search_start,'String')) : ...
                        str2double(get(handles.gui_decay_energy_spread_search_step,'String')) : ...
                        str2double(get(handles.gui_decay_energy_spread_search_stop,'String'))];
end;

Search.mask_inner_radius=str2double(get(handles.gui_mask_inner_radius,'String'));
Search.mask_outer_radius=str2double(get(handles.gui_mask_outer_radius,'String'));

Search.mask_inner_radius_bg=str2double(get(handles.gui_background_mask_inner_radius,'String'));
Search.mask_outer_radius_bg=str2double(get(handles.gui_background_mask_outer_radius,'String'));

str=get(handles.gui_filter,'String');
tmp_filt=str(get(handles.gui_filter,'Value'));

if (strcmp(tmp_filt,'real space 3x3'))
    Search.ps_filter=3;
end;

str=get(handles.gui_Size,'String');
tmp_filt=str(get(handles.gui_Size,'Value'));

Search.ps_size=str2num(tmp_filt{1});



guidata(hObject, handles);






% --- Executes on button press in gui_defocus_delta_fixed_button.
function gui_defocus_delta_fixed_button_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_delta_fixed_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_defocus_delta_fixed_button


% --- Executes on button press in gui_defocus_delta_search_button.
function gui_defocus_delta_search_button_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_delta_search_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_defocus_delta_search_button



function gui_defocus_delta_fixed_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_delta_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_defocus_delta_fixed as text
%        str2double(get(hObject,'String')) returns contents of gui_defocus_delta_fixed as a double


% --- Executes during object creation, after setting all properties.
function gui_defocus_delta_fixed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_defocus_delta_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_defocus_delta_search_start_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_delta_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_defocus_delta_search_start as text
%        str2double(get(hObject,'String')) returns contents of gui_defocus_delta_search_start as a double


% --- Executes during object creation, after setting all properties.
function gui_defocus_delta_search_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_defocus_delta_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_defocus_delta_search_step_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_delta_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_defocus_delta_search_step as text
%        str2double(get(hObject,'String')) returns contents of gui_defocus_delta_search_step as a double


% --- Executes during object creation, after setting all properties.
function gui_defocus_delta_search_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_defocus_delta_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_defocus_delta_search_stop_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_delta_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_defocus_delta_search_stop as text
%        str2double(get(hObject,'String')) returns contents of gui_defocus_delta_search_stop as a double


% --- Executes during object creation, after setting all properties.
function gui_defocus_delta_search_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_defocus_delta_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gui_defocus_delta_center.
function gui_defocus_delta_center_Callback(hObject, eventdata, handles)
% hObject    handle to gui_defocus_delta_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mid=str2num(get(handles.gui_defocus_delta_fixed,'String'));
step=str2num(get(handles.gui_defocus_delta_search_step,'String'));
mid_diff=mid./10;
set(handles.gui_defocus_delta_search_start,'String',sprintf('%0.5g',round((mid-mid_diff.*2)./step).*step));
set(handles.gui_defocus_delta_search_stop,'String',sprintf('%0.5g',round((mid+mid_diff.*2)./step).*step));



% --- Executes on button press in gui_Phi_0_fixed_button.
function gui_Phi_0_fixed_button_Callback(hObject, eventdata, handles)
% hObject    handle to gui_Phi_0_fixed_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_Phi_0_fixed_button


% --- Executes on button press in gui_Phi_0_search_button.
function gui_Phi_0_search_button_Callback(hObject, eventdata, handles)
% hObject    handle to gui_Phi_0_search_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_Phi_0_search_button



function gui_Phi_0_fixed_Callback(hObject, eventdata, handles)
% hObject    handle to gui_Phi_0_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_Phi_0_fixed as text
%        str2double(get(hObject,'String')) returns contents of gui_Phi_0_fixed as a double


% --- Executes during object creation, after setting all properties.
function gui_Phi_0_fixed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_Phi_0_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_Phi_0_search_start_Callback(hObject, eventdata, handles)
% hObject    handle to gui_Phi_0_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_Phi_0_search_start as text
%        str2double(get(hObject,'String')) returns contents of gui_Phi_0_search_start as a double


% --- Executes during object creation, after setting all properties.
function gui_Phi_0_search_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_Phi_0_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_Phi_0_search_step_Callback(hObject, eventdata, handles)
% hObject    handle to gui_Phi_0_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_Phi_0_search_step as text
%        str2double(get(hObject,'String')) returns contents of gui_Phi_0_search_step as a double


% --- Executes during object creation, after setting all properties.
function gui_Phi_0_search_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_Phi_0_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_Phi_0_search_stop_Callback(hObject, eventdata, handles)
% hObject    handle to gui_Phi_0_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_Phi_0_search_stop as text
%        str2double(get(hObject,'String')) returns contents of gui_Phi_0_search_stop as a double


% --- Executes during object creation, after setting all properties.
function gui_Phi_0_search_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_Phi_0_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gui_Phi_0_center.
function gui_Phi_0_center_Callback(hObject, eventdata, handles)
% hObject    handle to gui_Phi_0_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mid=str2num(get(handles.gui_Phi_0_fixed,'String'));
step=str2num(get(handles.gui_Phi_0_search_step,'String'));
mid_diff=mid./10;
set(handles.gui_Phi_0_search_start,'String',sprintf('%0.5g',round((mid-mid_diff.*2)./step).*step));
set(handles.gui_Phi_0_search_stop,'String',sprintf('%0.5g',round((mid+mid_diff.*2)./step).*step));


% --- Executes on button press in gui_amplitude_contrast_button.
function gui_amplitude_contrast_button_Callback(hObject, eventdata, handles)
% hObject    handle to gui_amplitude_contrast_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_amplitude_contrast_button



function gui_amplitude_contrast_fixed_Callback(hObject, eventdata, handles)
% hObject    handle to gui_amplitude_contrast_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_amplitude_contrast_fixed as text
%        str2double(get(hObject,'String')) returns contents of gui_amplitude_contrast_fixed as a double


% --- Executes during object creation, after setting all properties.
function gui_amplitude_contrast_fixed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_amplitude_contrast_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_ampltude_contrast_search_start_Callback(hObject, eventdata, handles)
% hObject    handle to gui_ampltude_contrast_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_ampltude_contrast_search_start as text
%        str2double(get(hObject,'String')) returns contents of gui_ampltude_contrast_search_start as a double


% --- Executes during object creation, after setting all properties.
function gui_ampltude_contrast_search_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_ampltude_contrast_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_ampltude_contrast_search_step_Callback(hObject, eventdata, handles)
% hObject    handle to gui_ampltude_contrast_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_ampltude_contrast_search_step as text
%        str2double(get(hObject,'String')) returns contents of gui_ampltude_contrast_search_step as a double


% --- Executes during object creation, after setting all properties.
function gui_ampltude_contrast_search_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_ampltude_contrast_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_ampltude_contrast_search_stop_Callback(hObject, eventdata, handles)
% hObject    handle to gui_ampltude_contrast_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_ampltude_contrast_search_stop as text
%        str2double(get(hObject,'String')) returns contents of gui_ampltude_contrast_search_stop as a double


% --- Executes during object creation, after setting all properties.
function gui_ampltude_contrast_search_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_ampltude_contrast_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gui_ampltude_contrast_center.
function gui_ampltude_contrast_center_Callback(hObject, eventdata, handles)
% hObject    handle to gui_ampltude_contrast_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mid=str2num(get(handles.gui_amplitude_contrast_fixed,'String'));
step=str2num(get(handles.gui_ampltude_contrast_search_step,'String'));
mid_diff=mid./10;
set(handles.gui_ampltude_contrast_search_start,'String',sprintf('%0.5g',round((mid-mid_diff.*2)./step).*step));
set(handles.gui_ampltude_contrast_search_stop,'String',sprintf('%0.5g',round((mid+mid_diff.*2)./step).*step));


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in gui_decay_button.
function gui_decay_button_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_decay_button



function gui_decay_fixed_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_decay_fixed as text
%        str2double(get(hObject,'String')) returns contents of gui_decay_fixed as a double


% --- Executes during object creation, after setting all properties.
function gui_decay_fixed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_decay_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_decay_search_start_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_decay_search_start as text
%        str2double(get(hObject,'String')) returns contents of gui_decay_search_start as a double


% --- Executes during object creation, after setting all properties.
function gui_decay_search_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_decay_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_decay_search_step_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_decay_search_step as text
%        str2double(get(hObject,'String')) returns contents of gui_decay_search_step as a double


% --- Executes during object creation, after setting all properties.
function gui_decay_search_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_decay_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_decay_search_stop_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_decay_search_stop as text
%        str2double(get(hObject,'String')) returns contents of gui_decay_search_stop as a double


% --- Executes during object creation, after setting all properties.
function gui_decay_search_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_decay_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gui_decay_center.
function gui_decay_center_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mid=str2num(get(handles.gui_decay_fixed,'String'));
step=str2num(get(handles.gui_decay_search_step,'String'));
mid_diff=mid./10;
set(handles.gui_decay_search_start,'String',sprintf('%0.5g',round((mid-mid_diff.*2)./step).*step));
set(handles.gui_decay_search_stop,'String',sprintf('%0.5g',round((mid+mid_diff.*2)./step).*step));


% --- Executes on button press in gui_decay_part_coh_ill_button.
function gui_decay_part_coh_ill_button_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_part_coh_ill_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_decay_part_coh_ill_button



function gui_decay_part_coh_ill_fixed_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_part_coh_ill_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_decay_part_coh_ill_fixed as text
%        str2double(get(hObject,'String')) returns contents of gui_decay_part_coh_ill_fixed as a double


% --- Executes during object creation, after setting all properties.
function gui_decay_part_coh_ill_fixed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_decay_part_coh_ill_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_decay_part_coh_ill_search_start_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_part_coh_ill_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_decay_part_coh_ill_search_start as text
%        str2double(get(hObject,'String')) returns contents of gui_decay_part_coh_ill_search_start as a double


% --- Executes during object creation, after setting all properties.
function gui_decay_part_coh_ill_search_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_decay_part_coh_ill_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_decay_part_coh_ill_search_step_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_part_coh_ill_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_decay_part_coh_ill_search_step as text
%        str2double(get(hObject,'String')) returns contents of gui_decay_part_coh_ill_search_step as a double


% --- Executes during object creation, after setting all properties.
function gui_decay_part_coh_ill_search_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_decay_part_coh_ill_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_decay_part_coh_ill_search_stop_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_part_coh_ill_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_decay_part_coh_ill_search_stop as text
%        str2double(get(hObject,'String')) returns contents of gui_decay_part_coh_ill_search_stop as a double


% --- Executes during object creation, after setting all properties.
function gui_decay_part_coh_ill_search_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_decay_part_coh_ill_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gui_decay_part_coh_ill_center.
function gui_decay_part_coh_ill_center_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_part_coh_ill_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mid=str2num(get(handles.gui_decay_part_coh_ill_fixed,'String'));
step=str2num(get(handles.gui_decay_part_coh_ill_search_step,'String'));
mid_diff=mid./10;
set(handles.gui_decay_part_coh_ill_search_start,'String',sprintf('%0.5g',round((mid-mid_diff.*2)./step).*step));
set(handles.gui_decay_part_coh_ill_search_stop,'String',sprintf('%0.5g',round((mid+mid_diff.*2)./step).*step));


% --- Executes on button press in copy_Fit_to_workspace.
function copy_Fit_to_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to copy_Fit_to_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[image EM Search]=get_gui_values(hObject, eventdata, handles);


ctf.Fit=handles.Fit;
ctf.Search=Search;
assignin('base','ctf_st',ctf);
assignin('base','Fit',handles.Fit);
assignin('base','EM',EM);
assignin('base','Search',Search);
assignin('base','ps',handles.ps);
assignin('base','ps_orig',handles.untreated_ps);


% --- Executes on button press in gui_decay_energy_spread_button.
function gui_decay_energy_spread_button_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_energy_spread_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_decay_energy_spread_button



function gui_energy_spread_fixed_Callback(hObject, eventdata, handles)
% hObject    handle to gui_energy_spread_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_energy_spread_fixed as text
%        str2double(get(hObject,'String')) returns contents of gui_energy_spread_fixed as a double


% --- Executes during object creation, after setting all properties.
function gui_energy_spread_fixed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_energy_spread_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_decay_energy_spread_search_start_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_energy_spread_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_decay_energy_spread_search_start as text
%        str2double(get(hObject,'String')) returns contents of gui_decay_energy_spread_search_start as a double


% --- Executes during object creation, after setting all properties.
function gui_decay_energy_spread_search_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_decay_energy_spread_search_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_decay_energy_spread_search_step_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_energy_spread_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_decay_energy_spread_search_step as text
%        str2double(get(hObject,'String')) returns contents of gui_decay_energy_spread_search_step as a double


% --- Executes during object creation, after setting all properties.
function gui_decay_energy_spread_search_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_decay_energy_spread_search_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_decay_energy_spread_search_stop_Callback(hObject, eventdata, handles)
% hObject    handle to gui_decay_energy_spread_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_decay_energy_spread_search_stop as text
%        str2double(get(hObject,'String')) returns contents of gui_decay_energy_spread_search_stop as a double


% --- Executes during object creation, after setting all properties.
function gui_decay_energy_spread_search_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_decay_energy_spread_search_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gui_energy_spread_center.
function gui_energy_spread_center_Callback(hObject, eventdata, handles)
% hObject    handle to gui_energy_spread_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mid=str2num(get(handles.gui_energy_spread_fixed,'String'));
step=str2num(get(handles.gui_decay_energy_spread_search_step,'String'));
mid_diff=mid./10;
set(handles.gui_decay_energy_spread_search_start,'String',sprintf('%0.5g',round((mid-mid_diff.*2)./step).*step));
set(handles.gui_decay_energy_spread_search_stop,'String',sprintf('%0.5g',round((mid+mid_diff.*2)./step).*step));

% --- Executes on button press in adapt_contrast_button.
function adapt_contrast_button_Callback(hObject, eventdata, handles)
% hObject    handle to adapt_contrast_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
imcontrast(handles.main_axes);
warning on;

function [im_handle h_amplitude h_phase defocus_label_points_x defocus_label_points_y fit_ps]=show_ctf_fit(in, EM, Fit)

[im_handle h_amplitude h_phase defocus_label_points_x defocus_label_points_y fit_ps]=tom_show_ctf_fit(in, EM, Fit);



% --- Executes on button press in gui_perform_fit.
function gui_perform_fit_Callback(hObject, eventdata, handles)
% hObject    handle to gui_perform_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=fit_it(handles,eventdata,hObject);

guidata(hObject, handles);


% --- Executes on button press in show_minima_red_button.
function show_minima_red_button_Callback(hObject, eventdata, handles)
% hObject    handle to show_minima_red_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of show_minima_red_button

if get(hObject,'Value')==1
    set(handles.h_phase,'Visible','on');
else
    set(handles.h_phase,'Visible','off');
end;

% --- Executes on button press in show_maxima_green_button.
function show_maxima_green_button_Callback(hObject, eventdata, handles)
% hObject    handle to show_maxima_green_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of show_maxima_green_button

if get(hObject,'Value')==1
    set(handles.h_amplitude,'Visible','on');
else
    set(handles.h_amplitude,'Visible','off');
end;


% --- Executes on button press in show_zeros_yellow_button.
function show_zeros_yellow_button_Callback(hObject, eventdata, handles)
% hObject    handle to show_zeros_yellow_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_zeros_yellow_button
if get(hObject,'Value')==1
    set(handles.label_x,'Visible','on');
    set(handles.label_y,'Visible','on');
    set(handles.label_x_dots,'Visible','on');
    set(handles.label_y_dots,'Visible','on');
else
    set(handles.label_x,'Visible','off');
    set(handles.label_x_dots,'Visible','off');
    set(handles.label_y_dots,'Visible','off');
    set(handles.label_y,'Visible','off');
end;


% --- Executes on slider movement.
function Defocus_Slider_Callback(hObject, eventdata, handles)
% hObject    handle to Defocus_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider



slider_defocus=get(hObject,'Value');

%defocus=str2double(get(handles.gui_defocus_fixed,'String'));
%set(hObject,'Value',defocus);
%increment=str2double(get(handles.gui_defocus_search_step,'String'));

set(handles.gui_defocus_fixed,'String',num2str(slider_defocus));
[image EM Search]=get_gui_values(hObject, eventdata, handles);

Fit=handles.Fit;
Fit.Dz_det=slider_defocus.*1e-06;
axes(handles.main_axes);
[im_handle h_amplitude h_phase defocus_label_points_x defocus_label_points_y fit_ps]=show_ctf_fit(image, EM, Fit);


handles.h_amplitude=h_amplitude;
handles.h_phase=h_phase;

%gui_perform_fit_Callback(hObject, eventdata, handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Defocus_Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Defocus_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in gui_robust_fit_button.
function gui_robust_fit_button_Callback(hObject, eventdata, handles)
% hObject    handle to gui_robust_fit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

txt=get(handles.gui_robust_fit_button,'String');
set(handles.gui_robust_fit_button,'String','Searching ...');drawnow;
%set(handles.gui_perform_fit,'Enable','inactive')
[image EM Search]=get_gui_values(hObject, eventdata, handles);

Dz_search_increment=2e-6;
Dz_delta_search_increment=2e-6;
Phi_0_search_increment=30;

Search.Dz_search=[-40e-6:Dz_search_increment:0]; 
Search.Dz_delta_search=[0:Dz_delta_search_increment:40e-6];
Search.Phi_0_search=[0:Phi_0_search_increment:150];
tic;
for i=1:5
    [Fit]=tom_fit_ctf(image,EM,Search);
    Search.Dz_search=[Fit.Dz_det-Dz_search_increment:Dz_search_increment./2:Fit.Dz_det+Dz_search_increment];
    if Fit.Dz_delta_det-Dz_delta_search_increment>=0
        Search.Dz_delta_search=[Fit.Dz_delta_det-Dz_delta_search_increment:Dz_delta_search_increment./2:Fit.Dz_delta_det+Dz_delta_search_increment];
    else
        Search.Dz_delta_search=[0:Dz_delta_search_increment./2:Fit.Dz_delta_det+Dz_delta_search_increment];
    end;        
%    Search.Phi_0_search=[Fit.Phi_0_det-Phi_0_search_increment:Phi_0_search_increment./2:Fit.Phi_0_det+Phi_0_search_increment];
    Search.Phi_0_search=[0:Phi_0_search_increment:150];
    Dz_search_increment=Dz_search_increment./2;
    Dz_delta_search_increment=Dz_delta_search_increment./2;
%    Phi_0_search_increment=Phi_0_search_increment./2;
    if Dz_search_increment<.1e-6 break; end;
    
end;
toc    
[im_handle h_amplitude h_phase defocus_label_points_x defocus_label_points_y]=show_ctf_fit(image, EM, Fit);

handles.EM=EM;
handles.Search=Search;
handles.Fit=Fit;


hold on;

defocus_vals_x=sqrt((defocus_label_points_x(:,1)-(size(handles.ps,1)./2)).^2+(defocus_label_points_x(:,2)-(size(handles.ps,1)./2)).^2)./(size(handles.ps,1)./2)./(handles.pixelsize.*2);
defocus_vals_y=sqrt((defocus_label_points_y(:,1)-(size(handles.ps,1)./2)).^2+(defocus_label_points_y(:,2)-(size(handles.ps,1)./2)).^2)./(size(handles.ps,1)./2)./(handles.pixelsize.*2);
handles.label_x=0;handles.label_y=0;
for i=1:size(defocus_vals_x,1)
    lbl_x=sprintf('%0.3g',1./defocus_vals_x(i));
    lbl_x=[' ' lbl_x];
    th=text(defocus_label_points_x(i,1),defocus_label_points_x(i,2),lbl_x);
    set(th,'Color','white');
    set(th,'Fontsize',16);
    handles.label_x(i)=th;
end;
for i=1:size(defocus_vals_y,1)
    lbl_y=sprintf('%0.3g',1./defocus_vals_y(i));
    lbl_y=[' ' lbl_y];
    th2=text(defocus_label_points_y(i,1),defocus_label_points_y(i,2),lbl_y);
    set(th2,'Color','white');
    set(th2,'Fontsize',16);
    handles.label_y(i)=th2;
end;

handles.label_x_dots=plot(defocus_label_points_x(:,1),defocus_label_points_x(:,2),'yo','MarkerSize',3,'Linewidth',5);
handles.label_y_dots=plot(defocus_label_points_y(:,1),defocus_label_points_y(:,2),'yo','MarkerSize',3,'Linewidth',5);

hold off;
handles.h_amplitude=h_amplitude;
handles.h_phase=h_phase;

set(handles.gui_defocus_fixed,'String',sprintf('%0.5g',Fit.Dz_det.*1e6));
set(handles.Defocus_Slider,'Value',Fit.Dz_det.*1e6);
set(handles.gui_text_long_axis,'String',['Defocus long axis: ' sprintf('%0.5g',Fit.Dz_det.*1e6+(Fit.Dz_delta_det.*1e6)./2) ' mue.']);
set(handles.gui_text_short_axis,'String',['Defocus short axis: ' sprintf('%0.5g',Fit.Dz_det.*1e6-(Fit.Dz_delta_det.*1e6)./2) ' mue.']);
set(handles.gui_text_angle_axis,'String',['Angle (long&short axis): ' sprintf('%0.5g',Fit.Phi_0_det) ' deg.']);

set(handles.gui_defocus_delta_fixed,'String',sprintf('%0.5g',Fit.Dz_delta_det.*1e6));
set(handles.gui_Phi_0_fixed,'String',sprintf('%0.5g',Fit.Phi_0_det));
set(handles.gui_amplitude_contrast_fixed,'String',sprintf('%0.3g',Fit.amplitude_contrast_det));
set(handles.gui_decay_fixed,'String',sprintf('%0.3g',Fit.decay_det));
set(handles.gui_decay_part_coh_ill_fixed,'String',sprintf('%0.3g',Fit.decay_part_coh_ill_det));
set(handles.gui_energy_spread_fixed,'String',sprintf('%0.3g',Fit.decay_energy_spread_det));

set(handles.gui_robust_fit_button,'String','Robust CTF Search');
set(handles.gui_robust_fit_button,'Enable','on');
drawnow;
set(handles.show_minima_red_button,'Value',1);
set(handles.show_maxima_green_button,'Value',1);
set(handles.show_zeros_yellow_button,'Value',1);
axes(handles.ccc_plot_axes); 
r=reshape(Fit.corr_all,[numel(Fit.corr_all) 1]);
[n,xout] = hist(r,15); h=bar(xout,n,'red');
xlim([-1 1]);
ylim([0 max(n)]);
%set(handles.ccc_plot_axes,'XTickLabel',[-1:0.5:1]);
set(handles.ccc_plot_axes,'XGrid','on');
set(handles.ccc_plot_axes,'Color',[0 0 0])
set(handles.ccc_plot_axes,'YColor',[1 1 1]); 
set(handles.ccc_plot_axes,'XColor',[1 1 1]); 
hold on;
[me ma mi st va]=tom_dev(Fit.corr_all,'noinfo');
ht=text(max(xout),n(end),num2str(sprintf('%0.2g',max(r))));
if ma>me+st
    set(ht,'FontSize',11,'Color','green','FontWeight','Bold');
end;
if (ma<me+st./2 || ma<0.5)
        set(ht,'FontSize',11,'Color','yellow','FontWeight','Bold');
end;
if ma < 0.3
            set(ht,'FontSize',11,'Color','red','FontWeight','Bold');
            set(handles.ccc_plot_axes,'Color',[1 1 0])
end;

hold off;

guidata(hObject, handles);


% --- Executes on slider movement.
function gui_file_slider_Callback(hObject, eventdata, handles)
% hObject    handle to gui_file_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val=round(get(handles.gui_file_slider,'Value'));

[a b c]=fileparts(handles.folder_filenames{val});
set(handles.gui_filename,'String',[b c]);

 br_file=get(handles.fl_sl_browse,'Value');
 br_fit=get(handles.fl_sl_fit,'Value');
 br_mon=get(handles.fl_sl_monitor,'Value');

 
 if (br_file)
    handles=reload_it_man(handles,handles.pathname,[b c],1);  drawnow;
 end;
     
 if (br_mon)
     set(handles.mon_checked,'Value',1);
     if (exist([handles.folder_filenames{val} '.mat'],'file'))
       if (isempty(handles.old_name)==0)
             [image EM Search]=get_gui_values(hObject, eventdata, handles);
             st_out.img=handles.ps;
             st_out.Fit=handles.Fit;
             st_out.Fit.Dz_det=str2num(get(handles.gui_defocus_fixed,'String')).*1e-06;
             st_out.Fit.EM=EM;
             st_out.Search=Search;
             st_out.sel.selected=get(handles.mon_accepted,'Value');
             st_out.sel.accepted=get(handles.mon_checked,'Value');
             st_out.sel.checked=get(handles.mon_selected,'Value');
             save(handles.old_name,'st_out');
         end;
         load([handles.folder_filenames{val} '.mat']);
         handles.old_name=[handles.folder_filenames{val} '.mat'];
         try
            set(handles.mon_accepted,'Value',st_out.sel.selected);
            set(handles.mon_checked,'Value',st_out.sel.accepted);
            set(handles.mon_selected,'Value',1);
         catch
            set(handles.mon_accepted,'Value',1);
            set(handles.mon_checked,'Value',0);
            set(handles.mon_selected,'Value',1);
         end;
         
         handles.ctf_st=st_out;
         try
            handles=set_gui_values(handles,st_out.Search,st_out.Fit.EM,st_out.img,st_out.Fit);
           %chang me fb!
            handles.ps=st_out.img;
            handles.background_corrected_ps=st_out.img;
         catch
            
         end;
      end;
 end;

if (br_fit)
    handles=reload_it_man(handles,handles.pathname,[b c],1);  drawnow;
    handles=create_ps(handles);drawnow;
    handles=correct_bg(handles); drawnow;
    handles=apply_mask_bg(handles);  drawnow;
    handles=fit_it(handles,eventdata,hObject);  drawnow;
end;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function gui_file_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_file_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider contreload_it_manrols usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in gui_fit_folder.
function gui_fit_folder_Callback(hObject, eventdata, handles)
% hObject    handle to gui_fit_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[image EM Search]=get_gui_values(hObject, eventdata, handles);


ctf.Fit=handles.Fit;
ctf.Search=Search;

tom_fit_ctf_folder(handles.pathname,handles.file_ext,handles.pathname,ctf,'none');



% --- Executes on button press in gui_fit_and_correct_folder.
function gui_fit_and_correct_folder_Callback(hObject, eventdata, handles)
% hObject    handle to gui_fit_and_correct_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[image EM Search]=get_gui_values(hObject, eventdata, handles);

h=tom_reademheader([handles.pathname handles.filename]);
im_sz=h.Header.Size(1);
pixs=EM.Objectpixelsize;
ny=2e10.*pixs;

hstr=get(handles.pop_corr_meth,'String');

method=hstr{get(handles.pop_corr_meth,'Value')};
if (strcmp(method,'flip&default mtf'))
    method='flip&mtf';
end;

corr_cut_off=str2num(get(handles.gui_cut_off_edit,'String'));
corr_cut_off_pix=round((im_sz./2).*(ny./corr_cut_off));

ctf.Fit=handles.Fit;
ctf.Search=Search;

if isempty(strfind(handles.pathname,'/high/') )
    target_path=strrep(handles.pathname,'/low/','/low_corr/');
else
    target_path=strrep(handles.pathname,'/high/','/high_corr/');
end;

if (strcmp(handles.pathname,target_path))
    Answ=questdlg('source path == target path ');
    if (strcmp(Answ,'Yes')==0)
        return;
    end;
end;

tom_fit_ctf_folder(handles.pathname,handles.file_ext,target_path,ctf,method,corr_cut_off_pix);




% --- Executes on button press in gui_corr_folder.
function gui_corr_folder_Callback(hObject, eventdata, handles)
% hObject    handle to gui_corr_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[image EM Search]=get_gui_values(hObject, eventdata, handles);

h=tom_reademheader([handles.pathname handles.filename]);
im_sz=h.Header.Size(1);
pixs=EM.Objectpixelsize;
ny=2e10.*pixs;

hstr=get(handles.pop_corr_meth,'String');

method=hstr{get(handles.pop_corr_meth,'Value')};
if (strcmp(method,'flip&default mtf'))
    method='flip&mtf';
end;

corr_cut_off=str2num(get(handles.gui_cut_off_edit,'String'));
corr_cut_off_pix=round((im_sz./2).*(ny./corr_cut_off));

ctf.Fit=handles.Fit;
ctf.Search=Search;

[a b c]=fileparts(handles.pathname);


% if isempty(strfind(handles.pathname,'/high/') )
%     target_path=strrep(handles.pathname,'/low/','/low_corr/');
% else
%     target_path=strrep(handles.pathname,'/high/','/high_corr/');
% end;

target_path=[a '_corr'];

if (strcmp(handles.pathname,target_path))
    Answ=questdlg('source path == target path ');
    if (strcmp(Answ,'Yes')==0)
        return;
    end;
end;

ctf.no_fit=1;

disp(['Corrected Images 2 Folder: ' target_path]);

tom_fit_ctf_folder(handles.pathname,handles.file_ext,target_path,ctf,method,corr_cut_off_pix);



function gui_cut_off_Callback(hObject, eventdata, handles)
% hObject    handle to gui_cut_off_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_cut_off_text as text
%        str2double(get(hObject,'String')) returns contents of gui_cut_off_text as a double


% --- Executes during object creation, after setting all properties.
function gui_cut_off_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_cut_off_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function handles=reload_it_man(handles,pathname,filename,coming_from_slider)

    handles.pathname=pathname;
    handles.filename=filename;
    [a b handles.file_ext]=fileparts(filename);
    set(handles.gui_filename,'String',handles.filename);
    set(handles.gui_file_path,'String',handles.pathname);
    [header] = tom_reademheader([pathname filename]);
    handles.magnification=header.Header.Magnification;
    handles.image_fullsize=header.Header.Size(1:2)';
    if header.Header.Objectpixelsize==0 && coming_from_slider==0
        
        prompt={'Objectpixelsize [Angstrom]:','Voltage [kV]:','Cs [mm]:','Intended Defocus [muem]:'};
        name='EM data';
        defaultanswer={'0','200','2.0','-3'};
        numlines=1;
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';
        answer=inputdlg(prompt,name,numlines,defaultanswer,options);
        
        if isempty(answer)
            disp('no EM data given.');
        else
            handles.pixelsize=str2num(answer{1});
            handles.voltage=str2num(answer{2}).*1000;
            handles.cs=str2num(answer{3});
            handles.intended_defocus=str2num(answer{4});
        end;
        
    else
        handles.pixelsize=header.Header.Objectpixelsize;
        handles.intended_defocus=header.Header.Defocus./10000;
        handles.cs=header.Header.Cs;
        handles.voltage=header.Header.Voltage;
    end;
    
    handles.defocus_fixed=handles.intended_defocus;
    handles.cc=2.2;
    in=tom_emreadc([pathname filename],'resample',[header.Header.Size(1)./512 header.Header.Size(2)./512 1]);
    axes(handles.main_axes);imagesc(in.Value');axis image;colormap gray;
    handles.preview_image_size=size(in.Value);
    set(handles.gui_magnification,'String',handles.magnification);
    set(handles.gui_pixelsize,'String',handles.pixelsize);
    set(handles.gui_voltage,'String',handles.voltage./1000);
    set(handles.gui_cs,'String',handles.cs);
    set(handles.gui_cc,'String',handles.cc);
    
    % Defocus
%    

    
    set(handles.gui_text_right_nyquist,'String',sprintf('%0.3g',handles.pixelsize.*2));
    set(handles.gui_text_right_2nyquist,'String',sprintf('%0.3g',handles.pixelsize.*4));
    set(handles.gui_text_left_nyquist,'String',sprintf('%0.3g',handles.pixelsize.*2));
    set(handles.gui_text_left_2nyquist,'String',sprintf('%0.3g',handles.pixelsize.*4));
    set(handles.gui_manual_text,'Visible','off');

    
    
    
    function handles=create_ps(handles)
        
        contents = get(handles.gui_Size,'String');
        sz=contents{get(handles.gui_Size,'Value')};
        if ~isempty(findstr(sz,'User'))
            handles.size=str2num(sz(2:8));
            handles.region(1)=str2num(sz(12:15));
            handles.region(2)=handles.region(1);
        else
            handles.size=str2num(sz);
            handles.region=handles.image_fullsize;
        end;
        
        
        
        if handles.region(1)<handles.image_fullsize(1)
            pos=round(ginput(1));
            
            pos_full=pos.*(handles.image_fullsize./handles.preview_image_size);
            if pos_full(1)<handles.region(1); pos_full(1)=1; end;
            if pos_full(2)<handles.region(2); pos_full(2)=1; end;
            if pos_full(1)+handles.region(1)>handles.image_fullsize(1); pos_full(1)=handles.image_fullsize(1)-(handles.region(1)+1); end;
            if pos_full(2)+handles.region(2)>handles.image_fullsize(2); pos_full(2)=handles.image_fullsize(2)-(handles.region(2)+1); end;

            in=tom_emreadc([handles.pathname handles.filename],'subregion',[pos_full 1],[handles.region-1 0]);
            disp(['calculate PS of subregion, position: ' num2str(pos_full) ' to: ' num2str(pos_full+handles.region)]);
            figure;
            tom_imagesc(single(in.Value));
        else        
            in=tom_emreadc([handles.pathname handles.filename]);        
        end;
        
        ps=tom_calc_periodogram_parallel(single(in.Value),handles.size(1),0,handles.size(1)./16);
        handles.untreated_ps=fftshift(ps);
        ps=(log(fftshift(ps)));
        %ps=((fftshift(ps)));
        
        %[decay decay_img]=calc_decay(ps,15,92,32);
        %ps=(double(ps-decay_img));
        handles.ps=ps;
        handles.mask_outer_radius=size(ps,1)./2-1;
        handles.correctbackground_outer_radius=size(ps,1)./2;
        set(handles.gui_background_mask_outer_radius,'String',handles.correctbackground_outer_radius);
        
        axes(handles.main_axes);imagesc(ps');axis image;colormap gray;
        
        
        function handles=correct_bg(handles)
            ps=handles.ps;
            if strcmp(handles.correctbackground_method,'calc_decay')==1
                if handles.mask_inner_radius-1>0
                    [handles.decay handles.decay_image]=calc_decay(ps,handles.correctbackground_inner_radius,handles.correctbackground_outer_radius,32);
                else
                    [handles.decay handles.decay_image]=calc_decay(ps,handles.correctbackground_inner_radius,handles.correctbackground_outer_radius,32);
                end;
            end
            handles.background_corrected_ps=double(ps-handles.decay_image);
            axes(handles.main_axes);imagesc(handles.background_corrected_ps');axis image;colormap gray;
            
            
            
            function handles=apply_mask_bg(handles)
                img_size=handles.size;
                mask_in_radius=handles.mask_inner_radius;
                mask_out_radius=handles.mask_outer_radius;
                mask_in = tom_spheremask(ones(img_size),mask_in_radius,0,[img_size(1)./2+1 img_size(2)./2+1 1]);
                mask_out = tom_spheremask(ones(img_size),mask_out_radius,0,[img_size(1)./2+1 img_size(2)./2+1 1]);
                mask=mask_out-mask_in;
                
                if isempty('handles.background_corrected_ps')
                    handles.ps=handles.untreated_ps.*mask;
                else
                    handles.ps=handles.background_corrected_ps.*mask;
                end;
                
                handles.mask=mask;
                drawnow;
                axes(handles.main_axes);imagesc(handles.ps');axis image;colormap gray;
                drawnow;
                
                
                function handles=fit_it(handles,eventdata,hObject)
                    
                    %nickell it!
                    handles.mask_outer_radius=str2double(get(handles.gui_mask_outer_radius,'String'));
                    handles=apply_mask_bg(handles);
                    
                    txt=get(handles.gui_perform_fit,'String');
                    set(handles.gui_perform_fit,'String','Searching ...');drawnow;
                    %set(handles.gui_perform_fit,'Enable','inactive')
                    [image EM Search]=get_gui_values(hObject, eventdata, handles);
                    
                    [Fit]=tom_fit_ctf(image,EM,Search);
                    axes(handles.main_axes);
                    [im_handle h_amplitude h_phase defocus_label_points_x defocus_label_points_y fit_ps]=show_ctf_fit(image, EM, Fit);
                    
                    handles.EM=EM;
                    handles.Search=Search;
                    handles.Fit=Fit;
                    try
                        handles.Fit.decay_image=handles.decay_image;
                        handles.Fit.untreated_ps=handles.untreated_ps;
                    catch
                    end;
                    handles.Fit.ps=handles.ps;
                    handles.Fit.fit_ps=fit_ps;
                    handles.Fit.decay_mask_inner_radius=handles.correctbackground_inner_radius;
                    handles.Fit.decay_mask_outer_radius=handles.correctbackground_outer_radius;
                    handles.Fit.EM=EM;
                    handles.Fit.defocus_label_points_x=defocus_label_points_x;
                    handles.Fit.defocus_label_points_y=defocus_label_points_y;
                    
                    EM
                    Search
                    Fit
                    hold on;
                    
                    defocus_vals_x=sqrt((defocus_label_points_x(:,1)-(size(handles.ps,1)./2)).^2+(defocus_label_points_x(:,2)-(size(handles.ps,1)./2)).^2)./(size(handles.ps,1)./2)./(handles.pixelsize.*2);
                    defocus_vals_y=sqrt((defocus_label_points_y(:,1)-(size(handles.ps,1)./2)).^2+(defocus_label_points_y(:,2)-(size(handles.ps,1)./2)).^2)./(size(handles.ps,1)./2)./(handles.pixelsize.*2);
                    handles.label_x=0;handles.label_y=0;
                    for i=1:size(defocus_vals_x,1)
                        lbl_x=sprintf('%0.3g',1./defocus_vals_x(i));
                        lbl_x=[' ' lbl_x];
                        th=text(defocus_label_points_x(i,1),defocus_label_points_x(i,2),lbl_x);
                        set(th,'Color','white');
                        set(th,'Fontsize',16);
                        handles.label_x(i)=th;
                    end;
                    for i=1:size(defocus_vals_y,1)
                        lbl_y=sprintf('%0.3g',1./defocus_vals_y(i));
                        lbl_y=[' ' lbl_y];
                        th2=text(defocus_label_points_y(i,1),defocus_label_points_y(i,2),lbl_y);
                        set(th2,'Color','white');
                        set(th2,'Fontsize',16);
                        handles.label_y(i)=th2;
                    end;
                    handles.label_x_dots=plot(defocus_label_points_x(:,1),defocus_label_points_x(:,2),'yo','MarkerSize',3,'Linewidth',5);
                    handles.label_y_dots=plot(defocus_label_points_y(:,1),defocus_label_points_y(:,2),'yo','MarkerSize',3,'Linewidth',5);
                    
                    hold off;
                    handles.h_amplitude=h_amplitude;
                    handles.h_phase=h_phase;
                    
                    set(handles.gui_defocus_fixed,'String',sprintf('%0.5g',Fit.Dz_det.*1e6));
                    
                    set(handles.gui_text_long_axis,'String',['Defocus long axis: ' sprintf('%0.5g',Fit.Dz_det.*1e6+(Fit.Dz_delta_det.*1e6)./2) ' mue.']);
                    set(handles.gui_text_short_axis,'String',['Defocus short axis: ' sprintf('%0.5g',Fit.Dz_det.*1e6-(Fit.Dz_delta_det.*1e6)./2) ' mue.']);
                    set(handles.gui_text_angle_axis,'String',['Angle (long&short axis): ' sprintf('%0.5g',Fit.Phi_0_det) ' deg.']);
                    
                    set(handles.Defocus_Slider,'Value',Fit.Dz_det.*1e6);
                    
                    set(handles.gui_defocus_delta_fixed,'String',sprintf('%0.5g',Fit.Dz_delta_det.*1e6));
                    set(handles.gui_Phi_0_fixed,'String',sprintf('%0.5g',Fit.Phi_0_det));
                    set(handles.gui_amplitude_contrast_fixed,'String',sprintf('%0.3g',Fit.amplitude_contrast_det));
                    set(handles.gui_decay_fixed,'String',sprintf('%0.3g',Fit.decay_det));
                    set(handles.gui_decay_part_coh_ill_fixed,'String',sprintf('%0.3g',Fit.decay_part_coh_ill_det));
                    set(handles.gui_energy_spread_fixed,'String',sprintf('%0.3g',Fit.decay_energy_spread_det));
                    
                    set(handles.gui_perform_fit,'String','Perform CTF Search');
                    set(handles.gui_perform_fit,'Enable','on');
                    drawnow;
                    set(handles.show_minima_red_button,'Value',1);
                    set(handles.show_maxima_green_button,'Value',1);
                    set(handles.show_zeros_yellow_button,'Value',1);
                    axes(handles.ccc_plot_axes);
                    r=reshape(Fit.corr_all,[numel(Fit.corr_all) 1]);
                    [n,xout] = hist(r,15); h=bar(xout,n,'red');
                    xlim([-1 1]);
                    ylim([0 max(n)]);
                    %set(handles.ccc_plot_axes,'XTickLabel',[-1:0.5:1]);
                    set(handles.ccc_plot_axes,'XGrid','on');
                    set(handles.ccc_plot_axes,'Color',[0 0 0])
                    set(handles.ccc_plot_axes,'YColor',[1 1 1]);
                    set(handles.ccc_plot_axes,'XColor',[1 1 1]);
                    hold on;
                    [me ma mi st va]=tom_dev(Fit.corr_all,'noinfo');
                    ht=text(max(xout),n(end),num2str(sprintf('%0.2g',max(r))));
                    if ma>me+st
                        set(ht,'FontSize',11,'Color','green','FontWeight','Bold');
                    end;
                    if (ma<me+st./2 || ma<0.5)
                        set(ht,'FontSize',11,'Color','yellow','FontWeight','Bold');
                    end;
                    if ma < 0.3
                        set(ht,'FontSize',11,'Color','red','FontWeight','Bold');
                        set(handles.ccc_plot_axes,'Color',[1 1 0])
                    end;
                    
                    hold off;
                    
                    
  function handles=set_gui_values(handles,Search,EM,image,Fit)                  
                 
      if (exist('Search','var'))
          %Defocus
          if (length(Search.Dz_search)>1)
              set(handles.gui_defocus_fixed_button,'Value',0);
              start=Search.Dz_search(1).*1e06;
              step=Search.Dz_search(2).*1e06-Search.Dz_search(1).*1e06;
              stop=Search.Dz_search(length(Search.Dz_search)).*1e06;
              set(handles.gui_defocus_search_start,'String',num2str(start)); 
              set(handles.gui_defocus_search_step,'String',num2str(step));    
              set(handles.gui_defocus_search_stop,'String',num2str(stop)); 
          else
              set(handles.gui_defocus_fixed_button,'Value',1);
              set(handles.gui_defocus_fixed,'String',num2str(Search.Dz_search(1).*1e06));
          end;
          
          %Astig Delta
          if (length(Search.Dz_delta_search)>1)
              set(handles.gui_defocus_delta_fixed_button,'Value',0);
              start=Search.Dz_delta_search(1).*1e06;
              step=Search.Dz_delta_search(2).*1e06-Search.Dz_delta_search(1).*1e06;
              stop=Search.Dz_delta_search(length(Search.Dz_delta_search)).*1e06;
              set(handles.gui_defocus_delta_search_start,'String',num2str(start)); 
              set(handles.gui_defocus_delta_search_step,'String',num2str(step));    
              set(handles.gui_defocus_delta_search_stop,'String',num2str(stop)); 
          else
              set(handles.gui_defocus_delta_fixed_button,'Value',1);
              set(handles.gui_defocus_delta_fixed,'String',num2str(Search.Dz_search(1).*1e06));
          end;
          
          %Astig Angle
          if (length(Search.Phi_0_search)>1)
              set(handles.gui_Phi_0_fixed_button,'Value',0);
              start=Search.Phi_0_search(1);
              step=Search.Phi_0_search(2)-Search.Phi_0_search(1);
              stop=Search.Phi_0_search(length(Search.Phi_0_search));
              set(handles.gui_Phi_0_search_start,'String',num2str(start)); 
              set(handles.gui_Phi_0_search_step,'String',num2str(step));    
              set(handles.gui_Phi_0_search_stop,'String',num2str(stop)); 
          else
              set(handles.gui_defocus_fixed_button,'Value',1);
              set(handles.gui_Phi_0_fixed,'String',num2str(Search.Phi_0_search(1)));
          end;
          
          % % Amplitude contrast, Ice is 0.07, negative stain is 0.15.
          if (length(Search.amplitude_contrast_search)>1)
              set(handles.gui_amplitude_contrast_button,'Value',0);
              start=Search.amplitude_contrast_search(1);
              step=Search.amplitude_contrast_search(2)-Search.amplitude_contrast_search(1);
              stop=Search.amplitude_contrast_search(length(Search.amplitude_contrast_search));
              set(handles.gui_ampltude_contrast_search_start,'String',num2str(start)); 
              set(handles.gui_ampltude_contrast_search_step,'String',num2str(step));    
              set(handles.gui_ampltude_contrast_search_stop,'String',num2str(stop)); 
          else
              set(handles.gui_amplitude_contrast_button,'Value',1);
              set(handles.gui_amplitude_contrast_fixed,'String',num2str(Search.amplitude_contrast_search(1)));
          end;
          
           % Decay k^2
          if (length(Search.decay_search)>1)
              set(handles.gui_decay_button,'Value',0);
              start=Search.decay_search(1);
              step=Search.decay_search(2)-Search.decay_search(1);
              stop=Search.decay_search(length(Search.decay_search));
              set(handles.gui_decay_search_start,'String',num2str(start)); 
              set(handles.gui_decay_search_step,'String',num2str(step));    
              set(handles.gui_decay_search_stop,'String',num2str(stop)); 
          else
              set(handles.gui_decay_button,'Value',1);
              set(handles.gui_decay_fixed,'String',num2str(Search.decay_search(1)));
          end;
          
         
           % Decay partially coherent illumination
          if (length(Search.decay_part_coh_ill_search)>1)
              set(handles.gui_decay_part_coh_ill_button,'Value',0);
              start=Search.decay_part_coh_ill_search(1);
              step=Search.decay_part_coh_ill_search(2)-Search.decay_part_coh_ill_search(1);
              stop=Search.decay_part_coh_ill_search(length(Search.decay_part_coh_ill_search));
              set(handles.gui_decay_part_coh_ill_search_start,'String',num2str(start)); 
              set(handles.gui_decay_part_coh_ill_search_step,'String',num2str(step));    
              set(handles.gui_decay_part_coh_ill_search_stop,'String',num2str(stop)); 
          else
              set(handles.gui_decay_part_coh_ill_button,'Value',1);
              set(handles.gui_decay_part_coh_ill_fixed,'String',num2str(Search.decay_part_coh_ill_search(1)));
          end;
          
       
           % Decay due to energy spread
          if (length(Search.decay_part_coh_ill_search)>1)
              set(handles.gui_decay_energy_spread_button,'Value',0);
              start=Search.decay_energy_spread_search(1);
              step=Search.decay_energy_spread_search(2)-Search.decay_energy_spread_search(1);
              stop=Search.decay_energy_spread_search(length(Search.decay_energy_spread_search));
              set(handles.gui_decay_energy_spread_search_start,'String',num2str(start)); 
              set(handles.gui_decay_energy_spread_search_step,'String',num2str(step));    
              set(handles.gui_decay_energy_spread_search_stop,'String',num2str(stop)); 
          else
              set(handles.gui_decay_energy_spread_button,'Value',1);
              set(handles.gui_energy_spread_fixed,'String',num2str(Search.decay_part_coh_ill_search(1)));
          end;
          
          
          set(handles.gui_mask_inner_radius,'String',num2str(Search.mask_inner_radius));
          set(handles.gui_mask_outer_radius,'String',num2str(Search.mask_outer_radius));
          set(handles.gui_background_mask_inner_radius,'String',num2str(Search.mask_inner_radius_bg));
          set(handles.gui_background_mask_outer_radius,'String',num2str(Search.mask_outer_radius_bg));
          
%           if (Search.ps_filter==3)
%               set(handles.gui_filter,'Value',2);  
%           else
%               set(handles.gui_filter,'Value',1);  
%           end;
          
          if (Search.ps_size==256)
              set(handles.gui_Size,'Value',1);
          else
              set(handles.gui_Size,'Value',2);
          end;
 end;
 
 
 if (exist('EM','var'))
     set(handles.gui_pixelsize,'String',num2str(EM.Objectpixelsize.*1e+10));
     set(handles.gui_voltage,'String',num2str(EM.Voltage./1000.0)); 
     set(handles.gui_cs,'String',num2str(EM.Cs.*1e+3));
     set(handles.gui_cc,'String',num2str(EM.Cc.*1e+3));
end;
 

if (exist('image','var') && isempty(image)==0 )
    handles.ps=image;
    axes(handles.main_axes);imagesc(image');axis image;colormap gray;
end;

if (exist('Fit','var'))
    set(handles.gui_defocus_fixed,'String',num2str(Fit.Dz_det.*1e06));
    set(handles.gui_defocus_delta_fixed,'String',num2str(Fit.Dz_delta_det.*1e06));
    set(handles.gui_Phi_0_fixed,'String',num2str(Fit.Phi_0_det));
    set(handles.gui_amplitude_contrast_fixed,'String',num2str(Fit.amplitude_contrast_det));
    set(handles.gui_decay_fixed,'String',num2str(Fit.decay_det));
    set(handles.gui_decay_part_coh_ill_fixed,'String',num2str(Fit.decay_part_coh_ill_det));
    set(handles.gui_energy_spread_fixed,'String',num2str(Fit.decay_energy_spread_det));
    if ((exist('image','var') &&  isempty(image)==0) && exist('EM','var') )
        [handles.im_handle handles.h_amplitude handles.h_phase handles.defocus_label_points_x handles.defocus_label_points_y handles.fit_ps]=show_ctf_fit(image, EM, Fit);
    end;
    if get(handles.show_minima_red_button,'Value')==1
        set(handles.h_phase,'Visible','on');
    else
        set(handles.h_phase,'Visible','off');
    end;
    
    if get(handles.show_maxima_green_button,'Value')==1
        set(handles.h_amplitude,'Visible','on');
    else
        set(handles.h_amplitude,'Visible','off');
    end;
    handles.Fit=Fit;
end;





                    
                    



function high2lowdiff_Callback(hObject, eventdata, handles)
% hObject    handle to high2lowdiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of high2lowdiff as text
%        str2double(get(hObject,'String')) returns contents of high2lowdiff as a double


% --- Executes during object creation, after setting all properties.
function high2lowdiff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to high2lowdiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_corr_meth.
function pop_corr_meth_Callback(hObject, eventdata, handles)
% hObject    handle to pop_corr_meth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_corr_meth contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_corr_meth


% --- Executes during object creation, after setting all properties.
function pop_corr_meth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_corr_meth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mon_checked.
function mon_checked_Callback(hObject, eventdata, handles)
% hObject    handle to mon_checked (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mon_checked


% --- Executes on button press in mon_accepted.
function mon_accepted_Callback(hObject, eventdata, handles)
% hObject    handle to mon_accepted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mon_accepted


% --- Executes on button press in mon_selected.
function mon_selected_Callback(hObject, eventdata, handles)
% hObject    handle to mon_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mon_selected


% --- Executes on button press in filter_fit.
function filter_fit_Callback(hObject, eventdata, handles)
% hObject    handle to filter_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)








% --- Executes on button press in fl_sl_browse.
function fl_sl_browse_Callback(hObject, eventdata, handles)
% hObject    handle to fl_sl_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fl_sl_browse
set(handles.mon_accepted,'Visible','off');
set(handles.mon_checked,'Visible','off');
%set(handles.mon_selected,'Visible','off');

val=round(get(handles.gui_file_slider,'Value'));

[a b c]=fileparts(handles.folder_filenames{val});
set(handles.gui_filename,'String',[b c]);

handles=reload_it_man(handles,handles.pathname,[b c],0);  drawnow;

guidata(hObject, handles);


% --- Executes on button press in fl_sl_fit.
function fl_sl_fit_Callback(hObject, eventdata, handles)
% hObject    handle to fl_sl_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fl_sl_fit
set(handles.mon_accepted,'Visible','off');
set(handles.mon_checked,'Visible','off');
%set(handles.mon_selected,'Visible','off');

val=round(get(handles.gui_file_slider,'Value'));

[a b c]=fileparts(handles.folder_filenames{val});
set(handles.gui_filename,'String',[b c]);

handles=reload_it_man(handles,handles.pathname,[b c],0);  drawnow;
%handles=create_ps(handles);drawnow;
%handles=correct_bg(handles); drawnow;
%handles.mask_outer_radius=str2double(get(handles.gui_mask_outer_radius,'String'));
%handles=apply_mask_bg(handles);  drawnow;
%handles=fit_it(handles,eventdata,hObject);  drawnow;


guidata(hObject, handles);


% --- Executes on button press in fl_sl_monitor.
function fl_sl_monitor_Callback(hObject, eventdata, handles)
% hObject    handle to fl_sl_monitor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fl_sl_monitor

set(handles.mon_accepted,'Visible','on');
set(handles.mon_checked,'Visible','on');
%set(handles.mon_selected,'Visible','on');

val=round(get(handles.gui_file_slider,'Value'));

[a b c]=fileparts(handles.folder_filenames{val});
set(handles.gui_filename,'String',[b c]);


if (exist([handles.folder_filenames{val} '.mat'],'file'))
    
    if (isempty(handles.old_name)==0)
        [image EM Search]=get_gui_values(hObject, eventdata, handles);
        st_out.img=handles.ps;
        st_out.Fit=handles.Fit;
        st_out.Fit.Dz_det=str2num(get(handles.gui_defocus_fixed,'String')).*1e-06;
        st_out.Fit.EM=EM;
        st_out.Search=Search;
        save(handles.old_name,'st_out');
    end;
    load([handles.folder_filenames{val} '.mat']);
    handles.old_name=[handles.folder_filenames{val} '.mat'];
    handles.ctf_st=st_out;
    handles=set_gui_values(handles,st_out.Search,st_out.Fit.EM,st_out.img,st_out.Fit);
end;


guidata(hObject, handles);



function gui_cut_off_edit_Callback(hObject, eventdata, handles)
% hObject    handle to gui_cut_off_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_cut_off_edit as text
%        str2double(get(hObject,'String')) returns contents of gui_cut_off_edit as a double


% --- Executes during object creation, after setting all properties.
function gui_cut_off_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_cut_off_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on fl_sl_browse and none of its controls.
function fl_sl_browse_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to fl_sl_browse (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
