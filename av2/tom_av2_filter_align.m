function varargout = tom_av2_filter_align(varargin)
%  tom_av2_filter_align filters a alignment strucure and displays the
%  result
%  
%    tom_av2_filter_align(align2d,field,filter,binning) 
%  
%  PARAMETERS
%  
%    INPUT
%     align2d             (opt.) filename or struct  
%     field               (ccc) field 2 filter
%     filter              (3) filter kernel
%     binning             (2) binning of the displayed image
%
%    
%    OUTPUT
%
%  
%  EXAMPLE
%    
%
%  
%  REFERENCES
%  
%  SEE ALSO
%     ...
%  
%     created by FB 12/14/11
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_av2_filter_align_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_filter_align_OutputFcn, ...
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


% --- Executes just before tom_av2_filter_align is made visible.
function tom_av2_filter_align_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_av2_filter_align (see VARARGIN)

% Choose default command line output for tom_av2_filter_align
handles.output = hObject;

%init
handles.data.points_h=-1;

%transfer values
if (length(varargin) < 1)
    [file path]=uigetfile('*.mat');
    handles.data.f_align2d=[path '/' file];
else
    handles.data.f_align2d=varargin{1};
end;

if (length(varargin) < 2)
    ini_field='ccc';
else
    ini_field=varargin{2};
end;

if (length(varargin) < 3)
    handles.data.filter_v=3;
else
    handles.data.filter_v=varargin{3};
end;

if (length(varargin) < 4)
   handles.data.bin=2;
else
   handles.data.bin=varargin{4};
end;


set(handles.edit_filter_ker,'String',num2str(handles.data.filter_v));

if (isstruct(handles.data.f_align2d)==0)
    load(handles.data.f_align2d);
   
else
    align2d=handles.data.f_align2d;
end;

handles.data.align2d=align2d;
handles.data.fields=fieldnames(align2d);
handles.data.filenames_unique=unique({align2d(1,:).filename});
handles.data.filenames_all={align2d(1,:).filename};

handles=change_field(handles,ini_field);

[a b c]=fileparts(handles.data.filenames_unique{1});
set(handles.edit_file_path,'String',a);
set(handles.edit_img,'String',[b c]);

[tmp_val tmp_coord points_h]=render_image(handles);
handles.data.tmp_coord=tmp_coord;
handles.data.tmp_value=tmp_val;
handles.data.points_h=points_h;

num_of_files=length(handles.data.filenames_unique);
set(handles.sl_img,'Max',num_of_files);
set(handles.sl_img,'Min',1);
set(handles.sl_img,'Value',1);
set(handles.sl_img,'SliderStep',[(1./(num_of_files-1)) (1./(num_of_files-1))]);

% set(handles.sl_img,'SliderStep',[1 1]);
% set(handles.sl_img,'Max',length(handles.data.filenames_unique)-1);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_av2_filter_align wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_av2_filter_align_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in pop_fileds.
function pop_fileds_Callback(hObject, eventdata, handles)
% hObject    handle to pop_fileds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_fileds contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_fileds

all_fields=get(handles.pop_fileds,'String');
sel_field=all_fields{get(handles.pop_fileds,'Value')};

handles=change_field(handles,sel_field);
[tmp_val tmp_coord points_h]=render_image(handles);
handles.data.tmp_coord=tmp_coord;
handles.data.tmp_value=tmp_val;
handles.data.points_h=points_h;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pop_fileds_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_fileds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sl_img_Callback(hObject, eventdata, handles)
% hObject    handle to sl_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%nr=1+((length(handles.data.filenames_unique)-1).*get(handles.sl_img,'Value'));

nr=get(handles.sl_img,'Value');

f_name_tmp=handles.data.filenames_unique{nr};

[a b c]=fileparts(f_name_tmp);

set(handles.edit_img,'String',[b c]);
set(handles.edit_file_path,'String',a);

[tmp_val tmp_coord points_h]=render_image(handles);
handles.data.tmp_coord=tmp_coord;
handles.data.tmp_value=tmp_val;
handles.data.points_h=points_h;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sl_img_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_val_Callback(hObject, eventdata, handles)
% hObject    handle to edit_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_val as text
%        str2double(get(hObject,'String')) returns contents of edit_val as a double

handles.data.points_h=render_points(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sl_val_Callback(hObject, eventdata, handles)
% hObject    handle to sl_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

act_val=get(handles.sl_val,'Value').*(handles.data.val_min_max(2)- handles.data.val_min_max(1))+handles.data.val_min_max(1);

set(handles.edit_val,'String',num2str(act_val));

handles.data.points_h=render_points(handles);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sl_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_img_Callback(hObject, eventdata, handles)
% hObject    handle to edit_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_img as text
%        str2double(get(hObject,'String')) returns contents of edit_img as a double


% --- Executes during object creation, after setting all properties.
function edit_img_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sel_field_Callback(hObject, eventdata, handles)
% hObject    handle to sel_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sel_field as text
%        str2double(get(hObject,'String')) returns contents of sel_field as a double


% --- Executes during object creation, after setting all properties.
function sel_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sel_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_filter_ker_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filter_ker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filter_ker as text
%        str2double(get(hObject,'String')) returns contents of edit_filter_ker as a double


[tmp_val tmp_coord points_h]=render_image(handles);
handles.data.tmp_coord=tmp_coord;
handles.data.tmp_value=tmp_val;
handles.data.points_h=points_h;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edit_filter_ker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filter_ker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%helper functions

function v_min_max=min_max(align2d,field)

call=['vect=[align2d(1,:).' field '];'];
eval(call);

v_min_max(1)=min(vect);
v_min_max(2)=max(vect);
v_min_max(3)=mean(vect);



function [tmp_val tmp_coord points_h]=render_image(handles)


f_path=get(handles.edit_file_path,'String');
f_name=get(handles.edit_img,'String');

im=tom_emreadc([f_path '/' f_name],'binning',handles.data.bin);
im=tom_filter(im.Value,str2double(get(handles.edit_filter_ker,'String')));
axes(handles.ax_img);
imagesc(im'); colormap gray; axis image; 
set(handles.ax_img,'XTick',[]); set(handles.ax_img,'YTick',[]);
idx_tmp=find(ismember(handles.data.filenames_all,[f_path '/' f_name]));

tmp_field=get(handles.sel_field,'String');
tmp_coord=zeros(2,length(idx_tmp));
tmp_val=zeros(1,length(idx_tmp));
scal=2^handles.data.bin;
for i=1:length(idx_tmp)
    tmp_coord(:,i)=[(handles.data.align2d(1,idx_tmp(i)).position.x./scal) (handles.data.align2d(1,idx_tmp(i)).position.y./scal)]; 
    tmp_val(i)=getfield(handles.data.align2d(1,idx_tmp(i)),tmp_field);
end;
handles.data.tmp_coord=tmp_coord;
handles.data.tmp_value=tmp_val;
handles.data.points_h=-1;
points_h=render_points(handles);



function points_h=render_points(handles)

if (handles.data.points_h~=-1)
    delete(handles.data.points_h);
end;

thr=str2double(get(handles.edit_val,'String'));
tmp_idx=find(handles.data.tmp_value>thr);
hold on; points_h=plot(handles.data.tmp_coord(1,tmp_idx),handles.data.tmp_coord(2,tmp_idx),'ro'); hold off;
set(handles.num_part,'String',num2str(length(tmp_idx)));
set(handles.mean_cc,'String',num2str(mean(handles.data.tmp_value(tmp_idx))));



function num_part_Callback(hObject, eventdata, handles)
% hObject    handle to num_part (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_part as text
%        str2double(get(hObject,'String')) returns contents of num_part as a double


% --- Executes during object creation, after setting all properties.
function num_part_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_part (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mean_cc_Callback(hObject, eventdata, handles)
% hObject    handle to mean_cc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_cc as text
%        str2double(get(hObject,'String')) returns contents of mean_cc as a double


% --- Executes during object creation, after setting all properties.
function mean_cc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_cc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_file_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_file_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_file_path as text
%        str2double(get(hObject,'String')) returns contents of edit_file_path as a double


% --- Executes during object creation, after setting all properties.
function edit_file_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_file_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles=change_field(handles,ini_field)

selected_idx=find(ismember(handles.data.fields,ini_field));
set(handles.pop_fileds,'String',handles.data.fields);
set(handles.pop_fileds,'Value',selected_idx);
set(handles.sel_field,'String',handles.data.fields{selected_idx});
handles.data.val_min_max=min_max(handles.data.align2d,ini_field);
set(handles.edit_val,'String',num2str(handles.data.val_min_max(3)));
set(handles.sl_val,'Value',0.5);
