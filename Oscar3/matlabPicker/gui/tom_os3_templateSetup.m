function varargout = tom_os3_templateSetup(varargin)
% TOM_OS3_TEMPLATESETUP M-file for tom_os3_templateSetup.fig
%      TOM_OS3_TEMPLATESETUP, by itself, creates a new TOM_OS3_TEMPLATESETUP or raises the existing
%      singleton*.
%
%      H = TOM_OS3_TEMPLATESETUP returns the handle to a new TOM_OS3_TEMPLATESETUP or the handle to
%      the existing singleton*.
%
%      TOM_OS3_TEMPLATESETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_OS3_TEMPLATESETUP.M with the given input arguments.
%
%      TOM_OS3_TEMPLATESETUP('Property','Value',...) creates a new TOM_OS3_TEMPLATESETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_os3_templateSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_os3_templateSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_os3_templateSetup

% Last Modified by GUIDE v2.5 19-Nov-2008 14:29:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_os3_templateSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_os3_templateSetup_OutputFcn, ...
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


% --- Executes just before tom_os3_templateSetup is made visible.
function tom_os3_templateSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_os3_templateSetup (see VARARGIN)

% Choose default command line output for tom_os3_templateSetup
handles.output = hObject;
handles.single_stack_out_flag=0;

if length(varargin) > 0
    if (strcmp(varargin{1},'single_stack'))
        handles.single_stack_out_flag=1;
    end;
    if  length(varargin) > 1 
        handles.data.size=varargin{2};
    end;
end;


 
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_os3_templateSetup wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_os3_templateSetup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

close(handles.figure1);

% --- Executes on button press in mod_br_path.
function mod_br_path_Callback(hObject, eventdata, handles)
% hObject    handle to mod_br_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file path]=uigetfile('*.em;*.pdb');

try
    em_dens=tom_emread([path file]);
    em_dens=em_dens.Value;
    
    sz_org=size(em_dens);
    prompt = {'Enter Model pixelsize:','Enter micrograph pixelsize:','Invert Contrast (-1 for invert)'};
    dlg_title = 'Scaling of model';
    num_lines = 1;
    def = {'1','1','1'}; %try to pre fill
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    %calculate scaling factor
    fact=str2double(answer{1})./str2double(answer{2});
    inf_fact=str2double(answer{3});
    sz_rss=round(sz_org.*fact);
    
    prompt = {'Cubic length'};
    dlg_title = 'Enter Cubic size';
    num_lines = 1;
    sz_new=ceil(sz_rss./16).*16;
    def = {num2str(sz_new(1))}; %try to pre fill
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    modl=tom_rescale3d(em_dens,sz_rss).*inf_fact;
    %modl=tom_norm(modl,'mean0+1std');
    sz_new=[str2num(answer{1}) str2num(answer{1}) str2num(answer{1})];
    modl=tom_paste2(zeros(sz_new),modl,floor((sz_new-sz_rss)./2)+1);
    tom_emwrite([path '/' file '_rs.em'],modl);
    name=[path '/' file '_rs.em'];
catch
   em_dens=tom_pdbread([path file]);
   prompt = {'Enter Micrograph pixelsize: in A','Enter Cube size:'};
   dlg_title = 'Scaling of model';
   num_lines = 1;
   def = {'1','64 64 64'}; %try to pre fill
   answer = inputdlg(prompt,dlg_title,num_lines,def);
   em_dens=tom_pdb2em(em_dens,pixels,cube_size);
   [path file ext]=fileparts(file);
   tom_emwrite([path '/' file '.em'],em_dens);
   name=[path file ext];
end;

set(handles.mod_path,'String',name);

function mod_path_Callback(hObject, eventdata, handles)
% hObject    handle to mod_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mod_path as text
%        str2double(get(hObject,'String')) returns contents of mod_path as a double


% --- Executes during object creation, after setting all properties.
function mod_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mod_theta_end_Callback(hObject, eventdata, handles)
% hObject    handle to mod_theta_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mod_theta_end as text
%        str2double(get(hObject,'String')) returns contents of mod_theta_end as a double


% --- Executes during object creation, after setting all properties.
function mod_theta_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_theta_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mod_psi_end_Callback(hObject, eventdata, handles)
% hObject    handle to mod_psi_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mod_psi_end as text
%        str2double(get(hObject,'String')) returns contents of mod_psi_end as a double


% --- Executes during object creation, after setting all properties.
function mod_psi_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_psi_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mod_phi_end_Callback(hObject, eventdata, handles)
% hObject    handle to mod_phi_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mod_phi_end as text
%        str2double(get(hObject,'String')) returns contents of mod_phi_end as a double


% --- Executes during object creation, after setting all properties.
function mod_phi_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_phi_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mod_theta_incre_Callback(hObject, eventdata, handles)
% hObject    handle to mod_theta_incre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mod_theta_incre as text
%        str2double(get(hObject,'String')) returns contents of mod_theta_incre as a double


% --- Executes during object creation, after setting all properties.
function mod_theta_incre_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_theta_incre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mod_psi_incre_Callback(hObject, eventdata, handles)
% hObject    handle to mod_psi_incre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mod_psi_incre as text
%        str2double(get(hObject,'String')) returns contents of mod_psi_incre as a double


% --- Executes during object creation, after setting all properties.
function mod_psi_incre_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_psi_incre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mod_phi_incre_Callback(hObject, eventdata, handles)
% hObject    handle to mod_phi_incre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mod_phi_incre as text
%        str2double(get(hObject,'String')) returns contents of mod_phi_incre as a double


% --- Executes during object creation, after setting all properties.
function mod_phi_incre_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_phi_incre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mod_theta_start_Callback(hObject, eventdata, handles)
% hObject    handle to mod_theta_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mod_theta_start as text
%        str2double(get(hObject,'String')) returns contents of mod_theta_start as a double


% --- Executes during object creation, after setting all properties.
function mod_theta_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_theta_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mod_psi_start_Callback(hObject, eventdata, handles)
% hObject    handle to mod_psi_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mod_psi_start as text
%        str2double(get(hObject,'String')) returns contents of mod_psi_start as a double


% --- Executes during object creation, after setting all properties.
function mod_psi_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_psi_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mod_phi_start_Callback(hObject, eventdata, handles)
% hObject    handle to mod_phi_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mod_phi_start as text
%        str2double(get(hObject,'String')) returns contents of mod_phi_start as a double


% --- Executes during object creation, after setting all properties.
function mod_phi_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_phi_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in mod_proj_progr.
function mod_proj_progr_Callback(hObject, eventdata, handles)
% hObject    handle to mod_proj_progr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns mod_proj_progr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mod_proj_progr



% --- Executes during object creation, after setting all properties.
function mod_proj_progr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_proj_progr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mod_view.
function mod_view_Callback(hObject, eventdata, handles)
% hObject    handle to mod_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
st=get_gui_values(handles);
im=tom_emread(st.threed_mod.mod_path);
figure; tom_dspcub(im.Value);



% --- Executes on button press in part_stack_br.
function part_stack_br_Callback(hObject, eventdata, handles)
% hObject    handle to part_stack_br (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file path]=uigetfile('*.em');

set(handles.part_stack_path,'String',[path '/' file]);



function part_stack_path_Callback(hObject, eventdata, handles)
% hObject    handle to part_stack_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of part_stack_path as text
%        str2double(get(hObject,'String')) returns contents of part_stack_path as a double


% --- Executes during object creation, after setting all properties.
function part_stack_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to part_stack_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in part_stack_view.
function part_stack_view_Callback(hObject, eventdata, handles)
% hObject    handle to part_stack_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in part_stack_gen.
function part_stack_gen_Callback(hObject, eventdata, handles)
% hObject    handle to part_stack_gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    tmp=tom_av2_particlepickergui('stackname_out',round(handles.data.size(1)./2));
catch
    tmp=tom_av2_particlepickergui('stackname_out');
end;

set(handles.tmp_path,'String', tmp.stackname);


% --- Executes on button press in temp_tmp_br.
function temp_tmp_br_Callback(hObject, eventdata, handles)
% hObject    handle to temp_tmp_br (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
st=get_gui_values(handles);

[file path]=uigetfile('*.em');

set(handles.tmp_path,'String',[path '/' file]);


function tmp_path_Callback(hObject, eventdata, handles)
% hObject    handle to tmp_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmp_path as text
%        str2double(get(hObject,'String')) returns contents of tmp_path as a double





% --- Executes during object creation, after setting all properties.
function tmp_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmp_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in temp_tmp_view.
function temp_tmp_view_Callback(hObject, eventdata, handles)
% hObject    handle to temp_tmp_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
st=get_gui_values(handles);
stack=tom_emread(st.temp_st.temp_ini_path);
if (size(stack.Value,3)==1)
    figure; tom_imagesc(stack.Value);
else
    figure; tom_dspcub(stack);
end;

% --- Executes on button press in temp_tmp_gen.
function temp_tmp_gen_Callback(hObject, eventdata, handles)
% hObject    handle to temp_tmp_gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
st=get_gui_values(handles);
stack=tom_emread(st.temp_st.tmp_path);

mask=tom_spheremask(ones(stack.Header.Size(1),stack.Header.Size(2)),round(stack.Header.Size(1)./2)-1);
if (isempty(st.temp_st.temp_ini_path))   
    tmp_tmpl=stack.Value(:,:,1);
else
    tmp_tmpl=tom_emread(st.temp_st.temp_ini_path);
    tmp_tmpl=tmp_tmpl.Value;
end;
[stack_out all_tmpl cc_out]= tom_os3_alignStack2(stack.Value,tmp_tmpl,mask,0,'mean0+1std',[4 3],'default');
if (isempty(st.temp_st.temp_ini_path))   
    [file path]=uiputfile('*.em');
    tom_emwrite([path file],all_tmpl{4}); 
    set(handles.temp_ini_path,'String',[path file]);
else
    tom_emwrite(st.temp_st.temp_ini_path,all_tmpl{4}); 
end;


% --- Executes on button press in templ_reduce.
function templ_reduce_Callback(hObject, eventdata, handles)
% hObject    handle to templ_reduce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
st=get_gui_values(handles);
stack=tom_emread(st.temp_st.tmp_path);
stack=stack.Value;

if (isempty(st.temp_st.temp_ini_path))
    template='';
else
    template=tom_emread(st.temp_st.temp_ini_path); 
    template=template.Value;
end;

classes=tom_os3_classify_templStack(stack,template,str2num(st.temp_st.temp_num_of_cent),str2num(st.temp_st.temp_bin_alignment),str2num(st.temp_st.tmpl_binning_classification),st.temp_st.temp_programm,str2num(st.temp_st.tmp_num_of_rot));

[file path]=uiputfile('*.em');

set(handles.red_stack_path,'String',[path '/' file]);

tom_emwrite([path '/' file],classes);




% --- Executes on selection change in temp_programm.
function temp_programm_Callback(hObject, eventdata, handles)
% hObject    handle to temp_programm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns temp_programm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from temp_programm


% --- Executes during object creation, after setting all properties.
function temp_programm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temp_programm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function temp_num_of_cent_Callback(hObject, eventdata, handles)
% hObject    handle to temp_num_of_cent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of temp_num_of_cent as text
%        str2double(get(hObject,'String')) returns contents of temp_num_of_cent as a double


% --- Executes during object creation, after setting all properties.
function temp_num_of_cent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temp_num_of_cent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in red_br_stack.
function red_br_stack_Callback(hObject, eventdata, handles)
% hObject    handle to red_br_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file path]=uigetfile('*.em');

set(handles.red_stack_path,'String',[path '/' file]);



function red_stack_path_Callback(hObject, eventdata, handles)
% hObject    handle to red_stack_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of red_stack_path as text
%        str2double(get(hObject,'String')) returns contents of red_stack_path as a double


% --- Executes during object creation, after setting all properties.
function red_stack_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_stack_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in red_view.
function red_view_Callback(hObject, eventdata, handles)
% hObject    handle to red_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
st=get_gui_values(handles);
stack=tom_emread(st.red_stack.red_stack_path);
if (size(stack.Value,3)==1)
    figure; tom_imagesc(stack.Value);
else
    figure; tom_dspcub(stack);
end;


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in red_sort.
function red_sort_Callback(hObject, eventdata, handles)
% hObject    handle to red_sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);

stack=tom_emread(st.red_stack.red_stack_path);

[base filename ext]=fileparts(st.red_stack.red_stack_path);

outfilestruct.method='singleFiles'; 
outfilestruct.num_of_entries=1000; 
outfilestruct.path=[base '/part_fold'];
outfilestruct.folder_num_offset=1;
outfilestruct.filename='part_';
outfilestruct.ext='.em';
outfilestruct.fileformat='em';

tom_av2_stack_clean(outfilestruct,'outfilestruct');
tom_av2_stack_write(stack.Value,outfilestruct,1); 
save([outfilestruct.path '.mat'],'outfilestruct');
tmpp{1}=[outfilestruct.path '.mat'];
outfilestruct=tom_av2_stackbrowser(tmpp,'makerefmode');
load(outfilestruct);
pause2(0.5);
stack=tom_av2_stack_read(outfilestruct,outfilestruct.idx.written_index);
tom_av2_stack_clean(outfilestruct,'outfilestruct');
tom_emwrite(st.red_stack.red_stack_path,stack);



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function temp_ini_path_Callback(hObject, eventdata, handles)
% hObject    handle to temp_ini_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of temod_pathmp_ini_path as text
%        str2double(get(hObject,'String')) returns contents of temp_ini_path as a double


% --- Executes during object creation, after setting all properties.
function temp_ini_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temp_ini_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in temp_tmp_mod.
function temp_tmp_mod_Callback(hObject, eventdata, handles)
% hObject    handle to temp_tmp_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
st=get_gui_values(handles);
tmpl=tom_emread(st.temp_st.temp_ini_path);
new=tom_av2_alignref(tmpl.Value);
tom_emwrite(st.temp_st.temp_ini_path,new); 

% --- Executes on button press in temp_st_br.
function temp_st_br_Callback(hObject, eventdata, handles)
% hObject    handle to temp_st_br (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file path]=uigetfile('*.em');
set(handles.temp_ini_path,'String',[path '/' file]);

% --- Executes on button press in temp_st_view.
function temp_st_view_Callback(hObject, eventdata, handles)
% hObject    handle to temp_st_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);
stack=tom_emread(st.temp_st.tmp_path);
figure; tom_dspcub(stack);


% --- Executes on button press in mod_project.
function mod_project_Callback(hObject, eventdata, handles)
% hObject    handle to mod_project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)1

st=get_gui_values(handles);

mod=tom_emreadc(st.threed_mod.mod_path);
sz_mod=mod.Header.Size;
mask=tom_spheremask(ones(sz_mod),round(sz_mod(1)./2)-2,0);
%mod=tom_norm(mod.Value,'mean0+1std',mask).*mask;
mod=mod.Value;

phi=st.threed_mod.phi;
psi=st.threed_mod.psi;
theta=st.threed_mod.theta;


if (sum(phi==[0 0 0])==3)
    phi=[1 1 1];
end;

if (sum(psi==[0 0 0])==3)
    psi=[1 1 1];
end;

if (sum(theta==[0 0 0])==3)
    theta=[1 1 1];
end;



l_proj_st=length(phi(1):phi(2):phi(3))*length(psi(1):psi(2):psi(3))*length(theta(1):theta(2):theta(3));

proj=zeros(size(mod,1),size(mod,2),l_proj_st);
ang=zeros(3,l_proj_st);

% h = waitbar(0,'Projecting...');
% zz=1;
% for phi_iter=phi(1):phi(2):phi(3)
%    for psi_iter=psi(1):psi(2):psi(3)
%        for theta_iter=theta(1):theta(2):theta(3)
%             ang=tom_sum_rotation([0 psi_iter theta_iter; 270 90 phi_iter],[0 0 0; 0 0 0]);
%             proj(:,:,zz)=sum(tom_rotate(mod,ang),2);
%             zz=zz+1;
%             waitbar(zz/l_proj_st,h)
%         end;
%    end;
% end;
% close(h);



zz=1;
for phi_iter=phi(1):phi(2):phi(3)
   for psi_iter=psi(1):psi(2):psi(3)
       for theta_iter=theta(1):theta(2):theta(3)
            ang(:,zz)=tom_sum_rotation([0 psi_iter theta_iter; 270 90 phi_iter],[0 0 0; 0 0 0]);
            zz=zz+1;
        end;
   end;
end;

parfor iii=1:zz-1
    proj(:,:,iii)=sum(tom_rotate(mod,ang(:,iii)'),2);
    disp(iii);
end;


if (st.threed_mod.mod_append_proj==1 && isempty(st.temp_st.tmp_path)==0)
    proj_old=tom_emreadc3(st.temp_st.tmp_path); proj_old=proj_old.Value;
    proj=cat(3,proj_old,proj);
    out_name=st.temp_st.tmp_path;
else
    [file path]=uiputfile('*.em');
    out_name=[path '/' file];
end;

tom_emwrite(out_name,proj);
set(handles.tmp_path,'String',out_name);




function st=get_gui_values(handles)

st.threed_mod.mod_path=get(handles.mod_path,'String');
st.threed_mod.phi=[str2num(get(handles.mod_phi_start,'String')) str2num(get(handles.mod_phi_incre,'String')) str2num(get(handles.mod_phi_end,'String'))];
st.threed_mod.psi=[str2num(get(handles.mod_psi_start,'String')) str2num(get(handles.mod_psi_incre,'String')) str2num(get(handles.mod_psi_end,'String'))];
st.threed_mod.theta=[str2num(get(handles.mod_theta_start,'String')) str2num(get(handles.mod_theta_incre,'String')) str2num(get(handles.mod_theta_end,'String'))];
tmp=get(handles.mod_proj_progr,'String');
st.threed_mod.mod_proj_progr=tmp{get(handles.mod_proj_progr,'Value')};
st.threed_mod.mod_append_proj=get(handles.mod_append_proj,'Value');
%st.temp_st.part_stack_path=get(handles.part_stack_path,'String');

st.temp_st.temp_ini_path=get(handles.temp_ini_path,'String');
st.temp_st.tmp_path=get(handles.tmp_path,'String');
st.temp_st.temp_num_of_cent=get(handles.temp_num_of_cent,'String');
tmp=get(handles.temp_programm,'String');
st.temp_st.temp_programm=tmp{get(handles.temp_programm,'Value')};
st.temp_st.tmpl_binning_classification=get(handles.tmpl_binning_classification,'String');
st.temp_st.temp_bin_alignment=get(handles.temp_bin_alignment,'String');
st.temp_st.tmp_num_of_rot=get(handles.tmp_num_of_rot,'String');
st.red_stack.red_stack_path=get(handles.red_stack_path,'String');


% --- Executes on button press in red_mod.
function red_mod_Callback(hObject, eventdata, handles)
% hObject    handle to red_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB% eu=tom_sum_rotation([0 0 theta_iter; 270 90 psi_iter],[0 0 0 ; 0 0 0]);
% handles    structure with handles and user data (see GUIDATA)
st=get_gui_values(handles);
stack=tom_emread(st.red_stack.red_stack_path);
stack_out=tom_av2_alignref(stack.Value);
tom_emwrite(st.red_stack.red_stack_path,stack_out);



function tmp_num_of_rot_Callback(hObject, eventdata, handles)
% hObject    handle to tmp_num_of_rot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmp_num_of_rot as text
%        str2double(get(hObject,'String')) returns contents of tmp_num_of_rot as a double


% --- Executes during object creation, after setting all properties.
function tmp_num_of_rot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmp_num_of_rot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function temp_bin_alignment_Callback(hObject, eventdata, handles)
% hObject    handle to temp_bin_alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of temp_bin_alignment as text
%        str2double(get(hObject,'String')) returns contents of temp_bin_alignment as a double


% --- Executes during object creation, after setting all properties.
function temp_bin_alignment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temp_bin_alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in exit_push.
function exit_push_Callback(hObject, eventdata, handles)
% hObject    handle to exit_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


st=get_gui_values(handles);
if (handles.single_stack_out_flag==0)
    handles=rmfield(handles,'output');
    path=uigetdir('Create final template Stack');
    handles.output.path=path;
    
    tmp=tom_emreadc(st.red_stack.red_stack_path);
    outfilestruct.method='singleFiles';
    outfilestruct.num_of_entries=10000; %just for the show ...use 10000 for real application hurz!
    outfilestruct.path=path;
    outfilestruct.folder_num_offset=-1;
    outfilestruct.filename='template_';
    outfilestruct.ext='.em';
    outfilestruct.fileformat='em';
    % outfilestruct.size=[80 80 50];
    tom_av2_stack_clean(outfilestruct);
    tom_av2_stack_write(tmp.Value,outfilestruct,1);
    
else
   handles=rmfield(handles,'output');
   handles.output.path=st.red_stack.red_stack_path;
end;


guidata(hObject, handles);
uiresume(handles.figure1);



function tmpl_binning_classification_Callback(hObject, eventdata, handles)
% hObject    handle to tmpl_binning_classification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmpl_binning_classification as text
%        str2double(get(hObject,'String')) returns contents of tmpl_binning_classification as a double


% --- Executes during object creation, after setting all properties.
function tmpl_binning_classification_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmpl_binning_classification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mod_append_proj.
function mod_append_proj_Callback(hObject, eventdata, handles)
% hObject    handle to mod_append_proj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mod_append_proj


% --- Executes on button press in mod_volxyz.
function mod_volxyz_Callback(hObject, eventdata, handles)
% hObject    handle to mod_volxyz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


