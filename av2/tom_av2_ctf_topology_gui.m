function varargout = tom_av2_ctf_topology_gui(varargin)
% TOM_AV2_CTF_TOPOLOGY_GUI M-file for tom_av2_ctf_topology_gui.fig
%      TOM_AV2_CTF_TOPOLOGY_GUI, by itself, creates a new TOM_AV2_CTF_TOPOLOGY_GUI or raises the existing
%      singleton*.
%
%      H = TOM_AV2_CTF_TOPOLOGY_GUI returns the handle to a new TOM_AV2_CTF_TOPOLOGY_GUI or the handle to
%      the existing singleton*.
%
%      TOM_AV2_CTF_TOPOLOGY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_AV2_CTF_TOPOLOGY_GUI.M with the given input arguments.
%
%      TOM_AV2_CTF_TOPOLOGY_GUI('Property','Value',...) creates a new TOM_AV2_CTF_TOPOLOGY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_av2_ctf_topology_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_av2_ctf_topology_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% save('Fit_structures.mat','Fit','Search','EM')
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_av2_ctf_topology_gui

% Last Modified by GUIDE v2.5 19-Apr-2010 12:10:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_av2_ctf_topology_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_ctf_topology_gui_OutputFcn, ...
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


% --- Executes just before tom_av2_ctf_topology_gui is made visible.
function tom_av2_ctf_topology_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_av2_ctf_topology_gui (see VARARGIN)

% Choose default command line output for tom_av2_ctf_topology_gui
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_av2_ctf_topology_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_av2_ctf_topology_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browse_align2d_file.
function browse_align2d_file_Callback(hObject, eventdata, handles)
% hObject    handle to browse_align2d_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f p]=uigetfile('.mat', 'Select the alignment file');
handles.align2d_filename = f;
handles.align2d_pathname = p;
%[a b c]=fileparts(f);
%handles.directoryname_extension=c;
set(handles.align2d_file,'String',[handles.align2d_pathname handles.align2d_filename]);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in browse_particle_stack.
function browse_particle_stack_Callback(hObject, eventdata, handles)
% hObject    handle to browse_particle_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f p]=uiputfile('.em', 'Create particle stack file');
handles.particle_filename = f;
handles.particle_pathname = p;
set(handles.particle_file,'String',[handles.particle_pathname handles.particle_filename]);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in fit_correct_button.
function fit_correct_button_Callback(hObject, eventdata, handles)
% hObject    handle to fit_correct_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

start_number=str2double(get(handles.start_number_gui,'String'));
disp(['Start number: ' num2str(start_number) ]);
ps_area=str2double(get(handles.ps_area_gui,'String'));
disp(['PS area: ' num2str(ps_area) ]);
ps_size=str2double(get(handles.ps_size_gui,'String'));
disp(['PS size: ' num2str(ps_size) ]);
step_size=str2double(get(handles.step_size_gui,'String'));
disp(['Step size: ' num2str(step_size) ]);

disp(['Load alignment file: ' handles.align2d_pathname '/' handles.align2d_filename]);
try
    load ([handles.align2d_pathname '/' handles.align2d_filename ]);
catch
    disp(['Cannot load alignment file']);    return;
end;
disp(['Correct in total: ' num2str(length(align2d)) ' particles.']);
disp(['Load CTF fit structure file: ' handles.ctf_fit_structure_pathname '/' handles.ctf_fit_structure_filename]);
try
    load ([handles.ctf_fit_structure_pathname '/' handles.ctf_fit_structure_filename ]);
    Fit_local=Fit;
    EM_local=EM;
    Search_local=Search;
catch
    disp(['Cannot load CTF fit structure file']);
end;

disp(['Create ctf parameter folder: ' handles.particle_pathname(1:end-1) '_' 'ctf_param']);
try
    warning off;
    mkdir([handles.particle_pathname(1:end-1) '_' 'ctf_param']);
    warning on;
catch
    disp(['Cannot create directory: ' handles.particle_pathname '/' 'ctf_param']);
    error;
end;
ctf_param_pathname=[handles.particle_pathname(1:end-1) '_' 'ctf_param/'];
last_filename='';
i=19006;
s=0;
idx=19006;
try
    fp=fopen([ctf_param_pathname 'ctfDatFile.dat'],'wt');
catch
    error('cannot open ctfDatfile.');
end;
try
    fp_sel=fopen([handles.particle_pathname '/parts.sel'],'wt');
catch
    error('cannot open parts.sel file.');
end;
while i<=length(align2d)
    
    if strcmp(last_filename,align2d(i).filename)==0
        disp(['-----------------------------------------------------------------------------']);
        disp(['Load image file: ' align2d(i).filename]);
        img=tom_emreadc(align2d(i).filename); img=img.Value;
        disp(['Load corresponding CTF fit file: ' align2d(i).filename '.mat']);
        load([align2d(i).filename '.mat']);
        Dz_range=500e-9;
        Dz_step=10e-9;
        disp(['Dz global: ' num2str(st_out.Fit.Dz_det) ' m.']);
        disp(['Search Dz Range: +/-' num2str(Dz_range) ' m.']);
        disp(['Search Dz step: ' num2str(Dz_step) ' m.']);
        Dz_start=(st_out.Fit.Dz_det-Dz_range);
        Dz_stop=(st_out.Fit.Dz_det+Dz_range);
        if Dz_stop>-1e-6; Dz_stop=-1e-6;Dz_start=-1.7e-6;end;        
        Search.Dz_search=[Dz_start:Dz_step:Dz_stop];
        disp(['Search Dz start: ' num2str(min(Search.Dz_search)) ' stop: ' num2str(max(Search.Dz_search))]);
        [Fit Dz_det Dz_delta_det statistic x y z]=tom_fit_ctf_topology(img,Fit,EM,Search,ps_area,ps_size,step_size,'plane');
        a=Fit.Plane.Hesse_coeff;
        disp(['new fit']);
        last_filename=align2d(i).filename;
        sz_img=size(img,1);
    end;
    
    x=align2d(i).position.x; % particle position on micrograph
    y=align2d(i).position.y;
    Dz_fit_particle=-((a(1) + a(2)*x+a(3)*y)./a(4)); % plane fitted defocus of particle on micrograph
    sz_particle=(align2d(i).radius).*4;
    pos_x_start=x-sz_particle./2+1;
    pos_x_stop=x+sz_particle./2;
    pos_y_start=y-sz_particle./2+1;
    pos_y_stop=y+sz_particle./2;
    if pos_x_start<1; pos_x_start=1; pos_x_stop=pos_x_start+sz_particle-1;end;
    if pos_y_start<1; pos_y_start=1; pos_y_stop=pos_y_start+sz_particle-1;end;
    if pos_x_stop>sz_img; pos_x_stop=sz_img; pos_x_start=pos_x_stop-sz_particle+1;end;
    if pos_y_stop>sz_img; pos_y_stop=sz_img; pos_y_start=pos_y_stop-sz_particle+1;end;
    
    particle_cut=img(pos_x_start:pos_x_stop,pos_y_start:pos_y_stop);
    sz_particle_cut=size(particle_cut,1);
    Fit=st_out.Fit;
    Fit.Dz_det=Dz_fit_particle;
    img_corr=tom_correct_for_ctf_and_mtf_new(double(particle_cut),Fit,'flip',sz_particle_cut./2-4,1); % corr_cut_off in pixel
    disp(['particle nr.: ' num2str(i) ' corrected with Dz: ' num2str(Dz_fit_particle)]);
    particle_filename=[handles.particle_pathname '/parts_' num2str(idx) '.spi'];
    sz_center_start=sz_particle_cut./2-align2d(i).radius+1;
    sz_center_stop=sz_particle_cut./2+align2d(i).radius;
    % center and normalize
    img_corr_center=img_corr(sz_center_start:sz_center_stop,sz_center_start:sz_center_stop);
    img_corr_center=tom_xmipp_normalize(img_corr_center,'Ramp');
    img_corr_center=tom_norm(-img_corr_center,'mean0+1std');
    try
        tom_spiderwrite(particle_filename,img_corr_center);
    catch
        error('cannot write particle file.');
    end;
    
    ctf_out_name=[ctf_param_pathname '/ctf_' num2str(idx) '.ctfparam'];
    write_ctfparam(ctf_out_name,Fit);
    fprintf(fp,[particle_filename ' ' ctf_out_name ' \n']);
    fprintf(fp_sel,[particle_filename ' 1\n']);

    idx=idx+1;
    i=i+1;
    %s=s+img_corr; tom_imagesc(s); drawnow;
    
end;
fclose(fp);
fclose(fp_sel);

function write_ctfparam(ctf_out_name,Fit)

%transfer parameters
sampling_rate=Fit.EM.Objectpixelsize.*1e10;
Voltage=Fit.EM.Voltage./1000;
spherical_aberration=Fit.EM.Cs.*1000;
chromatic_aberration=Fit.EM.Cc.*1000;

dz_u=Fit.Dz_det.*1e10+((-Fit.Dz_delta_det.*1e10)./2);
dz_v=Fit.Dz_det.*1e10-((-Fit.Dz_delta_det.*1e10)./2);
ang=Fit.Phi_0_det+135;

fp_ctf=fopen(ctf_out_name,'wt');
fprintf(fp_ctf,['sampling_rate=        '  num2str(sampling_rate) '\n']);
fprintf(fp_ctf,['voltage=              '  num2str(Voltage) '\n']);
fprintf(fp_ctf,['defocusU=             '  num2str(dz_u) '\n']);
fprintf(fp_ctf,['defocusV=             '  num2str(dz_v) '\n']);
fprintf(fp_ctf,['azimuthal_angle=      '  num2str(ang) '\n']);
fprintf(fp_ctf,['spherical_aberration= '  num2str(spherical_aberration) '\n']);
fprintf(fp_ctf,['chromatic_aberration= '  num2str(chromatic_aberration) '\n']);
fprintf(fp_ctf,['energy_loss=          '  num2str(0) '\n']);
fprintf(fp_ctf,['lens_stability=       '  num2str(0) '\n']);
fprintf(fp_ctf,['convergence_cone=     '  num2str(0) '\n']);
fprintf(fp_ctf,['longitudinal_displace='  num2str(0) '\n']);
fprintf(fp_ctf,['transversal_displace= '  num2str(0) '\n']);
fprintf(fp_ctf,['Q0=                   '  num2str(-0.1) '\n']);
fprintf(fp_ctf,['K=                    '  num2str(1) '\n']);
fclose(fp_ctf);


function start_number_gui_Callback(hObject, eventdata, handles)
% hObject    handle to start_number_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_number_gui as text
%        str2double(get(hObject,'String')) returns contents of start_number_gui as a double
%handles.start_number_gui=str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function start_number_gui_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_number_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ps_area_gui_Callback(hObject, eventdata, handles)
% hObject    handle to ps_area_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ps_area_gui as text
%        str2double(get(hObject,'String')) returns contents of ps_area_gui as a double


% --- Executes during object creation, after setting all properties.
function ps_area_gui_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ps_area_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ps_size_gui_Callback(hObject, eventdata, handles)
% hObject    handle to ps_size_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ps_size_gui as text
%        str2double(get(hObject,'String')) returns contents of ps_size_gui as a double


% --- Executes during object creation, after setting all properties.
function ps_size_gui_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ps_size_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function step_size_gui_Callback(hObject, eventdata, handles)
% hObject    handle to step_size_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step_size_gui as text
%        str2double(get(hObject,'String')) returns contents of step_size_gui as a double


% --- Executes during object creation, after setting all properties.
function step_size_gui_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_size_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_ctf_fit_structures.
function browse_ctf_fit_structures_Callback(hObject, eventdata, handles)
% hObject    handle to browse_ctf_fit_structures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f p]=uigetfile('.mat', 'Select CTF fit structure file');
handles.ctf_fit_structure_filename = f;
handles.ctf_fit_structure_pathname = p;
set(handles.CTF_fit_structures,'String',[f p]);
% Update handles structure
guidata(hObject, handles);
