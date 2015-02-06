function varargout = tom_av2_rec(varargin)
%TOM_AV2_REC creates ...
%
%   varargout = tom_av2_rec(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout           ...
%
%EXAMPLE
%   ... = tom_av2_rec(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_av2_rec_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_rec_OutputFcn, ...
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


%% --- Executes just before tom_av2_rec is made visible.
function tom_av2_rec_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_av2_rec (see VARARGIN)

% Choose default command line output for tom_av2_rec
handles.output = hObject;


%define allowed filters 

%PROJECTION
handles.storage.projection.filter_def.mask.types = {'sphere3d','cylinder3d'};

%CLASSIFICATION
handles.storage.classify.filter_def.mask.classify1.types = {'sphere','cylinder'};
handles.storage.classify.filter_def.mask.classify2.types = {'sphere','cylinder'};
handles.storage.classify.filter_def.mask.align.types = {'sphere','cylinder'};
handles.storage.classify.filter_def.mask.cc_rot.types = {'sphere','cylinder'};
handles.storage.classify.filter_def.mask.cc_trans.types = {'sphere','cylinder'};
handles.storage.classify.filter_def.filter.classify.types={'bandpass'};
handles.storage.classify.filter_def.filter.align.types={'bandpass'};


tmp.jobmanager='';
tmp.packageloss=0;
tmp.number_of_tasks=1;
tmp.workers.min=1;
tmp.workers.max=100;
tmp.timeout=0;
tmp.restart_workers=1;
handles.storage.classify.parallel=tmp;



%BACKPROJ
handles.storage.BackProj.filter_def.mask.types = {'sphere','cylinder'};
handles.storage.BackProj.filter_def.filter.types = {'bandpass'};


%POSTP
handles.storage.postp.filter_def.mask.types = {'sphere3d','cylinder3d'};


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_av2_rec wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%% --- Outputs from this function are returned to the command line.
function varargout = tom_av2_rec_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

%% ************************************************************************
%**************************FILE*****************************************
%% ************************************************************************


%% --- Button BROWSE ALIGNMENT FILE ---
function File_browse_alignment_Callback(hObject, eventdata, handles)

[filename,pathname] = uigetfile('*.mat','browse alignment file');

if (filename==0)
    set(handles.File_alignment,'String','');
else
    set(handles.File_alignment,'String',[pathname filename]);
end;


%% --- Button BROWSE STACK ---
function File_browse_stack_Callback(hObject, eventdata, handles)

[filename,pathname] = uigetfile('*.em, *.vol','browse particle stack');
tmpf=[pathname filename];
set(handles.File_stack,'String',tmpf);

h=tom_reademheader(tmpf);
sz=round(h.Header.Size(3)./2);

tmp=tom_emreadc([tmpf],'subregion',[1 1 sz],[h.Header.Size(1)-1 h.Header.Size(2)-1 0]);
display_thumbnail(handles.File_axes_stack,tmp.Value);
display_thumbnail(handles.classify_axes_particle,tmp.Value);
handles.storage.org_part=tmp.Value;

guidata(hObject, handles);

%% --- Button BROWSE START MODEL ---
function File_browse_startmodel_Callback(hObject, eventdata, handles)

st=get_gui_values(handles);
[filename,pathname] = uigetfile('*.em, *.vol','browse start model');
tmp=[pathname filename];
set(handles.File_startmodel,'String',tmp);
vol=tom_emreadc(tmp,'binning',st.File.Overall_binning);
vol=vol.Value;
sz=round(size(vol,3)./2);

handles.storage.start_proj=sum(vol(:,:,:),3);
handles.storage.start_central_slice=sum(vol(:,:,sz-5:sz+5),3);
handles.storage.start_model=vol;

% display thumbnails
display_thumbnail(handles.File_axes_startmodel,handles.storage.start_central_slice);
display_thumbnail(handles.project_axes_model,handles.storage.start_central_slice);
display_thumbnail(handles.BackProject_axes_Model,handles.storage.start_proj);
display_thumbnail(handles.Postp_axes,handles.storage.start_central_slice);
display_thumbnail(handles.classify_axes_model,handles.storage.start_proj);

%initialize work volums
handles.storage.BackProj.vol=vol;
handles.storage.postp.vol=vol;
handles.storage.projection.vol=vol;


guidata(hObject, handles);

%%  --- Executes on button press in File_browse_output_dir.
function File_browse_output_dir_Callback(hObject, eventdata, handles)

dirname = uigetdir;
if (dirname==0)
    dirname='';
end;
set(handles.File_outputdir,'String',dirname);



%% --- Executes on button press in File_Run_apply.
function File_Run_apply_Callback(hObject, eventdata, handles)
st=get_gui_values(handles);

if (st.File.Run.Step_nr > st.File.Run.Num_of_steps)
    errordlg('Act Step Number is bigger than Number of Steps');
    return;
end;

st.File.Run.Iterations_per_step{st.File.Run.Step_nr}=st.File.Run.Iterations;

sz_box=max(size(st.File.Run.Iterations_per_step));
sz=sz_box;
if (sz_box > st.File.Run.Num_of_steps)
    sz=st.File.Run.Num_of_steps;
end;

for i=1:sz
    tmp{i}=st.File.Run.Iterations_per_step{i};
end;


set(handles.File_Run_step_nr,'String',tmp);
set(handles.File_Run_step_nr,'Value',1);

guidata(hObject, handles);


%% --- Executes on button press in File_run_view_protocol.
function File_run_view_protocol_Callback(hObject, eventdata, handles)
% hObject    handle to File_run_view_protocol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
st=get_gui_values(handles);

error_m=check_inputs(handles,st,'param-view_protocol');

if (error_m==1)
    return;
end;

tom_av2_build_file_struct(st.File.Outputdir,'rec-logs');
    
tmp_path=[st.File.Outputdir '/logs/show_protocol.m'];  

%build m-file for publishing
fid=fopen(tmp_path,'w');
if (st.File.Run.Num_of_steps==1)
    fprintf(fid,'%s \n',['%% Overview: ' num2str(st.File.Run.Num_of_steps) ' Step is defined']);
else
    fprintf(fid,'%s \n',['%% Overview: ' num2str(st.File.Run.Num_of_steps) ' Steps are defined']);
end;
    
for i=1:st.File.Run.Num_of_steps
 
     tmp_nut=''; tmp_rot='';
     for ii=1:3
        tmp_nut=[tmp_nut ' ' st.project.angle.nutation_per_step{i,ii}];
        tmp_rot=[tmp_rot ' ' st.project.angle.rotation_per_step{i,ii}];
     end;
    fprintf(fid,'%s \n',['%% Step' num2str(i)]);
    fprintf(fid,'%s \n',['% Iteration: ',st.File.Run.Iterations_per_step{i}]);
    fprintf(fid,'%s \n',['%' ]);
    fprintf(fid,'%s \n',['% Nutation: ',tmp_nut]);
    fprintf(fid,'%s \n',['%' ]);
    fprintf(fid,'%s \n',['% Rotation: ',tmp_rot]);
    fprintf(fid,'%s \n',['%' ]);
    fprintf(fid,'%s \n',['% Classify Binning: ',st.classify.binning_all{i}]);
    fprintf(fid,'%s \n',['%' ]);
    fprintf(fid,'%s \n',['% Classify Alignment: ' st.classify.alignment_all{i}]);
end;
fclose(fid);

%publish tom html
options.format='html';
options.outputDir=[st.File.Outputdir '/logs/'];

publish(tmp_path,options);

[path name ext]=fileparts(tmp_path);
tmp_path=[path '/' name '.html'];

%open with matlab browser
open(tmp_path);



%%  --- Executes on button press in settings_load.
function settings_load_Callback(hObject, eventdata, handles)
% hObject    handle to settings_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[name path]=uigetfile;

st=load([path name]);
st=st.st;

if (isfield(st.project,'filter'))
    handles.storage.projection.filter=st.project.filter;
end;

if (isfield(st.classify,'filter'))
    handles.storage.classify.filter=st.classify.filter;
end;

if (isfield(st.Backproj,'filter'))
    handles.storage.BackProj.filter=st.Backproj.filter;
end;

if (isfield(st.Postp,'filter'))
    handles.storage.postp.filter=st.Postp.filter;
end;

set_gui_values(st,handles);


tmpf=st.File.Stack_Path;

h=tom_reademheader(tmpf);
sz=round(h.Header.Size(3)./2);

tmp=tom_emreadc(tmpf,'subregion',[1 1 sz],[h.Header.Size(1)-1 h.Header.Size(2)-1 0]);
display_thumbnail(handles.File_axes_stack,tmp.Value);
display_thumbnail(handles.classify_axes_particle,tmp.Value);
handles.storage.org_part=tmp.Value;

tmp=st.File.Startmodel_Path;
vol=tom_emreadc(tmp);
vol=vol.Value;
sz=round(size(vol,3)./2);

handles.storage.start_proj=sum(vol(:,:,:),3);
handles.storage.start_central_slice=sum(vol(:,:,sz-5:sz+5),3);
handles.storage.start_model=vol;

% display thumbnails
display_thumbnail(handles.File_axes_startmodel,handles.storage.start_central_slice);
display_thumbnail(handles.project_axes_model,handles.storage.start_central_slice);
display_thumbnail(handles.BackProject_axes_Model,handles.storage.start_proj);
display_thumbnail(handles.Postp_axes,handles.storage.start_central_slice);
display_thumbnail(handles.classify_axes_model,handles.storage.start_proj);

%initialize work volums
handles.storage.BackProj.vol=vol;
handles.storage.postp.vol=vol;
handles.storage.projection.vol=vol;

clear tmp;

for i=1:st.File.Run.Num_of_steps
    tmp{i}=num2str(i);
end;

%File
set(handles.File_run_step_show_step_nr,'Value',1);
set(handles.File_run_step_show_step_nr,'String',tmp);
%Project Run
set(handles.project_angles_step_nr,'Value',1);
set(handles.project_angles_step_nr,'String',tmp);
%Project Test
set(handles.project_test_step_nr,'Value',1);
set(handles.project_test_step_nr,'String',tmp);
%Classify Binning
set(handles.classify_step_nr_1,'Value',1);
set(handles.classify_step_nr_1,'String',tmp);
%Classify Alignment
set(handles.classify_step_nr,'Value',1);
set(handles.classify_step_nr,'String',tmp);




guidata(hObject, handles);

%% --- Executes on button press in settings_save.
function settings_save_Callback(hObject, eventdata, handles)
% hObject    handle to settings_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);

[name path]=uiputfile;

if (isfield(handles.storage.projection,'filter'))
    st.project.filter=handles.storage.projection.filter;
end;

if (isfield(handles.storage.classify,'filter'))
    st.classify.filter=handles.storage.classify.filter;
end;

if (isfield(handles.storage.BackProj,'filter'))
    st.Backproj.filter=handles.storage.BackProj.filter;
end;

if (isfield(handles.storage.postp,'filter'))
    st.Postp.filter=handles.storage.postp.filter;
end;

save([path name],'st');



%% --- Executes on button press in File_Go.
function File_Go_Callback(hObject, eventdata, handles)

st=get_gui_values(handles);
error_m=check_inputs(handles,st,['filter-projection param-projection filter-classify param-classify filter-BackProj  filter-postp param-postp']);

h=tom_reademheader(st.File.Startmodel_Path);
h2=tom_reademheader(st.File.Stack_Path);

if (error_m==1)
    return;
end;

%built alignment structure
if (isempty(st.File.Alignment_Path) )
    align2d=tom_av2_create_alignfromstack(st.File.Stack_Path);
else
    s=load(st.File.Alignment_Path);
    align2d=s.align2d;
end;

iter_num=size(align2d,1)+1;

%extend alignment structure by gui values

%FILE
align2d(iter_num,1).rec.file.Alignment_Path=st.File.Alignment_Path;
align2d(iter_num,1).rec.file.Stack_Path=st.File.Stack_Path;
align2d(iter_num,1).rec.file.Outputdir=st.File.Outputdir;
align2d(iter_num,1).rec.file.Startmodel_Path=st.File.Startmodel_Path;
align2d(iter_num,1).rec.file.Run=st.File.Run;
align2d(iter_num,1).rec.file.DemoMode=st.File.DemoMode;

%PROJECTION
nut=st.project.angle.nutation_per_step;
rot=st.project.angle.rotation_per_step;
for i=1:st.File.Run.Num_of_steps
        align2d(iter_num,1).rec.project.angle(i,:)=str2num([nut{i,1} ' ' nut{i,2} ' ' nut{i,3} ' ' rot{i,1} ' ' rot{i,2} ' ' rot{i,3}]);
end;
align2d(iter_num,1).rec.project.filter=handles.storage.projection.filter;
align2d(iter_num,1).rec.project.scheme=st.project.scheme;

%CLASSIFICATION
align2d(iter_num,1).rec.classify.filter=handles.storage.classify.filter;
align2d(iter_num,1).rec.classify.alignment=st.classify.alignment_all;
if (strcmp(st.File.Stack_num_of_parts,'all'))
    align2d(iter_num,1).rec.classify.particles_vector=eval(['1:' num2str(h2.Header.Size(3))]);
else
    align2d(iter_num,1).rec.classify.particles_vector=st.File.Stack_num_of_parts;
end;
align2d(iter_num,1).rec.classify.parallel=handles.storage.classify.parallel;
if (st.classify.paraell.on==0)
    align2d(iter_num,1).rec.classify.parallel.number_of_tasks=1;
end;
align2d(iter_num,1).rec.classify.indexs=st.classify.indexs;


%BACKPROJ
align2d(iter_num,1).rec.backproj.filter=handles.storage.BackProj.filter;
align2d(iter_num,1).rec.backproj.min_num_of_proj=st.Backproj.min_num_of_proj;
align2d(iter_num,1).rec.backproj.weighting=st.Backproj.weighting;

%POSTPROCESSING
align2d(iter_num,1).rec.post_proc.filter=handles.storage.postp.filter;
align2d(iter_num,1).rec.post_proc.sym.on=st.Postp.symmetry.on;
align2d(iter_num,1).rec.post_proc.sym.rot_angle=st.Postp.symmetry.rot_angle;
align2d(iter_num,1).rec.post_proc.sym.symmetry=st.Postp.symmetry.sym;
align2d(iter_num,1).rec.post_proc.binarize.on=st.Postp.binarize.on;
align2d(iter_num,1).rec.post_proc.binarize.mass=st.Postp.mass;
align2d(iter_num,1).rec.post_proc.binarize.pixelsize=st.Postp.pixelsize;   
align2d(iter_num,1).rec.post_proc.binarize.mass_is_black=st.Postp.mass_is_black;                                                   

%start the reconstruction process
tom_av2_reconstruction(align2d);


%% --- Textfield File_Run_num_of_steps---
function File_run_num_of_steps_Callback(hObject, eventdata, handles)
% hObject    handle to File_run_num_of_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of File_run_num_of_steps as text
%        str2double(get(hObject,'String')) returns contents of File_run_num_of_steps as a double

st=get_gui_values(handles);

%udate all step number popups
for i=1:st.File.Run.Num_of_steps
    tmp{i}=num2str(i);
end;

%File
set(handles.File_run_step_show_step_nr,'Value',1);
set(handles.File_run_step_show_step_nr,'String',tmp);
%Project Run
set(handles.project_angles_step_nr,'Value',1);
set(handles.project_angles_step_nr,'String',tmp);
%Project Test
set(handles.project_test_step_nr,'Value',1);
set(handles.project_test_step_nr,'String',tmp);
%Classify Binning
set(handles.classify_step_nr_1,'Value',1);
set(handles.classify_step_nr_1,'String',tmp);
%Classify Alignment
set(handles.classify_step_nr,'Value',1);
set(handles.classify_step_nr,'String',tmp);

guidata(hObject, handles);

%% --- Button VOLXYZ STACK ---
function File_volxyz_stack_Callback(hObject, eventdata, handles)

st=get_gui_values(handles);
im=tom_emreadc(st.File.Stack_Path);
tom_volxyz(im);

%% --- Button STACKBROWSER STACK ---
function File_stackbrowser_Callback(hObject, eventdata, handles)

st=get_gui_values(handles);
if (isempty(st.File.Alignment_Path))
    tom_av2_stackbrowser(st.File.Stack_Path);
else
    tom_av2_stackbrowser(st.File.Stack_Path,st.File.Alignment_Path);
end;


%% --- Button VOLXYZ START MODEL ---
function File_volxyz_startmodel_Callback(hObject, eventdata, handles)

st=get_gui_values(handles);
im=tom_emreadc(st.File.Startmodel_Path);
tom_volxyz(im);



%% --- Button DSPCUB START MODEL ---
function File_dspcub_startmodel_Callback(hObject, eventdata, handles)

st=get_gui_values(handles);
im=tom_emreadc(st.File.Startmodel_Path);
figure; tom_dspcub(im.Value);


%% ******************************************
% **************PROJECT**************
%% ****************************************


%%  --- Executes on button press in project_adjust_filter.
function project_adjust_filter_Callback(hObject, eventdata, handles)

in_struct=handles.storage.projection.filter_def;
st=get_gui_values(handles);
flag=check_inputs(handles,st,'filter-projection-verbose');

if (flag==0)
    handles.storage.projection.filter = tom_filtergui('mask',in_struct,handles.storage.projection.filter);
else
    handles.storage.projection.filter = tom_filtergui('mask',in_struct);

end;

guidata(hObject, handles);


%%  --- Executes on button press in project_apply_filter.
function project_apply_filter_Callback(hObject, eventdata, handles)
% hObject    handle to project_apply_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);
vol=handles.storage.start_model;
sz=round(size(vol,3)./2);

error_m=check_inputs(handles,st,'filter-projection');
if (error_m==1)
    return;
end;

mask=tom_create_mask(handles.storage.projection.filter.mask);

vol=tom_apply_mask(vol,mask,'mask and volume for Projection differ');
if (vol==0)
    return;
end;

%display image
display_thumbnail(handles.project_axes_model,sum(vol(:,:,sz-5:sz+5),3));

handles.storage.projection.vol=vol;

guidata(hObject, handles);

%% --- Executes on button press in project_test_project.
function project_test_project_Callback(hObject, eventdata, handles)
% hObject    handle to project_test_project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);


%check inputs
error_m=check_inputs(handles,st,'filter-projection param-projection file-warning');

tom_av2_build_file_struct(st.File.Outputdir,'rec',1);

if (error_m==1)
    return;
end;

iter_num=1;

% read file
vol=tom_emreadc(st.File.Startmodel_Path);
vol=vol.Value;

%FILE
align2d(iter_num,1).rec.file.Alignment_Path=st.File.Alignment_Path;
align2d(iter_num,1).rec.file.Stack_Path=st.File.Stack_Path;
align2d(iter_num,1).rec.file.Outputdir=st.File.Outputdir;
align2d(iter_num,1).rec.file.Startmodel_Path=st.File.Startmodel_Path;
align2d(iter_num,1).rec.file.Run=st.File.Run;
align2d(iter_num,1).rec.file.DemoMode=st.File.DemoMode;


%PROJECTION
nut=st.project.angle.nutation_per_step;
rot=st.project.angle.rotation_per_step;
for i=1:st.File.Run.Num_of_steps
        align2d(iter_num,1).rec.project.angle(i,:)=str2num([nut{i,1} ' ' nut{i,2} ' ' nut{i,3} ' ' rot{i,1} ' ' rot{i,2} ' ' rot{i,3}]);
end;
align2d(iter_num,1).rec.project.filter=handles.storage.projection.filter;
align2d(iter_num,1).rec.project.scheme=st.project.scheme;

count_st.hist=1; count_st.step=st.project.test.step_nr; count_st.iteration=1;
demo=st.File.DemoMode;

[align2d error_m]=tom_av2_create_projections(vol,align2d,count_st,demo);

%display only if demo mode is not running
if (demo.developer==0 & demo.presentation==0) 
    tmp=tom_emreadc([align2d(1,1).rec.project.path '/proj_' num2str(1) '.em']); 
    figure; tom_dspcub(tmp.Value);
end;

handles.storage.align2d=align2d;

guidata(hObject, handles);

%% --- Executes on button press in project_angle_apply.
function project_angle_apply_Callback(hObject, eventdata, handles)
% hObject    handle to project_angle_apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


st=get_gui_values(handles);

if (st.File.Run.Step_nr > st.File.Run.Num_of_steps)
    errordlg('Act Step Number is bigger than Number of Steps');
    return;
end;

error_m=check_inputs(handles,st,'param-Apply_Angle');

if (error_m==1)
    return;
end;

for i=1:3
    %exception handling
    st.project.angle.rotation_per_step{st.project.angle.step_nr,i}=str2num(st.project.angle.rotation{i});
    st.project.angle.nutation_per_step{st.project.angle.step_nr,i}=str2num(st.project.angle.nutation{i});
end;

sz_box=size(st.project.angle.rotation_per_step,1);
sz=sz_box;

if (sz_box > st.File.Run.Num_of_steps)
    sz=st.File.Run.Num_of_steps;
end;

for i=1:sz
    tmp_rot1{i}=st.project.angle.rotation_per_step{i,1};
    tmp_rot2{i}=st.project.angle.rotation_per_step{i,2};
    tmp_rot3{i}=st.project.angle.rotation_per_step{i,3};
    tmp_nut1{i}=st.project.angle.nutation_per_step{i,1};
    tmp_nut2{i}=st.project.angle.nutation_per_step{i,2};
    tmp_nut3{i}=st.project.angle.nutation_per_step{i,3};
end;

%set Values
set(handles.project_angles_rot_start_per_step,'Value',1);
set(handles.project_angles_rot_incre_per_step,'Value',1);
set(handles.project_angles_rot_stop_per_step,'Value',1);
set(handles.project_angles_nut_start_per_step,'Value',1);
set(handles.project_angles_nut_incre_per_step,'Value',1);
set(handles.project_angles_nut_stop_per_step,'Value',1);


%set strings
set(handles.project_angles_rot_start_per_step,'String',tmp_rot1);
set(handles.project_angles_rot_incre_per_step,'String',tmp_rot2);
set(handles.project_angles_rot_stop_per_step,'String',tmp_rot3);
set(handles.project_angles_nut_start_per_step,'String',tmp_nut1);
set(handles.project_angles_nut_incre_per_step,'String',tmp_nut2);
set(handles.project_angles_nut_stop_per_step,'String',tmp_nut3);

guidata(hObject, handles);


%% ---Textfield show iteration
function projection_angles_show_iteration_Callback(hObject, eventdata, handles)
st=get_gui_values(handles);

%nutation
set(handles.project_angles_nut_start_per_step,'Value',st.project.angle.show_step);
set(handles.project_angles_nut_incre_per_step,'Value',st.project.angle.show_step);
set(handles.project_angles_nut_stop_per_step,'Value',st.project.angle.show_step);

%rotation
set(handles.project_angles_rot_start_per_step,'Value',st.project.angle.show_step);
set(handles.project_angles_rot_incre_per_step,'Value',st.project.angle.show_step);
set(handles.project_angles_rot_stop_per_step,'Value',st.project.angle.show_step);





%% --- Button Project volxyz ---
function Project_volxyz_Callback(hObject, eventdata, handles)
tom_volxyz(handles.storage.projection.vol);


%% --- Button PROJECT_DSPCUB ---
function Project_dspcub_Callback(hObject, eventdata, handles)
figure; tom_dspcub(handles.storage.projection.vol);


%% **************************************************************************
%****CLASSIFY*************************************************************
%% *************************************************************************


%% --- Executes on button press in classify_adjust_filter.
function classify_adjust_filter_Callback(hObject, eventdata, handles)

st=get_gui_values(handles);
in_struct=handles.storage.classify.filter_def;
flag=check_inputs(handles,st,'filter-classify-verbose');

tmp.mask=in_struct.mask;
tmp2.filter=in_struct.filter;

if (flag==0)
    tmpp.mask=handles.storage.classify.filter.mask;
    tmpp2.filter=handles.storage.classify.filter.filter;
    handles.storage.classify.filter.mask = tom_filtergui('mask',tmp.mask,tmpp.mask);
    handles.storage.classify.filter.filter = tom_filtergui('filter',tmp2.filter,tmpp2.filter);
else
    tmp_out1=tom_filtergui('mask',tmp.mask);
    handles.storage.classify.filter.mask=tmp_out1;
    tmp_out2=tom_filtergui('filter',tmp2.filter);
    handles.storage.classify.filter.filter=tmp_out2; 
end;
%remove funny structure out of tom_filtergui


guidata(hObject, handles);


%% --- Executes on button press in classify_apply_filter.
function classify_apply_filter_Callback(hObject, eventdata, handles)

st=get_gui_values(handles);
flag=check_inputs(handles,st,'filter-classify');

if (flag==0)
    mask=tom_create_mask(handles.storage.classify.filter.mask.classify1);
    filt_st=handles.storage.classify.filter.filter.classify;
else
    errordlg('adjust filter first !'); return;
end;

im_tmp=handles.storage.start_proj;
im_tmp=tom_apply_filter(im_tmp,filt_st);
im_tmp=tom_apply_mask(im_tmp,mask,'mask size differs in classiy step');
if (im_tmp==0) 
    return;
end;
display_thumbnail(handles.classify_axes_model,im_tmp);

im_tmp=handles.storage.org_part;
im_tmp=tom_apply_filter(im_tmp,filt_st);
im_tmp=tom_apply_mask(im_tmp,mask,'mask size differs in classiy step');
if (im_tmp==0) 
    return;
end;
display_thumbnail(handles.classify_axes_particle,im_tmp);


%% --- Executes on button press in classify_correlate.
function classify_correlate_Callback(hObject, eventdata, handles)

st=get_gui_values(handles);

%process first image
im_tmp1=handles.storage.start_proj;

flag=check_inputs(handles,st,'filter-classify-verbose');
if (flag==0)
    im_tmp1=tom_apply_filter(im_tmp1,handles.storage.classify.filter.filter);
    im_tmp1=im_tmp1.*tom_create_mask(handles.storage.classify.filter.mask);    
end;

%process second image
im_tmp2=handles.storage.org_part;

flag=check_inputs(handles,st,'filter-classify-verbose');
if (flag==0)
    im_tmp2=tom_apply_filter(im_tmp2,handles.storage.classify.filter.filter);
    im_tmp2=im_tmp2.*tom_create_mask(handles.storage.classify.filter.mask);
end;

%correlate it
ccf=tom_corr(im_tmp1,im_tmp2,'norm');
display_thumbnail(handles.classify_axes_correlation,ccf);

%% --- Executes on button press in classify_test_classify.
function classify_test_classify_Callback(hObject, eventdata, handles)
% hObject    handle to classify_test_classify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles); 

%check for errors
error_m=check_inputs(handles,st,'filter-classify param-classify-align2d file-warning');

if (error_m==1)
    return;
end;

iter_num=1;
handles.storage.align2d(iter_num,1).rec.classify.filter=handles.storage.classify.filter;
handles.storage.align2d(iter_num,1).rec.classify.alignment=st.classify.alignment_all;
handles.storage.align2d(iter_num,1).rec.classify.particles_vector=st.classify.test.particle_vector;
handles.storage.align2d(iter_num,1).rec.classify.parallel=handles.storage.classify.parallel;
handles.storage.align2d(iter_num,1).rec.classify.indexs=st.classify.indexs;


if (st.classify.paraell.on==0)
    handles.storage.align2d(iter_num,1).rec.classify.parallel.number_of_tasks=1;
end;

count_st.hist=1; count_st.step=st.project.test.step_nr; count_st.iteration=1;
demo=st.File.DemoMode;

if (handles.storage.align2d(iter_num,1).rec.classify.indexs.on==0)
    [align2d error_m]=tom_av2_angular_classification(handles.storage.align2d,count_st,demo); % correlate projections and particles
else
    [align2d error_m]=tom_av2_index_angular_classification(handles.storage.align2d,count_st,demo); % correlate projections and particles
end;

if (error_m==0)
    tom_av2_plot_angleclasses(align2d,'show_all_classes',[1 1]);
else
    errordlg('classification error !!!');
    return;
end;

handles.storage.align2d=align2d;
guidata(hObject, handles);


%%  --- Executes on button press in classify_apply_binning.
function classify_apply_binning_Callback(hObject, eventdata, handles)
% hObject    handle to classify_apply_binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);

if (st.classify.binning_step > st.File.Run.Num_of_steps)
    errordlg('Act Step Number is bigger than Number of Steps');
    return;
end;

st.classify.binning_all{st.classify.binning_step}=num2str(st.classify.binning_tmp);

sz_box=max(size(st.classify.binning_all));
sz=sz_box;
if (sz_box > st.File.Run.Num_of_steps)
    sz=st.File.Run.Num_of_steps;
end;
for i=1:sz
    tmp{i}=st.classify.binning_all{i};
end;

set(handles.classify_binning_all,'String',tmp);
set(handles.classify_binning_all,'Value',1);


guidata(hObject, handles);


%% --- Executes on button press in classify_apply_alignment.
function classify_apply_alignment_Callback(hObject, eventdata, handles)
% hObject    handle to classify_apply_alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


st=get_gui_values(handles);

if (st.classify.binning_step > st.File.Run.Num_of_steps)
    errordlg('Act Step Number is bigger than Number of Steps');
    return;
end;

st.classify.alignment_all{st.classify.alignment_step}=st.classify.Alignment;
sz_box=max(size(st.classify.alignment_all));
sz=sz_box;
if (sz_box > st.File.Run.Num_of_steps)
    sz=st.File.Run.Num_of_steps;
end;

for i=1:sz
    tmp{i}=st.classify.alignment_all{i};
end;


set(handles.classify_alignment_all,'String',tmp);
set(handles.classify_alignment_all,'Value',1);

guidata(hObject, handles);



%% --- Executes on button press in classify_paralell_adjust.
function classify_paralell_adjust_Callback(hObject, eventdata, handles)

handles.storage.classify.parallel=tom_parallelsettings(handles.storage.classify.parallel);
guidata(hObject, handles);


%% **************************************************************************
% *******************BACKPROJ***********************************************
%% **************************************************************************


%% --- Executes on button press in BackProject_adjust_filter.
function BackProject_adjust_filter_Callback(hObject, eventdata, handles)


st=get_gui_values(handles);
in_struct=handles.storage.BackProj.filter_def;
flag=check_inputs(handles,st,'filter-BackProj-verbose');

tmp.mask=in_struct.mask;
tmp2.filter=in_struct.filter;

if (flag==0)
    tmpp.mask=handles.storage.BackProj.filter.mask;
    tmpp2.filter=handles.storage.BackProj.filter.filter;
    %handles.storage.BackProj.filter = tom_filtergui('mask',in_struct,handles.storage.BackProj.filter);
    handles.storage.BackProj.filter.mask = tom_filtergui('mask',tmp,tmpp.mask);
    handles.storage.BackProj.filter.filter = tom_filtergui('filter',tmp2,tmpp2.filter);
else
   handles.storage.BackProj.filter.mask = tom_filtergui('mask',tmp);
   handles.storage.BackProj.filter.filter = tom_filtergui('filter',tmp2);  
end

% tmp.mask=in_struct.mask;
% tmp2.filter=in_struct.filter;
% 
% if (flag==0)
%     tmpp.mask=handles.storage.classify.filter.mask;
%     tmpp2.filter=handles.storage.classify.filter.filter;
%     handles.storage.classify.filter.mask = tom_filtergui('mask',tmp.mask,tmpp.mask);
%     handles.storage.classify.filter.filter = tom_filtergui('filter',tmp2.filter,tmpp2.filter);
% else
%     tmp_out1=tom_filtergui('mask',tmp.mask);
%     handles.storage.classify.filter.mask=tmp_out1;
%     tmp_out2=tom_filtergui('filter',tmp2.filter);
%     handles.storage.classify.filter.filter=tmp_out2; 
% end;

guidata(hObject, handles);


%% --- Executes on button press in BackProject_apply_filter.
function BackProject_apply_filter_Callback(hObject, eventdata, handles)

st=get_gui_values(handles);
im_tmp=handles.storage.start_proj;

error_m=check_inputs(handles,st,'filter-BackProj');

if (error_m==1)
    return;
end;

im_tmp=tom_apply_filter(im_tmp,handles.storage.BackProj.filter.filter);
mask=tom_create_mask(handles.storage.BackProj.filter.mask);
im_tmp=tom_apply_mask(im_tmp,mask,'mask size differs from model size in backproj');

if (im_tmp==0)
    return;
end;
    
display_thumbnail(handles.BackProject_axes_Model,im_tmp);



%% --- Executes on button press in Backproject_test_backproj.
function Backproject_test_backproj_Callback(hObject, eventdata, handles)
% hObject    handle to Backproject_test_backproj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles); 
        

error_m=check_inputs(handles,st,'filter-BackProj param-Backproj file-warning');

if (error_m==1)
    return;
end;


iter_num=1;

%BACKPROJ
handles.storage.align2d(iter_num,1).rec.backproj.filter=handles.storage.BackProj.filter;
handles.storage.align2d(iter_num,1).rec.backproj.min_num_of_proj=st.Backproj.min_num_of_proj;
handles.storage.align2d(iter_num,1).rec.backproj.weighting=st.Backproj.weighting;

if (st.Backproj.test.use_startmod_proj)
    org_path=handles.storage.align2d(iter_num,1).rec.classify.avg.path;
    org_num_array=handles.storage.align2d(iter_num,1).rec.classify.avg.num_array;
    handles.storage.align2d(iter_num,1).rec.classify.avg.path=[handles.storage.align2d(iter_num,1).rec.project.path '/proj_1' handles.storage.align2d(iter_num,1).rec.project.ext];
    handles.storage.align2d(iter_num,1).rec.classify.avg.num_array=ones(max(size(org_num_array)),1)';
end;

count_st.hist=1; count_st.step=st.project.test.step_nr; count_st.iteration=1;
demo=st.File.DemoMode;

[vol error_m]=tom_av2_backproj(handles.storage.align2d,count_st,demo); 

sz=round(size(vol,3)./2);
axes(handles.BackProject_axes_Model);imagesc(sum(vol(:,:,sz-5:sz+5),3)'); colormap gray; set(gca,'XTick',[]); set(gca,'YTick',[]);

handles.storage.BackProj.vol=vol;

if (st.Backproj.test.use_startmod_proj)
    handles.storage.align2d(iter_num,1).rec.classify.avg.path=org_path;
    handles.storage.align2d(iter_num,1).rec.classify.avg.num_array=org_num_array;
end;

guidata(hObject, handles);



%% --- Executes on button press in BackProject_isosurf.
function BackProject_isosurf_Callback(hObject, eventdata, handles)
% hObject    handle to BackProject_isosurf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);
tom_isosurface(tom_norm(tom_bin(handles.storage.BackProj.vol,st.Backproj.iso_bin),1));


%% --- Executes on button press in Backproject_volxyz.
function Backproject_volxyz_Callback(hObject, eventdata, handles)
% hObject    handle to Backproject_volxyz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tom_volxyz(handles.storage.BackProj.vol); 

%% --- Executes on button press in Backproject_dspcub.
function Backproject_dspcub_Callback(hObject, eventdata, handles)
% hObject    handle to Backproject_dspcub (see GCBO)projection_angles_show_iteration
% eventdata  reserved - to be defined in a future version of MATLABf
% handles    structure with handles and user data (see GUIDATA)

figure; tom_dspcub(handles.storage.BackProj.vol); 



%% ***********************************************************************
%**************   POST PROCESSING   ************** 
%% *********************************************************************


%% --- Button VOLXYZ ---
function Postp_volxyz_Callback(hObject, eventdata, handles)

tom_volxyz(handles.storage.postp.vol);


%% --- Button DSPCUB ---
function Postp_dspcub_Callback(hObject, eventdata, handles)

figure; tom_dspcub(handles.storage.postp.vol);

%% --- Executes on button press in PostP_isosurface.
function PostP_isosurface_Callback(hObject, eventdata, handles)
% hObject    handle to PostP_isosurface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
st=get_gui_values(handles);

tom_isosurface(tom_norm(tom_bin(handles.storage.postp.vol,st.Postp.iso_binning),1));


%% --- Executes on button press in PostP_adjust_filter.
function PostP_adjust_filter_Callback(hObject, eventdata, handles)


st=get_gui_values(handles);
in_struct=handles.storage.postp.filter_def;
flag=check_inputs(handles,st,'filter-postp-verbose');

if (flag==0)
    
    
    handles.storage.postp.filter = tom_filtergui('mask',in_struct,handles.storage.postp.filter);

else


    handles.storage.postp.filter = tom_filtergui('mask',in_struct);
    handles.storage.postp.filter = tom_filtergui('mask',in_struct);
end;
 
guidata(hObject, handles);


%% --- Executes on button press in PostP_apply_filter.
function PostP_apply_filter_Callback(hObject, eventdata, handles)

st=get_gui_values(handles);
im_tmp=handles.storage.start_model;
sz=round(size(im_tmp,3)./2);
%projection_angles_show_iteration
error_m=check_inputs(handles,st,'filter-postp');

if (error_m==1)
    return;
end;

%im_tmp=tom_apply_filter(im_tmp,handles.storage.postp.filter.filter);
mask=tom_create_mask(handles.storage.postp.filter.mask);

im=tom_apply_mask(im_tmp,mask,'mask and model differ in size in postprocessing');

if (im==0)
    return;
end;

display_thumbnail(handles.Postp_axes,sum(im(:,:,sz-5:sz+5),3) );


guidata(hObject, handles);


%% --- Executes on button press in PostP_test_postp.
function PostP_test_postp_Callback(hObject, eventdata, handles)
% hObject    handle to PostP_test_postp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


st=get_gui_values(handles);

error_m=check_inputs(handles,st,'filter-postp param-postp file-warning');

if (error_m==1)
    return;
end;

%build file structure
tom_av2_build_file_struct(st.File.Outputdir,'rec',1);


%use start or tmp_model to test
if (st.Postp.test.use_startmod)
    vol=handles.storage.start_model;
else
    vol=handles.storage.BackProj.vol;
end;
sz=round(size(vol,3)./2);

iter_num=1;

%build up align structure to feed postp module
align2d(iter_num,1).rec.post_proc.filter=handles.storage.postp.filter;
align2d(iter_num,1).rec.post_proc.sym.on=st.Postp.symmetry.on;
align2d(iter_num,1).rec.post_proc.sym.rot_angle=st.Postp.symmetry.rot_angle;
align2d(iter_num,1).rec.post_proc.sym.symmetry=st.Postp.symmetry.sym;
align2d(iter_num,1).rec.post_proc.binarize.on=st.Postp.binarize.on;
align2d(iter_num,1).rec.post_proc.binarize.mass=st.Postp.mass;
align2d(iter_num,1).rec.post_proc.binarize.pixelsize=st.Postp.pixelsize;   
align2d(iter_num,1).rec.post_proc.binarize.mass_is_black=st.Postp.mass_is_black;                                                   
align2d(iter_num,1).rec.file.Outputdir=st.File.Outputdir;

if (isempty(st.File.Outputdir))
    errordlg('specify output dir before');
    return;
end;


count_st.hist=1; count_st.step=1; count_st.iteration=1;
[vol error_m]=tom_av2_post_processing(align2d,count_st,vol);

display_thumbnail(handles.Postp_axes,sum(vol(:,:,sz-5:sz+5),3)); 

handles.storage.postp.vol=vol;

guidata(hObject, handles);


%% *******************************
%*********HELPER FUNCTIONS*******
%% ********************************

%% ----get_gui_values------------------------------------------------------
function st=get_gui_values(handles)

%File
st.File.Alignment_Path=get(handles.File_alignment,'String');
st.File.Stack_Path=get(handles.File_stack,'String');
st.File.Startmodel_Path=get(handles.File_startmodel,'String');
st.File.Outputdir=get(handles.File_outputdir,'String');
st.File.Run.Num_of_steps=str2num(get(handles.File_run_num_of_steps,'String'));
tmp=get(handles.File_run_step_show_step_nr,'String'); 
num_tmp=get(handles.File_run_step_show_step_nr,'Value');
if (num_tmp > size(tmp,1))
    num_tmp=size(tmp,1);
end;

if (iscell(tmp));
    st.File.Run.Step_nr=str2num(tmp{num_tmp});
else
    st.File.Run.Step_nr=str2num(tmp);
end;

st.File.Run.Iterations=get(handles.File_run_num_of_iterations,'String');
st.File.Run.Iterations_per_step=get(handles.File_Run_step_nr,'String');
st.File.DemoMode.show_alignment=get(handles.File_demo_show_alignment,'Value');
st.File.DemoMode.presentation=get(handles.File_demo_presentation,'Value');
st.File.DemoMode.developer=get(handles.File_demo_developer,'Value');
tmp=get(handles.File_num_of_parts,'String');
if (isempty(findstr(tmp,':' )))
    st.File.Stack_num_of_parts=tmp;
else
    st.File.Stack_num_of_parts=eval(tmp);
end;
st.File.Overall_binning=str2num(get(handles.File_overall_binning,'String'));

%Project
tmp=get(handles.project_angles_step_nr,'String'); 
num_tmp=get(handles.project_angles_step_nr,'Value');
if (num_tmp > size(tmp,1))
    num_tmp=size(tmp,1);
end;
if (iscell(tmp))
    st.project.angle.step_nr=str2num(tmp{num_tmp});
else
    st.project.angle.step_nr=str2num(tmp);
end;

st.project.angle.show_step=str2num(get(handles.projection_angles_show_iteration,'String'));

%Values as String for further exception handling
st.project.angle.nutation{1}=get(handles.project_angles_nut_start,'String');
st.project.angle.nutation{2}=get(handles.project_angles_nut_incre,'String');
st.project.angle.nutation{3}=get(handles.project_angles_nut_stop,'String');
st.project.angle.rotation{1}=get(handles.project_angles_rot_start,'String');
st.project.angle.rotation{2}=get(handles.project_angles_rot_incre,'String');
st.project.angle.rotation{3}=get(handles.project_angles_rot_stop,'String');

tmp=get(handles.project_angles_scheme,'String');
st.project.scheme=tmp{get(handles.project_angles_scheme,'Value')};

tmp=get(handles.project_test_step_nr,'String'); 
num_tmp=get(handles.project_test_step_nr,'Value');
if (num_tmp > size(tmp,1))
    num_tmp=size(tmp,1);
end;

if (iscell(tmp))
    st.project.test.step_nr=str2num(tmp{num_tmp});
else
    st.project.test.step_nr=str2num(tmp);
end;


tmp=get(handles.project_angles_nut_start_per_step,'String');
for i=1:max(size(tmp));
    if (isempty(tmp))
        st.project.angle.nutation_per_step{i,1}='';
    else
        st.project.angle.nutation_per_step{i,1}=tmp{i};
    end;
end;
tmp=get(handles.project_angles_nut_incre_per_step,'String');
for i=1:max(size(tmp));
     if (isempty(tmp))
        st.project.angle.nutation_per_step{i,1}='';
     else
        st.project.angle.nutation_per_step{i,2}=tmp{i};
    end;
end;
tmp=get(handles.project_angles_nut_stop_per_step,'String');
for i=1:max(size(tmp))
    if (isempty(tmp))
        st.project.angle.nutation_per_step{i,3}='';
    else
       st.project.angle.nutation_per_step{i,3}=tmp{i}; 
    end;
        
end;
tmp=get(handles.project_angles_rot_start_per_step,'String');
for i=1:max(size(tmp))
    
    if (isempty(tmp))
        st.project.angle.rotation_per_step{i,3}='';
    else
        st.project.angle.rotation_per_step{i,1}=tmp{i};
    end;

end;
tmp=get(handles.project_angles_rot_incre_per_step,'String');
for i=1:max(size(tmp))
     if (isempty(tmp))
        st.project.angle.rotation_per_step{i,3}='';
     else
    
        st.project.angle.rotation_per_step{i,2}=tmp{i};
    end;
end;
tmp=get(handles.project_angles_rot_stop_per_step,'String');
for i=1:max(size(tmp))
    if (isempty(tmp))
        st.project.angle.rotation_per_step{i,3}='';
    else
        st.project.angle.rotation_per_step{i,3}=tmp{i};
    end;
end;

%Classify
st.classify.paraell.on=get(handles.classify_paralell_on,'Value');
tmp=get(handles.classify_alignment,'Value'); tmp_st=get(handles.classify_alignment,'String');
st.classify.Alignment=tmp_st{tmp};
st.classify.binning_tmp=str2num(get(handles.classify_binning,'String'));
tmp=get(handles.classify_step_nr_1,'String');
num_tmp=get(handles.classify_step_nr_1,'Value');
if (num_tmp > size(tmp,1))
    num_tmp=size(tmp,1);
end;
if (iscell(tmp))
    st.classify.binning_step=str2num(tmp{num_tmp});
else
    st.classify.binning_step=str2num(tmp);
end;

st.classify.binning_all=get(handles.classify_binning_all,'String');
tmp=get(handles.classify_step_nr,'String');
num_tmp=get(handles.classify_step_nr,'Value');
if (num_tmp > size(tmp,1))
    num_tmp=size(tmp,1);
end;
if (iscell(tmp))
    st.classify.alignment_step=str2num(tmp{num_tmp});
else
    st.classify.alignment_step=str2num(tmp);
end;

tmp=get(handles.classify_alignment,'String');
st.classify.alingment_tmp=tmp{get(handles.classify_alignment,'Value')};
st.classify.alignment_all=get(handles.classify_alignment_all,'String');
tmp=get(handles.classify_test_num_of_particles,'String');
if (isempty(findstr(tmp,':' )))
    st.classify.test.particle_vector=str2num(tmp);
else
    st.classify.test.particle_vector=eval(tmp);
end;

st.classify.indexs.on=get(handles.classify_indexs_on,'Value');
tmp=get(handles.classify_indexs_treesplit,'String');
st.classify.indexs.treesplit=tmp{get(handles.classify_indexs_treesplit,'Value')};
tmp=get(handles.classify_indexs_pca_eigs,'String');
[e1 r]=strtok(tmp,':'); e1=str2num(e1);
e2=strrep(r,':',''); e2=str2num(e2);
st.classify.indexs.pca_eigs=[e1 e2];
st.classify.indexs.pca_binning=str2num(get(handles.classify_indexs_pca_bin,'String'));


%Backproj
st.Backproj.min_num_of_proj=str2num(get(handles.BackProject_min_num_of_proj,'String'));
tmp=get(handles.BackProject_weighting,'Value'); tmp_st=get(handles.BackProject_weighting,'String');
st.Backproj.weighting=tmp_st{tmp};
st.Backproj.test.use_startmod_proj=get(handles.Backproject_test_use_start,'Value');
st.Backproj.test.use_tmpmod_poj=get(handles.Backproject_test_use_tmp,'Value');
st.Backproj.iso_bin=str2num(get(handles.BackProject_iso_bin,'String'));


%Postp
st.Postp.symmetry.on=get(handles.Postp_sym_on,'Value');
st.Postp.symmetry.rot_angle=str2num([get(handles.PostP_sym_rot_phi,'String') ' ' get(handles.PostP_sym_rot_psi,'String') ' ' get(handles.PostP_sym_rot_theta,'String')]);
st.Postp.symmetry.sym=str2num(get(handles.Postp_sym,'String'));
st.Postp.binarize.on=get(handles.Postp_bin_on,'Value');
st.Postp.mass=str2num(get(handles.Postp_mass,'String'));
st.Postp.mass_is_black=get(handles.Postp_mass_col,'Value');
st.Postp.test.use_startmod=get(handles.PostP_test_use_startmod,'Value');
st.Postp.test.use_tmp_mod=get(handles.PostP_test_use_tmp_model,'Value');
st.Postp.pixelsize=str2num(get(handles.PostP_pixelsize,'String'));
st.Postp.iso_binning=str2num(get(handles.postp_bin4iso,'String'));


%% set_gui_values
function set_gui_values(st,handles)
%File

set(handles.File_alignment,'String',st.File.Alignment_Path);
set(handles.File_stack,'String',st.File.Stack_Path)
set(handles.File_startmodel,'String',st.File.Startmodel_Path);
set(handles.File_outputdir,'String',st.File.Outputdir);
set(handles.File_run_num_of_steps,'String',st.File.Run.Num_of_steps);
set(handles.File_run_num_of_iterations,'String',st.File.Run.Iterations);
set(handles.File_Run_step_nr,'String',st.File.Run.Iterations_per_step);
set(handles.File_demo_show_alignment,'Value',st.File.DemoMode.show_alignment);
set(handles.File_demo_presentation,'Value',st.File.DemoMode.presentation);
set(handles.File_demo_developer,'Value',st.File.DemoMode.developer);
set(handles.File_num_of_parts,'String','all');
set(handles.File_overall_binning,'String',st.File.Overall_binning);

%Project
tmp=get(handles.project_angles_step_nr,'String'); 


num_tmp=get(handles.project_angles_step_nr,'Value');
if (num_tmp > size(tmp,1))
    num_tmp=size(tmp,1);
end;
if (iscell(tmp))
    st.project.angle.step_nr=str2num(tmp{num_tmp});
else
    st.project.angle.step_nr=str2num(tmp);
end;

st.project.angle.show_step=str2num(get(handles.projection_angles_show_iteration,'String'));


for i=1:size(st.project.angle.nutation_per_step,1)
    %nutation
    tmp_nut1{i}=st.project.angle.nutation_per_step{i,1};
    tmp_nut2{i}=st.project.angle.nutation_per_step{i,2};
    tmp_nut3{i}=st.project.angle.nutation_per_step{i,3};
    %rotation
    tmp_rot1{i}=st.project.angle.rotation_per_step{i,1};
    tmp_rot2{i}=st.project.angle.rotation_per_step{i,2};
    tmp_rot3{i}=st.project.angle.rotation_per_step{i,3};
end;

%nutation
set(handles.project_angles_nut_start_per_step,'String',tmp_nut1);
set(handles.project_angles_nut_incre_per_step,'String',tmp_nut2);
set(handles.project_angles_nut_stop_per_step,'String',tmp_nut3);
%rotation
set(handles.project_angles_rot_start_per_step,'String',tmp_rot1);
set(handles.project_angles_rot_incre_per_step,'String',tmp_rot2);
set(handles.project_angles_rot_stop_per_step,'String',tmp_rot3);

%Classify
set(handles.classify_paralell_on,'Value',st.classify.paraell.on);
set(handles.classify_binning_all,'String',st.classify.binning_all);
set(handles.classify_alignment_all,'String',st.classify.alignment_all);
 
%Backproj
set(handles.BackProject_min_num_of_proj,'String',num2str(st.Backproj.min_num_of_proj));
if (strcmp('Projection',st.Backproj.weighting))
    tmp=1;
else
    tmp=2;
end;
set(handles.BackProject_weighting,'Value',tmp);



%Postp
set(handles.Postp_sym_on,'Value',st.Postp.symmetry.on);
set(handles.Postp_sym,'Value',st.Postp.symmetry.sym);
set(handles.PostP_sym_rot_phi,'String',num2str(st.Postp.symmetry.rot_angle(1)));
set(handles.PostP_sym_rot_psi,'String',num2str(st.Postp.symmetry.rot_angle(2)));
set(handles.PostP_sym_rot_theta,'String',num2str(st.Postp.symmetry.rot_angle(3)));
set(handles.Postp_mass,'String',num2str(st.Postp.mass));
set(handles.PostP_pixelsize,'String',num2str(st.Postp.pixelsize));
if (st.Postp.mass_is_black==1)
    set(handles.Postp_mass_col,'Value',1);   
else
    set(handles.Postp_mass_col,'Value',2);
end;


%% ---check_inputs---------------------------------------------------
function error_m=check_inputs(handles,st,options)
%
%Performs exception handling for tom_av2_rec gui
%
%
%SYNTAX
%error_m=check_inputs(handles,st,options)
%
%DESCRIPTION
%Input
%   handles               :handles struct of gui
%   st                    :struct including gui values
%   options               :string with items to be checked                    
%Ouput
%   error_m               :error code 0 no errror 
%                                     1 error  
%
%EXAMPLE
%
% error_m=check_inputs(handles,st,'filter-projection param-projection file-warning');
%
%SEE ALSO
%
% get_gui_values
%
%Copyright (c) 2005
%TOM toolbox for Electron Tomography
%Max-Planck-Institute for Biochemistry
%Dept. Molecular Structural Biology
%82152 Martinsried, Germany
%http://www.biochem.mpg.de/tom
%
%Created 20/10/06 vom Isar Beder


error_m=0;

%check inputs

%% check Filter params
if(isempty(findstr(options,'filter'))==0)
    tmp=findstr(options,'filter');

    tmp_error=1;
    for i=1:max(size(tmp))
        task=strtok(options(tmp(i):end),' ');
        [t task]=strtok(task,'-');
        [task  verbose_flag]=strtok(task,'-');
        task=strrep(task,'-','');
        verbose_flag=strrep(verbose_flag,'-','');
         
        
        if (isfield(handles,'storage'))
            if (isfield(handles.storage,task))
                st_tmp=getfield(handles.storage,task);
                if (isfield(st_tmp,'filter'))
                    tmp_error=0;
                end;
            end;
        end;

        if (tmp_error==1)
            if (strcmp(verbose_flag,'verbose')==0)
                errordlg(['adjust filter for ' task]);
            end;
            error_m=1;
            break;
        end;
        tmp_error=1;

    end;
    
    if (error_m==1)
        return;
    end;
    
end;




%% check projection params
if (isempty(findstr(options,'param-projection'))==0 || isempty(findstr(options,'param-view_protocol'))==0 )
    if (isempty(st.File.Outputdir))
        errordlg('select a Output directory!');
        error_m=1;
        return;
    end;

    if (max(size(st.File.Run.Iterations_per_step)) < st.File.Run.Num_of_steps)
        errordlg('Iterations per step are not complete!');
        error_m=1;
        return
    end;

    if (size(st.project.angle.nutation_per_step,1) < st.File.Run.Num_of_steps)
        errordlg('Nutations per step are not complete!');
        error_m=1;
        return
    end;

    if (size(st.project.angle.rotation_per_step,1) < st.File.Run.Num_of_steps)
        errordlg('Rotations per step are not complete!');
        error_m=1;
        return
    end;
    
    nut=st.project.angle.nutation_per_step;
    rot=st.project.angle.rotation_per_step;
    for i=1:st.File.Run.Num_of_steps
        if (isempty(tom_create_angles_st(str2num([nut{i,1} ' ' nut{i,2} ' ' nut{i,3}]),str2num([rot{i,1} ' ' rot{i,2} ' ' rot{i,3}]),st.project.scheme))) 
            error_m=1;
            errordlg(['angular sampling fo step: ' num2str(i) ' corrupt']);
            return;
        end;
    end;
end;

%% check view protocol params
if (isempty(findstr(options,'param-view_protocol'))==0)
     if ( max(size(st.classify.alignment_all))  <  st.File.Run.Num_of_steps )
        errordlg('Alignment per step is not complete!');
        error_m=1;
        return
    end;
    
    
    tmp_flag=0;
    for i=1:max(size(st.classify.alignment_all))
        if (iscell((st.classify.alignment_all)))
            if (isempty(st.classify.alignment_all{i}))
                tmp_flag=1;
            end;
        else
            if (isempty(st.classify.alignment_all))
                tmp_flag=1;
            end;
        end;

    end;
    if (tmp_flag==1)
        errordlg('Alignment per step is not defined!');
        error_m=1;
        return;
    end;
    
    
    if ( max(size(st.classify.binning_all))  <  st.File.Run.Num_of_steps )
        errordlg('Classify Binning per step is not complete!');
        error_m=1;
        return
    end;
    
    
    tmp_flag=0;
    for i=1:max(size(st.classify.binning_all))
        if (iscell((st.classify.binning_all)))
            if (isempty(st.classify.binning_all{i}))
                tmp_flag=1;
            end;
        else
            if (isempty(st.classify.binning_all))
                tmp_flag=1;
            end;
        end;

    end;
    if (tmp_flag==1)
        errordlg('Classify Binning per step is not defined!');
        error_m=1;
        return;
    end;
    
end;


%% check classify params
if (isempty(findstr(options,'param-classify'))==0)

    
    if (findstr(options,'align2d'))

        if (isfield(handles.storage,'align2d')==0)
            errordlg('create projections first!');
            error_m=1;
            return;
        end;

    end;

    
    if ( max(size(st.classify.alignment_all))  <  st.File.Run.Num_of_steps )
        errordlg('Alignment per step is not complete!');
        error_m=1;
        return
    end;
    
    
    tmp_flag=0;
    for i=1:max(size(st.classify.alignment_all))
        if (iscell((st.classify.alignment_all)))
            if (isempty(st.classify.alignment_all{i}))
                tmp_flag=1;
            end;
        else
            if (isempty(st.classify.alignment_all))
                tmp_flag=1;
            end;
        end;

    end;
    if (tmp_flag==1)
        errordlg('Alignment per step is not defined!');
        error_m=1;
        return;
    end;
   
    
     if ( max(size(st.classify.binning_all))  <  st.File.Run.Num_of_steps )
        errordlg('Classify Binning per step is not complete!');
        error_m=1;
        return
    end;
    
    
    tmp_flag=0;
    for i=1:max(size(st.classify.binning_all))
        if (iscell((st.classify.binning_all)))
            if (isempty(st.classify.binning_all{i}))
                tmp_flag=1;
            end;
        else
            if (isempty(st.classify.binning_all))
                tmp_flag=1;
            end;
        end;

    end;
    if (tmp_flag==1)
        errordlg('Classify Binning per step is not defined!');
        error_m=1;
        return;
    end;
    
    
end;

 
 
%% check Backproj params
if (isempty(findstr(options,'param-Backproj'))==0)
    
    tmp_flag=0;
    if (isfield(handles.storage,'align2d'))
        if isfield(handles.storage.align2d(1,1).rec,'classify')
            tmp_flag=1;
        end;
    end;

    if (tmp_flag==0)
        errordlg('create class averages first');
        error_m=1;
        return;
    end;
    
end;



    
%% check postp params
if (isempty(findstr(options,'param-postp'))==0)
     if (isempty(st.File.Outputdir))
        errordlg('select a Output directory!');
        error_m=1;
        return;
    end;
end;


%% check apply angle params
if (isempty(findstr(options,'param-Apply_Angle'))==0)
    for i=1:3
        if (isempty(st.project.angle.rotation{i}) )
            errordlg('empty rotation!');
        end;
        if (isempty(st.project.angle.nutation{i}) )
            errordlg('empty nutation!');
        end;
    end;

end;

%% File Warning
if (isempty(findstr(options,'file-warning'))==0)
    button=questdlg(['Overwrite data in' st.File.Outputdir ' ?'], 'WARNING');
    
    if (strcmp(button,'Cancel') | strcmp(button,'No') )
        error_m=1; 
        return;
    end;
end;


%% display thumbnail
function display_thumbnail(axes_handle,image)

axes(axes_handle);
if (max(size(size(image)))==2)
    image=image';
end;

imagesc(image); colormap gray; set(gca,'XTick',[]); set(gca,'YTick',[]); 










%% **************   CALCULATION   ************** 

%--- Edit TASKS ---
function Calc_tasks_Callback(hObject, eventdata, handles)

%--- Edit FFTSIZE ---
function Calc_fftsize_Callback(hObject, eventdata, handles)


%**************   RECONSTRUCTION   ************** 

%--- Edit FILTER LOW ---
function Rec_filter_low_Callback(hObject, eventdata, handles)

%--- Edit FILTER HIGH ---
function Rec_filter_high_Callback(hObject, eventdata, handles)

%--- Radio Button WEIGHTING PROJECTION ---
function Rec_proj_Callback(hObject, eventdata, handles)

%--- Radio Button WEIGHTING VOLUME ---
function Rec_volume_Callback(hObject, eventdata, handles)


%**************   RUN PARAMETERS   **************



%--- Checkbox DISPLAY ---
function Runp_display_Callback(hObject, eventdata, handles)


%**************   DEGREES OF FREEDOM   **************

%--- Edit NUTATION START ---
function project_angles_nut_start_per_step_Callback(hObject, eventdata, handles)

%--- Edit NUTATION INCREMENT ---
function project_angles_nut_incre_per_step_Callback(hObject, eventdata, handles)

%--- Edit NUTATION STOP ---
function project_angles_nut_stop_per_step_Callback(hObject, eventdata, handles)

%--- Edit ROTATION START ---
function project_angles_rot_start_per_step_Callback(hObject, eventdata, handles)

%--- Edit ROTATION INCREMENT ---
function project_angles_rot_incre_per_step_Callback(hObject, eventdata, handles)

%--- Edit ROTATION STOP ---
function project_angles_rot_stop_per_step_Callback(hObject, eventdata, handles)

%--- Checkbox CORRECTION TRANSLATION ---
function Dof_translation_Callback(hObject, eventdata, handles)

%--- Checkbox CORRECTION ROTATION ---
function Dof_rotation_Callback(hObject, eventdata, handles)



%% --- Edit SYMMETRY ---
function Postp_sym_Callback(hObject, eventdata, handles)

%% --- Edit BINARIZE ---
function Postp_bin_Callback(hObject, eventdata, handles)

%% --- Edit MASS ---
function Postp_mass_Callback(hObject, eventdata, handles)

%% --- Edit MASK SIZE LENGHT ---
function Postp_masksize_length_Callback(hObject, eventdata, handles)

%% --- Edit MASK SIZE WIDTH ---
function Postp_masksize_high_Callback(hObject, eventdata, handles)

%% --- Edit SMOOTH ---
function Postp_smooth_Callback(hObject, eventdata, handles)


% --- Executes on selection change in classify_alignment.
function classify_alignment_Callback(hObject, eventdata, handles)
% hObject    handle to classify_alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns classify_alignment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from classify_alignment


% --- Executes during object creation, after setting all properties.
function classify_alignment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classify_alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end








% --- Executes during object creation, after setting all properties.
function File_run_num_of_steps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_run_num_of_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes during object creation, after setting all properties.
function File_Run_step_nr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_Run_step_nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function project_angles_step_nr_Callback(hObject, eventdata, handles)
% hObject    handle to project_angles_step_nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of project_angles_step_nr as text
%        str2double(get(hObject,'String')) returns contents of project_angles_step_nr as a double


% --- Executes during object creation, after setting all properties.
function project_angles_step_nr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to project_angles_step_nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function File_outputdir_Callback(hObject, eventdata, handles)
% hObject    handle to File_outputdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of File_outputdir as text
%        str2double(get(hObject,'String')) returns contents of File_outputdir as a double


% --- Executes during object creation, after setting all properties.
function File_outputdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_outputdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5




function BackProject_min_num_of_proj_Callback(hObject, eventdata, handles)
% hObject    handle to BackProject_min_num_of_proj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BackProject_min_num_of_proj as text
%        str2double(get(hObject,'String')) returns contents of BackProject_min_num_of_proj as a double


% --- Executes during object creation, after setting all properties.
function BackProject_min_num_of_proj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BackProject_min_num_of_proj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Postp_sym_on.
function Postp_sym_on_Callback(hObject, eventdata, handles)
% hObject    handle to Postp_sym_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Postp_sym_on


% --- Executes on button press in Postp_bin_on.
function Postp_bin_on_Callback(hObject, eventdata, handles)
% hObject    handle to Postp_bin_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Postp_bin_on




% --- Executes on button press in classify_paralell_on.
function classify_paralell_on_Callback(hObject, eventdata, handles)
% hObject    handle to classify_paralell_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of classify_paralell_on


% --- Executes on button press in File_demo_presentation.
function File_demo_presentation_Callback(hObject, eventdata, handles)
% hObject    handle to File_demo_presentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of File_demo_presentation


% --- Executes on button press in File_demo_developer.
function File_demo_developer_Callback(hObject, eventdata, handles)
% hObject    handle to File_demo_developer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of File_demo_developer


% --- Executes on button press in File_demo_show_alignment.
function File_demo_show_alignment_Callback(hObject, eventdata, handles)
% hObject    handle to File_demo_show_alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of File_demo_show_alignment




% --- Executes on button press in project_angles_view.
function project_angles_view_Callback(hObject, eventdata, handles)
% hObject    handle to project_angles_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function File_run_step_show_step_nr_Callback(hObject, eventdata, handles)
% hObject    handle to File_run_step_show_step_nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of File_run_step_show_step_nr as text
%        str2double(get(hObject,'String')) returns contents of File_run_step_show_step_nr as a double


% --- Executes during object creation, after setting all properties.
function File_run_step_show_step_nr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_run_step_show_step_nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function projection_angles_show_iteration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to projection_angles_show_iteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in project_angles_nut_start.
function project_angles_nut_start_Callback(hObject, eventdata, handles)
% hObject    handle to project_angles_nut_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns project_angles_nut_start contents as cell array
%        contents{get(hObject,'Value')} returns selected item from project_angles_nut_start


% --- Executes during object creation, after setting all properties.
function project_angles_nut_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to project_angles_nut_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in project_angles_nut_incre.
function project_angles_nut_incre_Callback(hObject, eventdata, handles)
% hObject    handle to project_angles_nut_incre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns project_angles_nut_incre contents as cell array
%        contents{get(hObject,'Value')} returns selected item from project_angles_nut_incre


% --- Executes during object creation, after setting all properties.
function project_angles_nut_incre_CreateFcn(hObject, eventdata, handles)
% hObject    handle to project_angles_nut_incre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in project_angles_nut_stop.
function project_angles_nut_stop_Callback(hObject, eventdata, handles)
% hObject    handle to project_angles_nut_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns project_angles_nut_stop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from project_angles_nut_stop


% --- Executes during object creation, after setting all properties.
function project_angles_nut_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to project_angles_nut_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in project_angles_rot_start.
function project_angles_rot_start_Callback(hObject, eventdata, handles)
% hObject    handle to project_angles_rot_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns project_angles_rot_start contents as cell array
%        contents{get(hObject,'Value')} returns selected item from project_angles_rot_start


% --- Executes during object creation, after setting all properties.
function project_angles_rot_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to project_angles_rot_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in project_angles_rot_incre.
function project_angles_rot_incre_Callback(hObject, eventdata, handles)
% hObject    handle to project_angles_rot_incre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns project_angles_rot_incre contents as cell array
%        contents{get(hObject,'Value')} returns selected item from project_angles_rot_incre


% --- Executes during object creation, after setting all properties.
function project_angles_rot_incre_CreateFcn(hObject, eventdata, handles)
% hObject    handle to project_angles_rot_incre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in project_angles_rot_stop.
function project_angles_rot_stop_Callback(hObject, eventdata, handles)
% hObject    handle to project_angles_rot_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns project_angles_rot_stop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from project_angles_rot_stop


% --- Executes during object creation, after setting all properties.
function project_angles_rot_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to project_angles_rot_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







% --- Executes on selection change in BackProject_weighting.
function BackProject_weighting_Callback(hObject, eventdata, handles)
% hObject    handle to BackProject_weighting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns BackProject_weighting contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BackProject_weighting


% --- Executes during object creation, after setting all properties.
function BackProject_weighting_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BackProject_weighting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on selection change in Postp_mass_col.
function Postp_mass_col_Callback(hObject, eventdata, handles)
% hObject    handle to Postp_mass_col (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Postp_mass_col contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Postp_mass_col


% --- Executes during object creation, after setting all properties.
function Postp_mass_col_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Postp_mass_col (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Edit RAD ---
function Prep_rad_Callback(hObject, eventdata, handles)

%--- Edit SIGMA ---
function Prep_sigma_Callback(hObject, eventdata, handles)



%--- Edit FILTER LOW ---
function Prep_filter_low_Callback(hObject, eventdata, handles)

%% --- Edit FILTER HIFH ---
function Prep_filter_high_Callback(hObject, eventdata, handles)



function PostP_sym_rot_phi_Callback(hObject, eventdata, handles)
% hObject    handle to PostP_sym_rot_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PostP_sym_rot_phi as text
%        str2double(get(hObject,'String')) returns contents of PostP_sym_rot_phi as a double


% --- Executes during object creation, after setting all properties.
function PostP_sym_rot_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PostP_sym_rot_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PostP_sym_rot_psi_Callback(hObject, eventdata, handles)
% hObject    handle to PostP_sym_rot_psi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PostP_sym_rot_psi as text
%        str2double(get(hObject,'String')) returns contents of PostP_sym_rot_psi as a double


% --- Executes during object creation, after setting all properties.
function PostP_sym_rot_psi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PostP_sym_rot_psi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PostP_sym_rot_theta_Callback(hObject, eventdata, handles)
% hObject    handle to PostP_sym_rot_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PostP_sym_rot_theta as text
%        str2double(get(hObject,'String')) returns contents of PostP_sym_rot_theta as a double


% --- Executes during object creation, after setting all properties.
function PostP_sym_rot_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PostP_sym_rot_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function PostP_pixelsize_Callback(hObject, eventdata, handles)
% hObject    handle to PostP_pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PostP_pixelsize as text
%        str2double(get(hObject,'String')) returns contents of PostP_pixelsize as a double


% --- Executes during object creation, after setting all properties.
function PostP_pixelsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PostP_pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in classify_alignment_all.
function classify_alignment_all_Callback(hObject, eventdata, handles)
% hObject    handle to classify_alignment_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns classify_alignment_all contents as cell array
%        contents{get(hObject,'Value')} returns selected item from classify_alignment_all


% --- Executes during object creation, after setting all properties.
function classify_alignment_all_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classify_alignment_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function classify_step_nr_Callback(hObject, eventdata, handles)
% hObject    handle to classify_step_nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of classify_step_nr as text
%        str2double(get(hObject,'String')) returns contents of classify_step_nr as a double


% --- Executes during object creation, after setting all properties.
function classify_step_nr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classify_step_nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in classify_test_stackbrowser.
function classify_test_stackbrowser_Callback(hObject, eventdata, handles)
% hObject    handle to classify_test_stackbrowser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function classify_test_num_of_particles_Callback(hObject, eventdata, handles)
% hObject    handle to classify_test_num_of_particles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of classify_test_num_of_particles as text
%        str2double(get(hObject,'String')) returns contents of classify_test_num_of_particles as a double


% --- Executes during object creation, after setting all properties.
function classify_test_num_of_particles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classify_test_num_of_particles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function project_test_step_nr_Callback(hObject, eventdata, handles)
% hObject    handle to project_test_step_nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of project_test_step_nr as text
%        str2double(get(hObject,'String')) returns contents of project_test_step_nr as a double


% --- Executes during object creation, after setting all properties.
function project_test_step_nr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to project_test_step_nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit45_Callback(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit45 as text
%        str2double(get(hObject,'String')) returns contents of edit45 as a double


% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function classify_test_part_nr_Callback(hObject, eventdata, handles)
% hObject    handle to classify_test_part_nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of classify_test_part_nr as text
%        str2double(get(hObject,'String')) returns contents of classify_test_part_nr as a double


% --- Executes during object creation, after setting all properties.
function classify_test_part_nr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classify_test_part_nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







% --- Executes on button press in PostP_test_use_startmod.
function PostP_test_use_startmod_Callback(hObject, eventdata, handles)
% hObject    handle to PostP_test_use_startmod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PostP_test_use_startmod


% --- Executes on button press in PostP_test_use_tmp_model.
function PostP_test_use_tmp_model_Callback(hObject, eventdata, handles)
% hObject    handle to PostP_test_use_tmp_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PostP_test_use_tmp_model





function classify_binning_Callback(hObject, eventdata, handles)
% hObject    handle to classify_binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of classify_binning as text
%        str2double(get(hObject,'String')) returns contents of classify_binning as a double


% --- Executes during object creation, after setting all properties.
function classify_binning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classify_binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on selection change in classify_binning_all.
function classify_binning_all_Callback(hObject, eventdata, handles)
% hObject    handle to classify_binning_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns classify_binning_all contents as cell array
%        contents{get(hObject,'Value')} returns selected item from classify_binning_all


% --- Executes during object creation, after setting all properties.
function classify_binning_all_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classify_binning_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function classify_step_nr_1_Callback(hObject, eventdata, handles)
% hObject    handle to classify_step_nr_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of classify_step_nr_1 as text
%        str2double(get(hObject,'String')) returns contents of classify_step_nr_1 as a double


% --- Executes during object creation, after setting all properties.
function classify_step_nr_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classify_step_nr_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in classify_test_advanced.
function classify_test_advanced_Callback(hObject, eventdata, handles)
% hObject    handle to classify_test_advanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function File_num_of_parts_Callback(hObject, eventdata, handles)
% hObject    handle to File_num_of_parts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of File_num_of_parts as text
%        str2double(get(hObject,'String')) returns contents of File_num_of_parts as a double


% --- Executes during object creation, after setting all properties.
function File_num_of_parts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_num_of_parts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function File_Run_step_nr_Callback(hObject, eventdata, handles)


function postp_bin4iso_Callback(hObject, eventdata, handles)
% hObject    handle to postp_bin4iso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of postp_bin4iso as text
%        str2double(get(hObject,'String')) returns contents of postp_bin4iso as a double


% --- Executes during object creation, after setting all properties.
function postp_bin4iso_CreateFcn(hObject, eventdata, handles)
% hObject    handle to postp_bin4iso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function File_overall_binning_Callback(hObject, eventdata, handles)
% hObject    handle to File_overall_binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of File_overall_binning as text
%        str2double(get(hObject,'String')) returns contents of File_overall_binning as a double


% --- Executes during object creation, after setting all properties.
function File_overall_binning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_overall_binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in Backproject_test_use_start.
function Backproject_test_use_start_Callback(hObject, eventdata, handles)
% hObject    handle to Backproject_test_use_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Backproject_test_use_start


% --- Executes on button press in Backproject_test_use_tmp.
function Backproject_test_use_tmp_Callback(hObject, eventdata, handles)
% hObject    handle to Backproject_test_use_tmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Backproject_test_use_tmp





function BackProject_iso_bin_Callback(hObject, eventdata, handles)
% hObject    handle to BackProject_iso_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BackProject_iso_bin as text
%        str2double(get(hObject,'String')) returns contents of BackProject_iso_bin as a double


% --- Executes during object creation, after setting all properties.
function BackProject_iso_bin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BackProject_iso_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%% --- Edit ITERATION ---
function File_run_num_of_iterations_Callback(hObject, eventdata, handles)



% --- Executes on button press in File_generate_startmodel.
function File_generate_startmodel_Callback(hObject, eventdata, handles)
% hObject    handle to File_generate_startmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in File_adjust_startmodel.
function File_adjust_startmodel_Callback(hObject, eventdata, handles)
% hObject    handle to File_adjust_startmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on selection change in BackProject_reconstruction_method.
function BackProject_reconstruction_method_Callback(hObject, eventdata, handles)
% hObject    handle to BackProject_reconstruction_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns BackProject_reconstruction_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BackProject_reconstruction_method


% --- Executes during object creation, after setting all properties.
function BackProject_reconstruction_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BackProject_reconstruction_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in project_angles_scheme.
function project_angles_scheme_Callback(hObject, eventdata, handles)
% hObject    handle to project_angles_scheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns project_angles_scheme contents as cell array
%        contents{get(hObject,'Value')} returns selected item from project_angles_scheme


% --- Executes during object creation, after setting all properties.
function project_angles_scheme_CreateFcn(hObject, eventdata, handles)
% hObject    handle to project_angles_scheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in classify_indexs_on.
function classify_indexs_on_Callback(hObject, eventdata, handles)
% hObject    handle to classify_indexs_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of classify_indexs_on


% --- Executes on selection change in classify_indexs_treesplit.
function classify_indexs_treesplit_Callback(hObject, eventdata, handles)
% hObject    handle to classify_indexs_treesplit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns classify_indexs_treesplit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from classify_indexs_treesplit


% --- Executes during object creation, after setting all properties.
function classify_indexs_treesplit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classify_indexs_treesplit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function classify_indexs_pca_eigs_Callback(hObject, eventdata, handles)
% hObject    handle to classify_indexs_pca_eigs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of classify_indexs_pca_eigs as text
%        str2double(get(hObject,'String')) returns contents of classify_indexs_pca_eigs as a double


% --- Executes during object creation, after setting all properties.
function classify_indexs_pca_eigs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classify_indexs_pca_eigs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function classify_indexs_pca_bin_Callback(hObject, eventdata, handles)
% hObject    handle to classify_indexs_pca_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of classify_indexs_pca_bin as text
%        str2double(get(hObject,'String')) returns contents of classify_indexs_pca_bin as a double


% --- Executes during object creation, after setting all properties.
function classify_indexs_pca_bin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classify_indexs_pca_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


