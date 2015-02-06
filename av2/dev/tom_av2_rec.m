function varargout = tom_av2_rec(varargin)

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


% --- Executes just before tom_av2_rec is made visible.
function tom_av2_rec_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_av2_rec (see VARARGIN)

% Choose default command line output for tom_av2_rec
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_av2_rec wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_av2_rec_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

%**************   FILE   **************

%--- Button BROWSE ALIGNMENT FILE ---
function File_browse_alignment_Callback(hObject, eventdata, handles)

%--- Edit PATH & NAME ALIGNMENT FILE ---
function File_alignment_Callback(hObject, eventdata, handles)

%--- Button BROWSE STACK ---
function File_browse_stack_Callback(hObject, eventdata, handles)

%--- Edit PATH & NAME STACK ---
function File_stack_Callback(hObject, eventdata, handles)

%--- Button BROWSE START MODEL ---
function File_browse_startmodel_Callback(hObject, eventdata, handles)

%--- Edit PATH & NAME START MODEL ---
function File_startmodel_Callback(hObject, eventdata, handles)

%--- Button VOLXYZ STACK ---
function File_Prep_volxyz_stack_Callback(hObject, eventdata, handles)

%--- Button DSPCUB STACK ---
function File_dspcub_stack_Callback(hObject, eventdata, handles)

%--- Button VOLXYZ START MODEL ---
function File_volxyz_startmodel_Callback(hObject, eventdata, handles)

%--- Button DSPCUB START MODEL ---
function File_dspcub_startmodel_Callback(hObject, eventdata, handles)

%**************   PRE PROCESSING   **************   

%--- Button VOLXYZ MASK ---
function Prep_stack_volxyz_Callback(hObject, eventdata, handles)

%--- Button DSPCUB MASK ---
function Prep_stack_dspcub_Callback(hObject, eventdata, handles)

%--- Edit RAD ---
function Prep_rad_Callback(hObject, eventdata, handles)

%--- Edit SIGMA ---
function Prep_sigma_Callback(hObject, eventdata, handles)

%--- Button VOLXYZ  MODEL---
function Prep_model_volxyz_Callback(hObject, eventdata, handles)

%--- Button DSPCUB MODEL ---
function Prep_model_dspcub_Callback(hObject, eventdata, handles)

%--- Edit FILTER LOW ---
function Prep_filter_low_Callback(hObject, eventdata, handles)

%--- Edit FILTER HIFH ---
function Prep_filter_high_Callback(hObject, eventdata, handles)

%**************   POST PROCESSING   ************** 

%--- Edit SYMMETRY ---
function Postp_sym_Callback(hObject, eventdata, handles)

%--- Edit BINARIZE ---
function Postp_bin_Callback(hObject, eventdata, handles)

%--- Edit MASS ---
function Postp_mass_Callback(hObject, eventdata, handles)

%--- Edit MASK SIZE LENGHT ---
function Postp_masksize_length_Callback(hObject, eventdata, handles)

%--- Edit MASK SIZE WIDTH ---
function Postp_masksize_high_Callback(hObject, eventdata, handles)

%--- Edit SMOOTH ---
function Postp_smooth_Callback(hObject, eventdata, handles)

%--- Button VOLXYZ ---
function Postp_volxyz_Callback(hObject, eventdata, handles)

%--- Button DSPCUB ---
function Postp_dspcub_Callback(hObject, eventdata, handles)


%**************   CALCULATION   ************** 

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

%--- Edit ITERATION ---
function Runp_iteration_Callback(hObject, eventdata, handles)

%--- Checkbox DISPLAY ---
function Runp_display_Callback(hObject, eventdata, handles)


%**************   DEGREES OF FREEDOM   **************


%--- Edit NUTATION START ---
function Dof_nutation_start_Callback(hObject, eventdata, handles)

%--- Edit NUTATION INCREMENT ---
function Dof_nutation_inc_Callback(hObject, eventdata, handles)

%--- Edit NUTATION STOP ---
function Dof_nutation_stop_Callback(hObject, eventdata, handles)

%--- Edit ROTATION START ---
function Dof_rotation_start_Callback(hObject, eventdata, handles)

%--- Edit ROTATION INCREMENT ---
function Dof_rotation_inc_Callback(hObject, eventdata, handles)

%--- Edit ROTATION STOP ---
function Dof_rotation_stop_Callback(hObject, eventdata, handles)

%--- Checkbox CORRECTION TRANSLATION ---
function Dof_translation_Callback(hObject, eventdata, handles)

%--- Checkbox CORRECTION ROTATION ---
function Dof_rotation_Callback(hObject, eventdata, handles)




