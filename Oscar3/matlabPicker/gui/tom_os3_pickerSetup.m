function varargout = tom_os3_pickerSetup(varargin)
% TOM_OS3_PICKERSETUP M-file for tom_os3_pickerSetup.fig
%      TOM_OS3_PICKERSETUP, by itself, creates a new TOM_OS3_PICKERSETUP or raises the existing
%      singleton*.
%
%      H = TOM_OS3_PICKERSETUP returns the handle to a new TOM_OS3_PICKERSETUP or the handle to
%      the existing singleton*.
%
%      TOM_OS3_PICKERSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_OS3_PICKERSETUP.M with the given input arguments.
%
%      TOM_OS3_PICKERSETUP('Property','Value',...) creates a new TOM_OS3_PICKERSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_os3_pickerSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_os3_pickerSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_os3_pickerSetup

% Last Modified by GUIDE v2.5 06-Dec-2007 15:47:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_os3_pickerSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_os3_pickerSetup_OutputFcn, ...
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


% --- Executes just before tom_os3_pickerSetup is made visible.
function tom_os3_pickerSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_os3_pickerSetup (see VARARGIN)

% Choose default command line output for tom_os3_pickerSetup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_os3_pickerSetup wait for user response (see UIRESUME)
% uiwait(handles.figure1);

loadOptionsFile = false;

if(~loadOptionsFile)
%process X with Y store in Z ...
options.job.jobType = '2d';
options.job.volumeDirectory= '';
options.job.templateDirectory= '';
options.job.resultDirectory= '';
options.job.wisdomDir = '/fs/pool/pool-nickell/tmp/T2/fft_wisdom/matlab';
options.job.mode = 'parallel';
options.job.filefilter = 'none';

%correlation- what type
options.correlation.type = 'FLCF';
options.correlation.maskFile = 'none';
options.correlation.angles.start =0;
options.correlation.angles.end = 359;
options.correlation.angles.increment =10;

options.psf.file = 'none';

%where does the job run and whats its name?
options.parallel.jobManager ='cluster09';
options.parallel.jobName ='test';
options.parallel.nodeCount =1;
options.parallel.subVolumeSize =[];

%some binning or bandpass?
options.modifications.binning =0;
options.modifications.bandpass.low = -1;
options.modifications.bandpass.high = [];
options.modifications.bandpass.smoothing = [];

% the logical combination (forget about that, but must be set anyway)
options.analysis.ccc = 1;
options.analysis.psr = 1;
options.analysis.autocorr = 1;
options.analysis.pce = 0;
options.analysis.confidence = 0;

%how are the results stored?
options.result.saveList = true;
options.result.saveStack = true;
options.result.saveResults = true;
options.result.returnResults = false;

%what will the correlation filter look like?
options.filter.xcf = 0.5;
options.filter.psr = 0.5;
options.filter.soc = 0.5;
options.filter.numberParticles = 200;

%classification?
options.classification.enabled= true;
options.classification.trainingStack = '';
options.classification.eigenStart = 1;
options.classification.eigenEnd = 20;
options.classification.sigmaScale = [7 4 -0.75];
options.classification.numberOfClusters = 1;

%something else...
set(handles.rad_saveList,'Value',1);
set(handles.rad_saveStack,'Value',1);
set(handles.rad_saveFilter,'Value',1);

else
    
    
    
    
end;
[valid optionsValid numberInvalidOptions] = tom_os3_checkOptions(options);

%fct_setCheckSetup(handles,optionsValid,numberInvalidOptions,valid);

handles.picker.options = options;

guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = tom_os3_pickerSetup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_selVolDir.
function btn_selVolDir_Callback(hObject, eventdata, handles)
% hObject    handle to btn_selVolDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[DirName] = uigetdir('.','Select volume directory');

if(~ischar(DirName))
    return;
end;

options = handles.picker.options;
options.job.volumeDirectory = [DirName];
handles.picker.options = options;
guidata(hObject, handles);

set(handles.lbl_volDir,'String',DirName);

% --- Executes on button press in btn_selTempl.
function btn_selTempl_Callback(hObject, eventdata, handles)
% hObject    handle to btn_selTempl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[DirName] = uigetdir('.','Select template directory');

if(~ischar(DirName))
    return;
end;

options = handles.picker.options;
options.job.templateDirectory = [DirName];
handles.picker.options = options;
guidata(hObject, handles);

set(handles.lbl_tempDir,'String',DirName);

% --- Executes on button press in btn_selResDir.
function btn_selResDir_Callback(hObject, eventdata, handles)
% hObject    handle to btn_selResDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[DirName] = uigetdir('.','Select result directory');

if(~ischar(DirName))
    return;
end;

options = handles.picker.options;
options.job.resultDirectory = [DirName];
handles.picker.options = options;
guidata(hObject, handles);

set(handles.lbl_resDir,'String',DirName);

function JobName_Callback(hObject, eventdata, handles)
% hObject    handle to JobName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of JobName as text
%        str2double(get(hObject,'String')) returns contents of JobName as a double


% --- Executes during object creation, after setting all properties.
function JobName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to JobName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_corrMnu.
function pop_corrMnu_Callback(hObject, eventdata, handles)
% hObject    handle to pop_corrMnu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_corrMnu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_corrMnu
string = get(hObject,'String');
string= string{get(hObject,'Value')};
if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.correlation.type= string;
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pop_corrMnu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_corrMnu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_clustMnu.
function pop_clustMnu_Callback(hObject, eventdata, handles)
% hObject    handle to pop_clustMnu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_clustMnu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_clustMnu
string = get(hObject,'String');
if(iscell(string))
    string= string{get(hObject,'Value')};
end;

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.parallel.jobManager= string;
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pop_clustMnu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_clustMnu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_jobName_Callback(hObject, eventdata, handles)
% hObject    handle to txt_jobName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_jobName as text
%        str2double(get(hObject,'String')) returns contents of txt_jobName as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.parallel.jobName= string;
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_jobName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_jobName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_selFilFilt.
function btn_selFilFilt_Callback(hObject, eventdata, handles)
% hObject    handle to btn_selFilFilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName PathName] = uigetfile('.','Select volume directory');

if(~ischar(FileName))
    return;
end;

options = handles.picker.options;
options.job.filefilter = [PathName '/' FileName];
handles.picker.options = options;
guidata(hObject, handles);

set(handles.lbl_fileFilt,'String',FileName);


function txt_angStartPhi_Callback(hObject, eventdata, handles)
% hObject    handle to txt_angStartPhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_angStartPhi as text
%        str2double(get(hObject,'String')) returns contents of txt_angStartPhi as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.correlation.angles.start(1)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_angStartPhi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_angStartPhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in rad_2D.
function rad_2D_Callback(hObject, eventdata, handles)
% hObject    handle to rad_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rad_2D
Value2d= get(hObject,'Value');
Value2dOnl= get(handles.rad_2DOnline,'Value');
Value3d= get(handles.rad_3D,'Value');

if(Value2dOnl == 1 || Value3d == 1)
    warndlg('You are chaning the processing mode to 2D!');
end;

options = handles.picker.options;
options.job.jobType = '2d';
handles.picker.options = options;
guidata(hObject, handles);

set(handles.rad_2DOnline,'Value',0);
set(handles.rad_3D,'Value',0);

set(handles.txt_angStartPsi,'String','-');
set(handles.txt_angStartTheta,'String','-');
set(handles.txt_angIncPsi,'String','-');
set(handles.txt_angIncTheta,'String','-');
set(handles.txt_angEndPsi,'String','-');
set(handles.txt_angEndTheta,'String','-');

set(handles.txt_subVX,'String','1');
set(handles.txt_subVY,'String','-');
set(handles.txt_subVZ,'String','-');

% --- Executes on button press in rad_3D.
function rad_3D_Callback(hObject, eventdata, handles)
% hObject    handle to rad_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rad_3D
Value3d= get(hObject,'Value');
Value2dOnl= get(handles.rad_2DOnline,'Value');
Value2d= get(handles.rad_2D,'Value');

if(Value2dOnl == 1 || Value2d == 1)
    warndlg('You are chaning the processing mode to 3D!');
end;

options = handles.picker.options;
options.job.jobType = '3d';
handles.picker.options = options;
guidata(hObject, handles);

set(handles.rad_2D,'Value',0);
set(handles.rad_2DOnline,'Value',0);

set(handles.txt_angStartPsi,'String','');
set(handles.txt_angStartTheta,'String','');
set(handles.txt_angIncPsi,'String','');
set(handles.txt_angIncTheta,'String','');
set(handles.txt_angEndPsi,'String','');
set(handles.txt_angEndTheta,'String','');

set(handles.txt_subVX,'String','');
set(handles.txt_subVY,'String','');
set(handles.txt_subVZ,'String','');

% --- Executes on button press in btn_selPSF.
function btn_selPSF_Callback(hObject, eventdata, handles)
% hObject    handle to btn_selPSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName PathName] = uigetfile('.','Select psf file');

if(~ischar(FileName))
    return;
end;

options = handles.picker.options;
options.psf.file = [PathName '/' FileName];
handles.picker.options = options;
guidata(hObject, handles);

set(handles.lbl_psf,'String',FileName);

function txt_subVX_Callback(hObject, eventdata, handles)
% hObject    handle to txt_subVX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_subVX as text
%        str2double(get(hObject,'String')) returns contents of txt_subVX as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.parallel.subVolumeSize(1)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_subVX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_subVX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_subVY_Callback(hObject, eventdata, handles)
% hObject    handle to txt_subVY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_subVY as text
%        str2double(get(hObject,'String')) returns contents of txt_subVY as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.parallel.subVolumeSize(2)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_subVY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_subVY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_subVZ_Callback(hObject, eventdata, handles)
% hObject    handle to txt_subVZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_subVZ as text
%        str2double(get(hObject,'String')) returns contents of txt_subVZ as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.parallel.subVolumeSize(3)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_subVZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_subVZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_binning_Callback(hObject, eventdata, handles)
% hObject    handle to txt_binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_binning as text
%        str2double(get(hObject,'String')) returns contents of txt_binning as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.modifications.binning= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_binning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_selMask.
function btn_selMask_Callback(hObject, eventdata, handles)
% hObject    handle to btn_selMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName PathName] = uigetfile('.','Select mask file');

if(~ischar(FileName))
    return;
end;

options = handles.picker.options;
options.correlation.maskFile = [PathName FileName];
handles.picker.options = options;
guidata(hObject, handles);

set(handles.lbl_maskFile,'String',FileName);

% --- Executes on button press in btnStart.
function btnStart_Callback(hObject, eventdata, handles)
% hObject    handle to btnStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fct_setCheckSetup(handles,optionsValid,numberInvalid,valid)

    if(valid)
    else
        invalidSelections = get(handles,'String','invalidSelections');
    end;

tom_os3_picker(handles.picker.options);

% --- Executes on button press in btn_trainingSet.
function btn_trainingSet_Callback(hObject, eventdata, handles)
% hObject    handle to btn_trainingSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName PathName] = uigetfile('.','Select trainingset');

if(~ischar(FileName))
    return;
end;

options = handles.picker.options;
options.classification.trainingStack= [PathName  FileName];
handles.picker.options = options;
guidata(hObject, handles);

set(handles.lbl_trainingSet,'String',FileName);
msgbox('Note: the first image of the trainingset must be the reference used for aligning the image stack!');

function txt_eigenStart_Callback(hObject, eventdata, handles)
% hObject    handle to txt_eigenStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_eigenStart as text
%        str2double(get(hObject,'String')) returns contents of txt_eigenStart as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.classification.eigenStart= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_eigenStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_eigenStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_eigenEnd_Callback(hObject, eventdata, handles)
% hObject    handle to txt_eigenEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_eigenEnd as text
%        str2double(get(hObject,'String')) returns contents of txt_eigenEnd as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.classification.eigenEnd= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_eigenEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_eigenEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_filterXcf_Callback(hObject, eventdata, handles)
% hObject    handle to txt_filterXcf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_filterXcf as text
%        str2double(get(hObject,'String')) returns contents of txt_filterXcf as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.filter.xcf= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_filterXcf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_filterXcf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_filterPsr_Callback(hObject, eventdata, handles)
% hObject    handle to txt_filterPsr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_filterPsr as text
%        str2double(get(hObject,'String')) returns contents of txt_filterPsr as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.filter.psr= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_filterPsr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_filterPsr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_filterSoc_Callback(hObject, eventdata, handles)
% hObject    handle to txt_filterSoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_filterSoc as text
%        str2double(get(hObject,'String')) returns contents of txt_filterSoc as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.filter.soc= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_filterSoc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_filterSoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_sigmaStart_Callback(hObject, eventdata, handles)
% hObject    handle to txt_sigmaStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_sigmaStart as text
%        str2double(get(hObject,'String')) returns contents of txt_sigmaStart as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.classification.sigmaScale(1)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_sigmaStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_sigmaStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_sigmaEnd_Callback(hObject, eventdata, handles)
% hObject    handle to txt_sigmaEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_sigmaEnd as text
%        str2double(get(hObject,'String')) returns contents of txt_sigmaEnd as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.classification.sigmaScale(2)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_sigmaEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_sigmaEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_sigmaStep_Callback(hObject, eventdata, handles)
% hObject    handle to txt_sigmaStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_sigmaStep as text
%        str2double(get(hObject,'String')) returns contents of txt_sigmaStep as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.classification.sigmaScale(3)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_sigmaStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_sigmaStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_angStartPsi_Callback(hObject, eventdata, handles)
% hObject    handle to txt_angStartPsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_angStartPsi as text
%        str2double(get(hObject,'String')) returns contents of txt_angStartPsi as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.correlation.angles.start(2)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_angStartPsi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_angStartPsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_angStartTheta_Callback(hObject, eventdata, handles)
% hObject    handle to txt_angStartTheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_angStartTheta as text
%        str2double(get(hObject,'String')) returns contents of txt_angStartTheta as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.correlation.angles.start(3)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_angStartTheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_angStartTheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_angIncPhi_Callback(hObject, eventdata, handles)
% hObject    handle to txt_angIncPhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_angIncPhi as text
%        str2double(get(hObject,'String')) returns contents of txt_angIncPhi as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.correlation.angles.increment(1)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_angIncPhi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_angIncPhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_angIncPsi_Callback(hObject, eventdata, handles)
% hObject    handle to txt_angIncPsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_angIncPsi as text
%        str2double(get(hObject,'String')) returns contents of txt_angIncPsi as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.correlation.angles.increment(2)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_angIncPsi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_angIncPsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_angIncTheta_Callback(hObject, eventdata, handles)
% hObject    handle to txt_angIncTheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_angIncTheta as text
%        str2double(get(hObject,'String')) returns contents of txt_angIncTheta as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.correlation.angles.increment(3)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_angIncTheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_angIncTheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_angEndPhi_Callback(hObject, eventdata, handles)
% hObject    handle to txt_angEndPhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_angEndPhi as text
%        str2double(get(hObject,'String')) returns contents of txt_angEndPhi as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.correlation.angles.end(1)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_angEndPhi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_angEndPhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_angEndPsi_Callback(hObject, eventdata, handles)
% hObject    handle to txt_angEndPsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_angEndPsi as text
%        str2double(get(hObject,'String')) returns contents of txt_angEndPsi as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.correlation.angles.end(2)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_angEndPsi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_angEndPsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_angEndTheta_Callback(hObject, eventdata, handles)
% hObject    handle to txt_angEndTheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_angEndTheta as text
%        str2double(get(hObject,'String')) returns contents of txt_angEndTheta as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.correlation.angles.end(3)= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_angEndTheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_angEndTheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_checkOptions.
function btn_checkOptions_Callback(hObject, eventdata, handles)
% hObject    handle to btn_checkOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_file_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_load_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_save_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uiputfile('*','Save options file as...');

options = handles.picker.options;
save([PathName FileName], 'options');

% --- Executes on button press in rad_2DOnline.
function rad_2DOnline_Callback(hObject, eventdata, handles)
% hObject    handle to rad_2DOnline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rad_2DOnline
Value2dOnl= get(hObject,'Value');
Value2d= get(handles.rad_2D,'Value');
Value3d= get(handles.rad_3D,'Value');

if(Value3d == 1 || Value2d == 1)
    warndlg('You are chaning the processing mode to 2D-online!');
end;

options = handles.picker.options;
options.job.jobType = '2d-online';
handles.picker.options = options;
guidata(hObject, handles);

set(handles.rad_2D,'Value',0);
set(handles.rad_3D,'Value',0);

set(handles.txt_angStartPsi,'String','-');
set(handles.txt_angStartTheta,'String','-');
set(handles.txt_angIncPsi,'String','-');
set(handles.txt_angIncTheta,'String','-');
set(handles.txt_angEndPsi,'String','-');
set(handles.txt_angEndTheta,'String','-');

set(handles.txt_subVX,'String','1');
set(handles.txt_subVY,'String','-');
set(handles.txt_subVZ,'String','-');

function txt_filterNuP_Callback(hObject, eventdata, handles)
% hObject    handle to txt_filterNuP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_filterNuP as text
%        str2double(get(hObject,'String')) returns contents of txt_filterNuP as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.filter.numberParticles= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_filterNuP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_filterNuP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rad_clasEnable.
function rad_clasEnable_Callback(hObject, eventdata, handles)
% hObject    handle to rad_clasEnable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rad_clasEnable
Value= get(hObject,'Value');

options = handles.picker.options;
options.classification.enabled = Value == 1;
handles.picker.options = options;
guidata(hObject, handles);


% --- Executes on button press in btn_alignMask.
function btn_alignMask_Callback(hObject, eventdata, handles)
% hObject    handle to btn_alignMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName PathName] = uigetfile('.','Select a mask for PCA calculation.');

if(~ischar(FileName))
    return;
end;

options = handles.picker.options;
options.classification.alignMask= [PathName '/' FileName];
handles.picker.options = options;
guidata(hObject, handles);

set(handles.lbl_alMask,'String',FileName);


function txt_NuClust_Callback(hObject, eventdata, handles)
% hObject    handle to txt_NuClust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_NuClust as text
%        str2double(get(hObject,'String')) returns contents of txt_NuClust as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.picker.options;
options.classification.numberOfClusters= str2double(string);
handles.picker.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_NuClust_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_NuClust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rad_parallel.
function rad_parallel_Callback(hObject, eventdata, handles)
% hObject    handle to rad_parallel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rad_parallel
Value= get(hObject,'Value');

options = handles.picker.options;
if(Value == 1)
    options.job.mode = 'parallel';
else
    options.job.mode = 'sequential';
end;
handles.picker.options = options;
guidata(hObject, handles);


% --- Executes on button press in rad_saveFilter.
function rad_saveFilter_Callback(hObject, eventdata, handles)
% hObject    handle to rad_saveFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rad_saveFilter
value = get(hObject,'Value');
options = handles.picker.options;
if(value == 1)
    options.result.saveResults = true;
else
    options.result.saveResults = false;
end;
handles.picker.options = options;
guidata(hObject, handles);

msgbox('Set true if you want to use filter results for set up. Set the other options - save picklist & stack , classification to off!','Save options');

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_filterSetup_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_filterSetup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tom_os3_filterSetup


% --- Executes on button press in rad_saveList.
function rad_saveList_Callback(hObject, eventdata, handles)
% hObject    handle to rad_saveList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rad_saveList
value = get(hObject,'Value');
options = handles.picker.options;
if(value == 1)
    options.result.saveList = true;
else
    options.result.saveList = false;
end;
handles.picker.options = options;
guidata(hObject, handles);

msgbox('Save picklist either after filtering or - if classifier is on - after classification','Save options');


% --- Executes on button press in rad_saveStack.
function rad_saveStack_Callback(hObject, eventdata, handles)
% hObject    handle to rad_saveStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rad_saveStack
value = get(hObject,'Value');
options = handles.picker.options;
if(value == 1)
    options.result.saveStack = true;
else
    options.result.saveStack = false;
end;
handles.picker.options = options;
guidata(hObject, handles);

msgbox('Save particle stack either after filtering or - if classifier is on - after classification','Save options');
