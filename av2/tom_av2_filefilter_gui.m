function varargout = tom_av2_filefilter_gui(varargin)
% TOM_AV2_FILEFILTER_GUI, graphical user interface for TOM_AV2_CTF_FIT_2_FILEFILTER
% which creates a 'particlepicker' structure with
% .filelist and .filefilter based on the fitted defocus (by tom_fit_ctf_gui).
%
%   [particlepicker in_range]=tom_av2_ctf_fit_2_filefilter(dirpath,percentage_off)
%
%PARAMETERS
%
%  INPUT
%   dirpath             path of mircorgraphs, .em and .em.mat files
%   percentage_off      difference to 'inteded defocus' from the micrograph
%                       header in percent
%  
%  OUTPUT
%   particlepicker      structure with .filelist and .filefilter
%   in_range            total number of micrographs in range of threshold
%
%EXAMPLE
%   %create 'particlepicker' based on fitted defocus, with a difference less
%   %than 30 percent of the intended defocus value:
%
%   [particlepicker in_range]=tom_av2_ctf_fit_2_filefilter('070410/07042010/low/low_3*.em',30); 
%
%REFERENCES
%
%  tom_fit_ctf_gui
%
%SEE ALSO
%   tom_av2_ctf_fit_2_filefilter
%
%   created by SN/FB 04/13/10
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

% Last Modified by GUIDE v2.5 14-Apr-2010 16:06:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_av2_filefilter_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_filefilter_gui_OutputFcn, ...
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


% --- Executes just before tom_av2_filefilter_gui is made visible.
function tom_av2_filefilter_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_av2_filefilter_gui (see VARARGIN)

% Choose default command line output for tom_av2_filefilter_gui
handles.output = hObject;
handles.defocus_percentage_off=20;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_av2_filefilter_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_av2_filefilter_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f p]=uigetfile('.em', 'Select a file in the directory');
handles.directoryname = p;
[a b c]=fileparts(f);
handles.directoryname_extension=c;
set(handles.dirpath,'String',[handles.directoryname '*' handles.directoryname_extension]);
% Update handles structure
guidata(hObject, handles);


function dirpath_Callback(hObject, eventdata, handles)
% hObject    handle to dirpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dirpath as text
%        str2double(get(hObject,'String')) returns contents of dirpath as a double


% --- Executes during object creation, after setting all properties.
function dirpath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dirpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in create.
function create_Callback(hObject, eventdata, handles)
% hObject    handle to create (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filter_file p] = uiputfile(pwd, 'Create Filterfile');
disp(['Create Filefilter: ' filter_file]);
disp(['with defocus difference of less than: ' num2str(handles.defocus_percentage_off) ' %']);
disp(['of folder: ' handles.directoryname '*' handles.directoryname_extension]);
[particlepicker in_range]=tom_av2_ctf_fit_2_filefilter([handles.directoryname '*' handles.directoryname_extension],handles.defocus_percentage_off);
for i=1:size(particlepicker.filefilter,2)
    particlepicker.filelist(i)=strrep(particlepicker.filelist(i),'low','high');
end;
save([p filter_file],'particlepicker');
disp(['Nr. of micrographs < threshold: ' num2str(in_range) ' of total nr.: ' num2str(size(particlepicker.filelist,2))]);
disp(['done.']);



function defocus_percentage_off_Callback(hObject, eventdata, handles)
% hObject    handle to defocus_percentage_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of defocus_percentage_off as text
%        str2double(get(hObject,'String')) returns contents of defocus_percentage_off as a double
handles.defocus_percentage_off=str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function defocus_percentage_off_CreateFcn(hObject, eventdata, handles)
% hObject    handle to defocus_percentage_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in low_to_high.
function low_to_high_Callback(hObject, eventdata, handles)
% hObject    handle to low_to_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of low_to_high

handles.low_to_high=get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);
