function varargout = tom_hpcl_status(varargin)
% TOM_HPCL_STATUS M-file for tom_hpcl_status.fig
%      TOM_HPCL_STATUS, by itself, creates a new TOM_HPCL_STATUS or raises the existing
%      singleton*.
%
%      H = TOM_HPCL_STATUS returns the handle to a new TOM_HPCL_STATUS or the handle to
%      the existing singleton*.
%
%      TOM_HPCL_STATUS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_HPCL_STATUS.M with the given input arguments.
%
%      TOM_HPCL_STATUS('Property','Value',...) creates a new TOM_HPCL_STATUS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_hpcl_status_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_hpcl_status_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_hpcl_status

% Last Modified by GUIDE v2.5 16-Sep-2010 11:01:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_hpcl_status_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_hpcl_status_OutputFcn, ...
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


% --- Executes just before tom_hpcl_status is made visible.
function tom_hpcl_status_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_hpcl_status (see VARARGIN)

% Choose default command line output for tom_hpcl_status
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_hpcl_status wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_hpcl_status_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in exit_button.
function exit_button_Callback(hObject, eventdata, handles)
% hObject    handle to exit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)


% --- Executes on button press in refresh_button.
function refresh_button_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[status hostname]=unix('hostname');

p=findstr(hostname,'hpcl');
if ~isempty(p)
    [status info]=unix('qstat');
    if status==1
        disp('qstat command not found.');return;
    end;
else
    [status info]=unix('ssh hpcl2001 qstat');
    if status==1
        disp('tried: ssh hpcl2001 qstat. Doesnt work.');return;
    end;
end;
p=findstr(info,'-');
info=info(p(end)+5:end);
idx=1;
while size(info,2)>100
    p=findstr(info,'.');    
    new_info(idx)=cellstr(info(p(1)-6:p(1)-6+99));
    info=info(p(1)-6+100:end);
    idx=idx+1;
end;
set(handles.listbox1,'String',new_info,...
	'Value',1)
handles.qstat_info=new_info;
guidata(hObject, handles);



% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
if strcmp(get(handles.figure1,'SelectionType'),'open')
    index_selected = get(handles.listbox1,'Value');
    job_info=handles.qstat_info{index_selected};
    job_ID=job_info(1:5);
    ButtonName = questdlg(['Delete job with ID: ' job_ID ' ?'], ...
                         ['Delete job ?'], ...
                         'Yes', 'No','No');
   switch ButtonName,
     case 'Yes',
      disp(['Delete job with ID: ' job_ID ' !']);
      [a b]=unix(['qdel ' job_ID]);
     case 'No',      
   end
end;

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
