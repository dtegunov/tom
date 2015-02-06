function varargout = tom_axes1_chooser(varargin)
%TOM_AXES1_CHOOSER by itself, creates a new TOM_AXES1_CHOOSER or raises the existing
%      singleton*.
%
%   varargout = tom_axes1_chooser(varargin)
%
%      H = TOM_AXES1_CHOOSER returns the handle to a new TOM_AXES1_CHOOSER or the handle to
%      the existing singleton*.
%
%      TOM_AXES1_CHOOSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_AXES1_CHOOSER.M with the given input arguments.
%
%      TOM_AXES1_CHOOSER('Property','Value',...) creates a new TOM_AXES1_CHOOSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_axes1_chooser_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_axes1_chooser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%EXAMPLE
%   ... = tom_axes1_chooser(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   GUIDE, GUIDATA, GUIHANDLES
%
% Copyright 2002-2003 The MathWorks, Inc.
%
% Edit the above text to modify the response to help tom_axes1_chooser
%
% Last Modified by GUIDE v2.5 09-Sep-2005 10:35:28


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_image_chooser_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_image_chooser_OutputFcn, ...
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


% --- Executes just before tom_axes1_chooser is made visible.
function tom_image_chooser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_axes1_chooser (see VARARGIN)

% Choose default command line output for tom_axes1_chooser
handles.output = hObject;
handles.two_images=get(handles.focal_pair,'Value');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_axes1_chooser wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_image_chooser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in yes.
function yes_Callback(hObject, eventdata, handles)
% hObject    handle to yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


 
disp('');


source_low_tmp=[handles.Source.path_low   handles.Source.basename  num2str(handles.Source.counter) handles.Source.ext];
source_high_tmp=[handles.Source.path_high     handles.Source.basename  num2str(handles.Source.counter) handles.Source.ext];
dest_low_tmp=[handles.Dest.path_low   handles.Dest.basename  num2str(handles.Dest.counter) handles.Dest.ext];
dest_high_tmp=[handles.Dest.path_high   handles.Dest.basename  num2str(handles.Dest.counter) handles.Dest.ext];


handles.cp_list{handles.Dest.counter,1}= source_low_tmp;
handles.cp_list{handles.Dest.counter,2} =dest_low_tmp;

handles.cp_list{handles.Dest.counter,3}=source_high_tmp;
handles.cp_list{handles.Dest.counter,4}=dest_high_tmp;


if  (handles.first == 1)
    handles.disp_high{handles.Dest.counter}=[handles.Dest.basename  num2str(handles.Dest.counter)];
    handles.disp_low{handles.Dest.counter}=[handles.Dest.basename  num2str(handles.Dest.counter)];
    handles.disp_high_us{handles.Used.counter}=[handles.Dest.basename  num2str(handles.Dest.counter)];
    handles.disp_low_us{handles.Used.counter}=[handles.Dest.basename  num2str(handles.Dest.counter)];
    handles.first=2;
else

    %insert new Values inverted  to have a better display in the list box

    help_high= handles.disp_high;
    help_low=handles.disp_low;

    handles.disp_high{handles.Dest.counter}='';
    handles.disp_low{handles.Dest.counter}='';

    for i=2:handles.Dest.counter
        handles.disp_high{i}=help_high{i-1};
        handles.disp_low{i}=help_low{i-1};
    end;

    handles.disp_high{1}=[handles.Dest.basename  num2str(handles.Dest.counter)];
    handles.disp_low{1}=[handles.Dest.basename  num2str(handles.Dest.counter)];

    help_high= handles.disp_high_us;
    help_low=handles.disp_low_us;

    handles.disp_high_us{handles.Used.counter}='';
    handles.disp_low_us{handles.Used.counter}='';

    for i=2:handles.Used.counter
        handles.disp_high_us{i}=help_high{i-1};
        handles.disp_low_us{i}=help_low{i-1};
    end;
    handles.disp_high_us{1}=[handles.Source.basename  num2str(handles.Source.counter)];
    handles.disp_low_us{1}=[handles.Source.basename  num2str(handles.Source.counter)];
end;



set(handles.File_Destination_low,'String',handles.disp_low ) ;
set(handles.File_Destination_High,'String',handles.disp_high );

set(handles.File_Used_low,'String',handles.disp_low_us ) ;
set(handles.File_Used_high,'String',handles.disp_high_us);



handles.Dest.counter=handles.Dest.counter+1;
handles.Used.counter=handles.Used.counter+1;
handles.Source.counter=handles.Source.counter+1;



set(handles.File_Source_High,'String',[handles.Source.basename  num2str(handles.Source.counter) handles.Source.ext  ]);
set(handles.File_Source_low,'String',[handles.Source.basename  num2str(handles.Source.counter) handles.Source.ext ]);
set(handles.Destination_counter,'String',num2str(handles.Dest.counter));


source_high_tmp=[handles.Source.path_high     handles.Source.basename  num2str(handles.Source.counter) handles.Source.ext];

try
    im=tom_emread(source_high_tmp);
catch
    msgbox(['image ' source_high_tmp ' missing']);
    handles.Source.counter=handles.Source.counter+1;
    guidata(hObject, handles);
    return;
end;

im=tom_bin(im.Value,3);
axes(handles.axes1);
tom_imagesc(im);
h=handles.cp_list;

save('cp_list','h');

guidata(hObject, handles);


% --- Executes on selection change in File_Destination_low.
function File_Destination_low_Callback(hObject, eventdata, handles)
% hObject    handle to File_Destination_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns File_Destination_low contents as cell array
%        contents{get(hObject,'Value')} returns selected item from File_Destination_low


% --- Executes during object creation, after setting all properties.
function File_Destination_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_Destination_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function File_Source_low_Callback(hObject, eventdata, handles)
% hObject    handle to File_Source_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of File_Source_low as text
%        str2double(get(hObject,'String')) returns contents of File_Source_low as a double


% --- Executes during object creation, after setting all properties.
function File_Source_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_Source_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Browse_Source.
function Browse_Source_Callback(hObject, eventdata, handles)
% hObject    handle to Browse_Source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox('select the first image in high!');
uiwait;
[myname, mypathname_high] = uigetfile('*.*')
myfile=[mypathname_high myname];
if (isunix==1)
    mypathname_low=[mypathname_high(1:findstr(mypathname_high,'high')-1) 'low/' ];
else
    mypathname_low=[mypathname_high(1:findstr(mypathname_high,'high')-1) 'low\' ];
end;
    

set(handles.Path_Source_low,'String',mypathname_low);
set(handles.File_Source_low,'String',myname);
set(handles.Path_Source_High,'String',mypathname_high);
set(handles.File_Source_High,'String',myname);


% --- Executes on button press in Browse_Destination.
function Browse_Destination_Callback(hObject, eventdata, handles)
% hObject    handle to Browse_Destination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('select the folder high to put the choosen images!');
uiwait;
[myname, mypathname_high] = uiputfile;
myfile=[mypathname_high myname];
if (isunix==1)
    mypathname_low=[mypathname_high(1:findstr(mypathname_high,'high')-1) 'low/' ];
else
    mypathname_low=[mypathname_high(1:findstr(mypathname_high,'high')-1) 'low\' ];
end;

set(handles.Path_Destination_low,'String',mypathname_low);
set(handles.Path_Destination_High,'String',mypathname_high);
[a,count]=strtok(myname,'_');
count=num2str(strrep(count,'_',''));
[count,ext]=strtok(count,'.');
extension=[ext];
set(handles.Destination_Ending,'String',extension);
set(handles.Destination_Filename,'String',[strtok(myname,'_')  '_']);
set(handles.Destination_counter,'String',count);
% set(handles.File_Destination_High,'String',myname);
% set(handles.File_Destination_low,'String',myname);

function Path_Destination_low_Callback(hObject, eventdata, handles)
% hObject    handle to Path_Destination_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Path_Destination_low as text
%        str2double(get(hObject,'String')) returns contents of Path_Destination_low as a double


% --- Executes during object creation, after setting all properties.
function Path_Destination_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Path_Destination_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Start.
function Start_Callback(hObject, eventdata, handles)
% hObject    handle to Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path=get(handles.Path_Source_High,'String');
name=get(handles.File_Source_High,'String');


[basename rest]=strtok(name,'_');
basename=[basename '_'];
rest=strrep(rest,'_','');
[count ext]=strtok(rest,'.');
count=str2num(count);
path_name=[path basename num2str(count) ext];

im=tom_emread(path_name);
im=tom_bin(im.Value,3);
axes(handles.axes1);
tom_imagesc(im);


handles.first=1;
set(handles.File_Source_low,'String',[basename num2str(count) ext]);

%initialize Variables Source
handles.Source.counter=count;
handles.Source.ext=ext;
handles.Source.basename=basename;
handles.Source.path_low=get(handles.Path_Source_low,'String');
handles.Source.path_high=get(handles.Path_Source_High,'String');


%initialize Variables Destination
handles.Dest.counter=str2num(get(handles.Destination_counter,'String'));
handles.Dest.ext=get(handles.Destination_Ending,'String');
handles.Dest.basename=get(handles.Destination_Filename,'String');
handles.Dest.path_low=get(handles.Path_Destination_low,'String');
handles.Dest.path_high=get(handles.Path_Destination_High,'String');


handles.Thrash.counter=1;
handles.Used.counter=1;

 guidata(hObject, handles);
    


function path_source_Callback(hObject, eventdata, handles)
% hObject    handle to path_source_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path_source_low as text
%        str2double(get(hObject,'String')) returns contents of path_source_low as a double


% --- Executes during object creation, after setting all properties.
function path_source_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path_source_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in no.
function no_Callback(hObject, eventdata, handles)
% hObject    handle to no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if  (handles.Thrash.counter == 1)
    handles.disp_high_th{handles.Thrash.counter}=[handles.Dest.basename  num2str(handles.Source.counter)];
    handles.disp_low_th{handles.Thrash.counter}=[handles.Dest.basename  num2str(handles.Source.counter)];
else
    
    %insert new Values inverted  to have a better display in the list box
    
    help_high= handles.disp_high_th;
    help_low=handles.disp_low_th;
    
    handles.disp_high_th{handles.Thrash.counter}='';
    handles.disp_low_th{handles.Thrash.counter}='';
    
    for i=2:handles.Thrash.counter
         handles.disp_high_th{i}=help_high{i-1};
         handles.disp_low_th{i}=help_low{i-1};
    end;
   
    
    handles.disp_high_th{1}=[handles.Source.basename  num2str(handles.Source.counter)];
    handles.disp_low_th{1}=[handles.Source.basename  num2str(handles.Source.counter)];
    
end;


handles.Source.counter=handles.Source.counter+1;   
handles.Thrash.counter=handles.Thrash.counter+1;

source_high_tmp=[handles.Source.path_high     handles.Source.basename  num2str(handles.Source.counter) handles.Source.ext];

set(handles.File_Source_High,'String',[handles.Source.basename  num2str(handles.Source.counter)]);
set(handles.File_Source_low,'String',[handles.Source.basename  num2str(handles.Source.counter)]);

set(handles.File_Thrash_low,'String',handles.disp_low_th ) ;
set(handles.File_Thrash_high,'String',handles.disp_high_th );


try
    im=tom_emread(source_high_tmp);
catch
    msgbox('last image clicked '),
end;

im=tom_bin(im.Value,3);
axes(handles.axes1);
tom_imagesc(im);

guidata(hObject, handles);
    

function Path_Source_low_Callback(hObject, eventdata, handles)
% hObject    handle to Path_Source_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Path_Source_low as text
%        str2double(get(hObject,'String')) returns contents of Path_Source_low as a double


% --- Executes during object creation, after setting all properties.
function Path_Source_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Path_Source_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function Path_Source_high_Callback(hObject, eventdata, handles)
% hObject    handle to Path_Source_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Path_Source_high as text
%        str2double(get(hObject,'String')) returns contents of Path_Source_high as a double


% --- Executes during object creation, after setting all properties.
function Path_Source_high_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Path_Source_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Path_Destination_High_Callback(hObject, eventdata, handles)
% hObject    handle to Path_Destination_High (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Path_Destination_High as text
%        str2double(get(hObject,'String')) returns contents of Path_Destination_High as a double


% --- Executes during object creation, after setting all properties.
function Path_Destination_High_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Path_Destination_High (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in File_Destination_High.
function File_Destination_High_Callback(hObject, eventdata, handles)
% hObject    handle to File_Destination_High (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns File_Destination_High contents as cell array
%        contents{get(hObject,'Value')} returns selected item from File_Destination_High


% --- Executes during object creation, after setting all properties.
function File_Destination_High_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_Destination_High (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function File_Source_High_Callback(hObject, eventdata, handles)
% hObject    handle to File_Source_High (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of File_Source_High as text
%        str2double(get(hObject,'String')) returns contents of File_Source_High as a double


% --- Executes during object creation, after setting all properties.
function File_Source_High_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_Source_High (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Path_Source_High_Callback(hObject, eventdata, handles)
% hObject    handle to Path_Source_High (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Path_Source_High as text
%        str2double(get(hObject,'String')) returns contents of Path_Source_High as a double


% --- Executes during object creation, after setting all properties.
function Path_Source_High_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Path_Source_High (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function Destination_Filename_Callback(hObject, eventdata, handles)
% hObject    handle to Destination_Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Destination_Filename as text
%        str2double(get(hObject,'String')) returns contents of Destination_Filename as a double


% --- Executes during object creation, after setting all properties.
function Destination_Filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Destination_Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Destination_counter_Callback(hObject, eventdata, handles)
% hObject    handle to Destination_counter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Destination_counter as text
%        str2double(get(hObject,'String')) returns contents of Destination_counter as a double


% --- Executes during object creation, after setting all properties.
function Destination_counter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Destination_counter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function Destination_Ending_Callback(hObject, eventdata, handles)
% hObject    handle to Destination_Ending (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Destination_Ending as text
%        str2double(get(hObject,'String')) returns contents of Destination_Ending as a double


% --- Executes during object creation, after setting all properties.
function Destination_Ending_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Destination_Ending (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in Stop.
function Stop_Callback(hObject, eventdata, handles)
% hObject    handle to Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in copy.
function copy_Callback(hObject, eventdata, handles)
% hObject    handle to copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=handles.cp_list;
handles.two_images=get(handles.focal_pair,'Value');

tmp=strrep(strrep(strrep(datestr(now),' ','_'),'-','_'),':','_');
path_t=get(handles.Path_Destination_High,'String');
if (handles.two_images==0)
    path_t=get(handles.Path_Destination_low,'String');
end;

%save([path_t tmp]','h');

h = waitbar(0,'Please wait...');



for i=1:size(handles.cp_list,1)
    if ((isempty(handles.cp_list{i,1})==0) )
        
        if (handles.two_images==1)
            copyfile(handles.cp_list{i,1},handles.cp_list{i,2});
        end;
        copyfile(handles.cp_list{i,3},handles.cp_list{i,4});
        
        waitbar(i/size(handles.cp_list,1));
        
        disp('---------------------------------------------------------------------------------------------------------------------------');
         if (handles.two_images==1)
            disp(['copy ' handles.cp_list{i,1} '    to    ' handles.cp_list{i,2}]);
        end;
        disp(['copy ' handles.cp_list{i,3} '    to    ' handles.cp_list{i,4}]);
        disp('---------------------------------------------------------------------------------------------------------------------------');
        fprintf('  \n ');
        
        
 
    end;
end;
close(h);


% --- Executes on selection change in File_Thrash_low.
function File_Thrash_low_Callback(hObject, eventdata, handles)
% hObject    handle to File_Thrash_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns File_Thrash_low contents as cell array
%        contents{get(hObject,'Value')} returns selected item from File_Thrash_low


% --- Executes during object creation, after setting all properties.
function File_Thrash_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_Thrash_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Path_Thrash_low_Callback(hObject, eventdata, handles)
% hObject    handle to Path_Thrash_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Path_Thrash_low as text
%        str2double(get(hObject,'String')) returns contents of Path_Thrash_low as a double


% --- Executes during object creation, after setting all properties.
function Path_Thrash_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Path_Thrash_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox5.
function listbox5_Callback(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox5


% --- Executes during object creation, after setting all properties.
function listbox5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Path_Thrash_High_Callback(hObject, eventdata, handles)
% hObject    handle to Path_Thrash_High (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Path_Thrash_High as text
%        str2double(get(hObject,'String')) returns contents of Path_Thrash_High as a double


% --- Executes during object creation, after setting all properties.
function Path_Thrash_High_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Path_Thrash_High (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in File_Thrash_high.
function File_Thrash_high_Callback(hObject, eventdata, handles)
% hObject    handle to File_Thrash_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns File_Thrash_high contents as cell array
%        contents{get(hObject,'Value')} returns selected item from File_Thrash_high


% --- Executes during object creation, after setting all properties.
function File_Thrash_high_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_Thrash_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Path_Thrash_high_Callback(hObject, eventdata, handles)
% hObject    handle to Path_Thrash_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Path_Thrash_high as text
%        str2double(get(hObject,'String')) returns contents of Path_Thrash_high as a double


% --- Executes during object creation, after setting all properties.
function Path_Thrash_high_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Path_Thrash_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in File_Used_low.
function File_Used_low_Callback(hObject, eventdata, handles)
% hObject    handle to File_Used_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns File_Used_low contents as cell array
%        contents{get(hObject,'Value')} returns selected item from File_Used_low


% --- Executes during object creation, after setting all properties.
function File_Used_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_Used_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in File_Used_high.
function File_Used_high_Callback(hObject, eventdata, handles)
% hObject    handle to File_Used_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns File_Used_high contents as cell array
%        contents{get(hObject,'Value')} returns selected item from File_Used_high


% --- Executes during object creation, after setting all properties.
function File_Used_high_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_Used_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)
% hObject    handle to back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.Dest.counter>1)
    handles.Dest.counter=handles.Dest.counter-1;
else
    handles.Dest.counter=1;
end;
if (handles.Used.counter>1)
    handles.Used.counter=handles.Used.counter-1;
else
    handles.Used.counter=1;
end;
if (handles.Source.counter>1)
    handles.Source.counter=handles.Source.counter-1;   
else
    handles.Source.counter=1;
end;
if (handles.Thrash.counter>1)
    handles.Thrash.counter=handles.Thrash.counter-1;   
else
     handles.Thrash.counter=1;
end;

set(handles.File_Source_High,'String',[handles.Source.basename  num2str(handles.Source.counter)]);
set(handles.File_Source_low,'String',[handles.Source.basename  num2str(handles.Source.counter)]);
set(handles.Destination_counter,'String',num2str(handles.Dest.counter));


source_high_tmp=[handles.Source.path_high     handles.Source.basename  num2str(handles.Source.counter) handles.Source.ext];

im=tom_emread(source_high_tmp);
im=tom_bin(im.Value,3);
axes(handles.axes1);
tom_imagesc(im);

 guidata(hObject, handles);



% --------------------------------------------------------------------
function save_cp_Callback(hObject, eventdata, handles)
% hObject    handle to save_cp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name path]=uiputfile;
cp_list.list=handles.cp_list;
cp_list.source_counter=handles.Source.counter;
save([path name],'cp_list');



% --------------------------------------------------------------------
function load_cp_Callback(hObject, eventdata, handles)
% hObject    handle to load_cp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[name path]=uigetfile;
load([path name]);


if (isfield(handles,'cp_list')==1)
    if (isempty(handles.cp_list)==0)
        button = questdlg('Overwrite existing copy list');
        if (strcmp(button,'No') | strcmp(button,'Cancel'))
            return;
        end;
    end;
end;
    
handles.cp_list=cp_list.list;
guidata(hObject, handles);

% parse information of Source
sz=size(handles.cp_list,1);

st=handles.cp_list{sz,1};

source_path=[st(1:findstr(st,'low')-1) ];
if (isunix==1)
    source_path_low=[source_path  'low/'];
    source_path_high=[source_path  'high/'];
else
    source_path_low=[source_path 'low\'];
    source_path_high=[source_path 'high\'];
end;

name=[st(findstr(st,'low')+4:size(st,2))];
name=strrep(name,'\','');
name=strrep(name,'/','');

[basename rest]=strtok(name,'_');
basename=[basename '_'];
rest=strrep(rest,'_','');
[count ext]=strtok(rest,'.');
count=str2num(count);

source_basename=basename;
source_startname=[basename num2str(count)];
source_ext=ext;
source_count=count;

% parse information of Source
st=handles.cp_list{sz,2};
dest_path=[st(1:findstr(st,'low')-1) ];
if (isunix==1)
    dest_path_low=[dest_path 'low/'];
    dest_path_high=[dest_path 'high/'];
else
    dest_path_low=[dest_path 'low\'];
    dest_path_high=[dest_path 'high\'];
end;

name=[st(findstr(st,'low')+4:size(st,2))];
name=strrep(name,'\','');
name=strrep(name,'/','');

[basename rest]=strtok(name,'_');
basename=[basename '_'];
rest=strrep(rest,'_','');
[count ext]=strtok(rest,'.');
count=str2num(count);

name=[st(findstr(st,'low')+4:size(st,2))];
name=strrep(name,'\','');
name=strrep(name,'/','');

[basename rest]=strtok(name,'_');
basename=[basename '_'];
rest=strrep(rest,'_','');
[count ext]=strtok(rest,'.');
count=str2num(count);

dest_basename=basename;
dest_startname=[basename num2str(count)];
dest_ext=ext;
dest_count=count+1;

%uptdate Gui source
set(handles.Path_Source_low,'String',source_path_low);
set(handles.Path_Source_High,'String',source_path_high);
set(handles.File_Source_low,'String',[source_basename num2str(cp_list.source_counter) dest_ext]);
set(handles.File_Source_High,'String',[source_basename num2str(cp_list.source_counter) dest_ext]);


%update Gui dest
set(handles.Path_Destination_low,'String',dest_path_low);
set(handles.Path_Destination_High,'String',dest_path_high);
set(handles.Destination_Ending,'String',dest_ext);
set(handles.Destination_Filename,'String',dest_basename);
set(handles.Destination_counter,'String',num2str(dest_count));


 

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in focal_pair.
function focal_pair_Callback(hObject, eventdata, handles)
% hObject    handle to focal_pair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of focal_pair


