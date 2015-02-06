function varargout = tom_picker(varargin)
% TOM_PICKER(VOLUME) a tool for interactive 3D particles picking 
%
%   motivelist = tom_picker(volume)
%
%   tom_picker(volume) This function can be used to pick particles from
%   tomograms interactively. The volumes is kept in memory and masked at 
%   coordinates of tentative particles. 
%   It is a GUI to scan through the volume in z-direction
%   and lets you adjust the contrast interactively. Additionally, a running
%   average in z-direction can be calculated to increase contrast.
%   All coordinates (x,y,z) are stored in a motivelist of dim 20xnparticles
%   in columns 8, 9, and 10.
%  
%   Example
%  ---------
%
%   a=tom_sphere([128 128 128],35,5);
%
%   motivelist = tom_picker(a);
%
%   motivelist     : matrix containing parameters of particles:
%
%   The following parameters are stored in the matrix MOTIVELIST of dimension (20, NPARTICLES):
%   column 
%      1         : Cross-Correlation Coefficient
%      2         : x-coordinate in full tomogram
%      3         : y-coordinate in full tomogram
%      4         : peak number%      5         : not occupied - e.g. use for # of tomogram of extracted
%                   particles
%      8         : x-coordinate in full tomogram
%      9         : y-coordinate in full tomogram
%      10        : z-coordinate in full tomogram
%      11        : x-shift in subvolume (AFTER rotation)
%      12        : y-shift in subvolume (AFTER rotation)
%      13        : z-shift in subvolume (AFTER rotation)
%      14        : x-shift in subvolume (BEFORE rotation)
%      15        : y-shift in subvolume (BEFORE rotation)
%      16        : z-shift in subvolume (BEFORE rotation)
%      17        : Phi
%      18        : Psi
%      19        : Theta 
%      20        : class no
%
%   See also   TOM_VOLXY TOM_CHOOSER TOM_INTERVOL
%
%   04/02/03   FF 
%   last revision 
%   07/22/03 FF
%    Copyright (c) 2004
%    TOM toolbox for Electron Tomography
%    Max-Planck-Institute for Biochemistry
%    Dept. Molecular Structural Biology
%    82152 Martinsried, Germany
%    http://www.biochem.mpg.de/tom   

%   PARAMETERS
%   

% TOM_VOLXY2 M-file for tom_volxy2.fig
%      TOM_VOLXY2, by itself, creates a new TOM_VOLXY2 or raises the existing
%      singleton*.
%
%      H = TOM_VOLXY2 returns the handle to a new TOM_VOLXY2 or the handle to
%      the existing singleton*.
%
%      TOM_VOLXY2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_VOLXY2.M with the given input arguments.
%
%      TOM_VOLXY2('Property','Value',...) creates a new TOM_VOLXY2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_volxy2_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_volxy2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_volxy2

% Last Modified by GUIDE v2.5 22-Dec-2009 14:56:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_picker_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_picker_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%End initialization code - DO NOT EDIT


% --- Executes just before tom_picker is made visible.
function tom_picker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_picker (see VARARGIN)
% Choose default command line output for tom_picker
handles.output = hObject;
%Update handles structure
%  error message ?




handles.data_cube=varargin{1};
%handles.radius=varargin{2};
handles.motl=zeros(20,1);
handles.ind = 1;
handles.slicexy = 1;
handles.sliceavg = 1;
clear varargin
guidata(hObject,handles);
[mean max min std]=tom_dev(handles.data_cube);
handles.DataScale=[mean-3*std mean+3*std];
dim_x=size(handles.data_cube,1);% handles.Size(1);
dim_y=size(handles.data_cube,2);%=handles.Size(2);
handles.Size=size(handles.data_cube);
handles.actualaxis=[1 dim_x 1 dim_y];
guidata(hObject,handles);
% generate all Min and Max
tmp_obj=findobj('Tag','XY_slider');
set(tmp_obj,'Min',1);
set(tmp_obj,'Max',handles.Size(3))
set(tmp_obj,'SliderStep',[1./handles.Size(3) 5./handles.Size(3)]);
tmp_obj=findobj('Tag','AVG_slider');
set(tmp_obj,'Min',1);
set(tmp_obj,'Max',handles.Size(3))
set(tmp_obj,'SliderStep',[1./handles.Size(3) 1./handles.Size(3)]);
updatehistogram_Callback(hObject, eventdata, handles);
%UPDATE_ALL_Callback(hObject, eventdata, handles);
% UIWAIT makes tom_picker wait for user response (see UIRESUME)
%uiwait(handles.figure1);
uiwait(gcf);

% --- Outputs from this function are returned to the command line.
function varargout = tom_picker_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout{1} = handles.motl;
delete(gcf);%close(gcf);

% Get default command line output from handles structure
%varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function XY_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XY_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject,'SliderStep',[1./16 1./16]);
set(hObject,'Value',[1]);


% --- Executes on slider movement.
function XY_slider_Callback(hObject, eventdata, handles)

%   slider callback for z-slice

% hObject    handle to XY_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%Data=handles.Data;

dim_z=handles.Size(3);
dim_x=handles.Size(1);
dim_y=handles.Size(2);
slicexy=round(get(hObject,'Value'));
tmp_obj=findobj('Tag','TEXT_XY');
set(tmp_obj,'String',num2str(slicexy));
tmp_obj=findobj('Tag','AVG_slider');
sliceavg=round(get(tmp_obj,'Value'));
rangex_min=handles.actualaxis(1);
rangex=handles.actualaxis(2)-rangex_min;
rangey_min=handles.actualaxis(3);
rangey=handles.actualaxis(4)-rangey_min;

if slicexy+sliceavg>dim_z; sliceavg=dim_z-slicexy;end;


tmp=handles.data_cube(rangex_min:rangex_min+rangex, rangey_min:rangey_min+rangey, slicexy:slicexy+sliceavg-1);



%tmp=double(tom_emreadinc(data_cube,[rangey_min rangex_min slicexy],[rangey rangex sliceavg-1]));

tmp=mean(tmp,3);
tmp_obj=findobj('Tag','XY_slice');
handles.disp = tmp_obj;
axes(tmp_obj);
imagesc(tmp',[handles.DataScale]);colormap gray;set(gcf,'DoubleBuffer','on');axis image; axis ij; 
clear tmp;
set(tmp_obj,'Tag','XY_slice');
handles.sliceavg = sliceavg;
handles.slicexy = slicexy;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function TEXT_XY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TEXT_YZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function TEXT_XY_Callback(hObject, eventdata, handles)
% hObject    handle to TEXT_YZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TEXT_YZ as text
%        str2double(get(hObject,'String')) returns contents of TEXT_YZ as a double

%Data=handles.Data;
dim_z=size(handles.Data,3); % bug - to be fixed

tmp_obj=findobj('Tag','TEXT_XY');
get(tmp_obj,'Value');
slicexy=round(str2double(eval(get(tmp_obj,'String'))));
if slicexy<1 slicexy=1; set(tmp_obj,'String',slicexy); end;
if slicexy>dim_z slicexy=dim_z; set(tmp_obj,'String',slicexy); end;
tmp_obj=findobj('Tag','XY_slider');
set(tmp_obj,'Value',slicexy);

tmp=squeeze(handles.Data(:,:,slicexy));
tmp_obj=findobj('Tag','XY_slice');
axes(tmp_obj);
tmp=imagesc(tmp',[handles.DataScale]);clear tmp;
colormap gray; 
%axis image;
handles.slicexy = slicexy;
set(tmp_obj,'Tag','XY_slice');
guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function AVG_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AVG_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function AVG_slider_Callback(hObject, eventdata, handles)
% hObject    handle to AVG_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliceavg=round(get(hObject,'Value'));
tmp_obj=findobj('Tag','TEXT_AVG');
set(tmp_obj,'String',num2str(sliceavg));
handles.sliceavg=sliceavg;
guidata(hObject,handles);
UPDATE_ALL_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function TEXT_AVG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TEXT_AVG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function TEXT_AVG_Callback(hObject, eventdata, handles)
% hObject    handle to TEXT_AVG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TEXT_AVG as text
%        str2double(get(hObject,'String')) returns contents of TEXT_AVG as a double

dim_z=handles.Size(3);

tmp_obj=findobj('Tag','TEXT_AVG');
get(tmp_obj,'Value');
sliceavg=round(eval(get(tmp_obj,'String')));
if sliceavg<1 sliceavg=1; set(tmp_obj,'String',sliceavg); end;
if sliceavg>dim_z sliceavg=dim_z; set(tmp_obj,'String',sliceavg); end;
tmp_obj=findobj('Tag','AVG_slider');
set(tmp_obj,'Value',sliceavg);
handles.sliceavg=sliceavg;
guidata(hObject,handles);
UPDATE_ALL_Callback(hObject, eventdata, handles);


% --- Executes on button press in updatehistogram.
function updatehistogram_Callback(hObject, eventdata, handles)
% hObject    handle to updatehistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


tmp_obj=findobj('Tag','Histogram');
axes(tmp_obj);
%[h,n]=tom_hist3d(handles.Histogram);
[h,n]=tom_hist3d(handles.data_cube);

handles.DataScale=[n(1)  n(size(n,2))];
h=200.*h./(100.*handles.Size(1).*handles.Size(2).*handles.Size(3));
bar(n,h);axis auto;

guidata(hObject,handles);
UPDATE_ALL_Callback(hObject, eventdata, handles);
set(tmp_obj,'Tag','Histogram');





% --- Executes on button press in sethistogram.
function sethistogram_Callback(hObject, eventdata, handles)
% hObject    handle to sethistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp_obj=findobj('Tag','Histogram');
axes(tmp_obj);
k = waitforbuttonpress;
     point1 = get(gca,'CurrentPoint');    % button down detected
     finalRect = rbbox;                   % return figure units
     point2 = get(gca,'CurrentPoint');    % button up detected
     point1 = point1(1,1:2);              % extract x and y
     point2 = point2(1,1:2);
     p1 = min(point1,point2);             % calculate locations
     offset = abs(point1-point2);         % and dimensions
     x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
     y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
     hold on;
     axis manual;
     
handles.DataScale=[x(1) x(2)];
guidata(hObject,handles);
set(gca,'Xlim',[x(1) x(2)]);
UPDATE_ALL_Callback(hObject, eventdata, handles);
set(tmp_obj,'Tag','Histogram');

function UPDATE_ALL_Callback(hObject, eventdata, handles)
% hObject    handle to TEXT_YZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TEXT_YZ as text
%        str2double(get(hObject,'String')) returns contents of TEXT_YZ as a double


dim_x=handles.Size(1);
dim_y=handles.Size(2);
dim_z=handles.Size(3);

tmp_obj=findobj('Tag','TEXT_XY');
get(tmp_obj,'Value');
slicexy=round(eval(get(tmp_obj,'String')));
tmp_obj=findobj('Tag','AVG_slider');
sliceavg=round(get(tmp_obj,'Value'));
if slicexy+sliceavg>dim_z; sliceavg=dim_z-slicexy;end;
tmp=handles.data_cube(1:dim_x, 1:dim_y, slicexy:slicexy+sliceavg-1);
%tmp=double(tom_emreadinc(data_cube,[1 1 slicexy],[dim_x-1 dim_y-1 sliceavg-1]));

tmp=mean(tmp,3);
tmp_obj=findobj('Tag','XY_slice');
axes(tmp_obj);
imagesc(tmp',[handles.DataScale]);axis image; axis ij;clear tmp; 
axis(handles.actualaxis);
colormap gray;
%axis image;

set(tmp_obj,'Tag','XY_slice');

tmp2_obj=findobj('Tag','XY_slider');
slicexy=round((get(tmp2_obj,'Value')));
axes(tmp_obj);

function Zoom_out_Callback(hObject, eventdata, handles)
dim_x=handles.Size(1);
dim_y=handles.Size(2);
dim_z=handles.Size(3);
tmp_obj=findobj('Tag','XY_slider');
slicexy=round(get(tmp_obj,'Value'));
tmp_obj=findobj('Tag','AVG_slider');
sliceavg=round(get(tmp_obj,'Value'));
if slicexy+sliceavg>dim_z; sliceavg=dim_z-slicexy;end;
tmp=handles.data_cube(1:dim_x, 1:dim_y, slicexy:slicexy+sliceavg-1);
%tmp=double(tom_emreadinc(data_cube,[1 1 slicexy],[dim_x-1 dim_y-1 sliceavg-1]));
tmp=mean(tmp,3);
tmp_obj=findobj('Tag','XY_slice');
axes(tmp_obj);
imagesc(tmp',[handles.DataScale]);set(gcf,'DoubleBuffer','on');clear tmp;
set(tmp_obj,'Tag','XY_slice');
handles.oldaxis=[1 dim_x 1 dim_y];
guidata(hObject,handles);



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp_obj=findobj('Tag','XY_slice');
axes(tmp_obj);
k = waitforbuttonpress;
     point1 = get(gca,'CurrentPoint');    % button down detected
     finalRect = rbbox;                   % return figure units
     point2 = get(gca,'CurrentPoint');    % button up detected
     point1 = point1(1,1:2);              % extract x and y
     point2 = point2(1,1:2);
     p1 = min(point1,point2);             % calculate locations
                                          % take care: x and y are
                                          % exchanged... FF 07/22/03
     offset = abs(point1-point2);         % and dimensions
     x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
     y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];

x=round(x);
y=round(y);
%     x=round(x); y=round(y);
handles.actualaxis=[x(1) x(2) y(1) y(3)];
axis([x(1) x(2) y(1) y(3)]);
guidata(hObject,handles);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp_obj=findobj('Tag','XY_slice');
axes(tmp_obj);
dim_z=handles.Size(3);
dim_x=handles.Size(1);
dim_y=handles.Size(2);
tmp_obj=findobj('Tag','XY_slider');
slicexy=round(get(tmp_obj,'Value'));
tmp_obj=findobj('Tag','AVG_slider');
sliceavg=round(get(tmp_obj,'Value'));
if slicexy+sliceavg>dim_z; sliceavg=dim_z-slicexy;end;
tmp=handles.data_cube(1:dim_x, 1:dim_y, slicexy:slicexy+sliceavg-1);
%tmp=double(tom_emreadinc(data_cube,[1 1 slicexy],[dim_x-1 dim_y-1
%sliceavg-1]));
tmp=mean(tmp,3);
tmp_obj=findobj('Tag','XY_slice');
axes(tmp_obj);
imagesc(tmp',[handles.DataScale]); set(gcf,'DoubleBuffer','on');clear tmp;
set(tmp_obj,'Tag','XY_slice');
handles.actualaxis=[1 dim_x 1 dim_y];
axis([1 dim_x 1 dim_y]);
guidata(hObject,handles);



% --- Executes on button press in pickit.
function pickit_Callback(hObject, eventdata, handles)
% hObject    handle to pickit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.radius = round(str2double(get(handles.editradius,'String')));
guidata(hObject,handles);


if isfield(handles,{'binningFactor'})
    binningFactor = handles.binningFactor; 
else
    binningFactor = 0;
end;
binningFactor = 2^binningFactor;

handles.motl(10,handles.ind)=handles.slicexy * binningFactor;
m=handles.output;
set(m,'Pointer','crosshair');
k = waitforbuttonpress;   
if k==0 %mouse button press
    point1 = get(gca,'CurrentPoint');    % button down detected
    pt = round(point1(1,1:2));
    handles.motl(9,handles.ind)= pt(1,2)*binningFactor;%fixed 07/07/03 FF
    handles.motl(8,handles.ind)= pt(1,1)*binningFactor;
    handles.motl(4,handles.ind)= handles.ind;
    dx = pt(1,1)-handles.radius; dy = pt(1,2)-handles.radius; dz =handles.slicexy-handles.radius;
    ux = pt(1,1)+handles.radius; uy = pt(1,2)+handles.radius; uz =handles.slicexy+handles.radius;
    % fixed FF 07/22/03
    handles.data_cube(dx:ux, dy:uy,dz:uz)= handles.data_cube(dx:ux, dy:uy,dz:uz) - ...
        tom_spheremask(handles.data_cube(dx:ux, dy:uy,dz:uz),handles.radius); %fixed 07/07/03 FF
    handles.ind = handles.ind + 1;
end;
guidata(hObject,handles);
set(m,'Pointer','arrow');


%UPDATE_ALL_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function editradius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
handles.radius = round(str2double(get(hObject,'String')));
guidata(hObject,handles);



function editradius_Callback(hObject, eventdata, handles)
% hObject    handle to editradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editradius as text
%        str2double(get(hObject,'String')) returns contents of editradius as a double
handles.radius = round(str2double(get(hObject,'String')));
guidata(hObject,handles);

% --- Executes on button press in byebye.
function byebye_Callback(hObject, eventdata, handles)
% hObject    handle to byebye (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(gcf);


%handles.ouput = handles.motl;
%guidata(hObject,handles);
%closereq;
%delete(gcf);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_motl_Callback(hObject, eventdata, handles)
% hObject    handle to load_motl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[motlFile motlPath] = uigetfile({'*.em'},'Select motif list stored as EM file.');

motl = tom_emread([motlPath motlFile]);
motl = motl.Value;

radius = inputdlg('Set masking radius.','Radius');
radius = str2double(radius);

motlLength = size(motl,2);

sphere = tom_spheremask(zeros(radius*2,radius*2,radius*2),radius);

if isfield(handles,{'binningFactor'})
    binningFactor = handles.binningFactor; 
else
    binningFactor = 0;
end;

binningFactor = 1 / (2^binningFactor);

for i = 1:motlLength
    
    x = motl(8,i)*binningFactor-radius;
    y = motl(9,i)*binningFactor-radius;
    z = motl(10,i)*binningFactor-radius;
    
    handles.data_cube = tom_paste(handles.data_cube,sphere,[x y z]);

end;

handles.ind = motlLength+1;
handles.motl = motl;

guidata(hObject,handles);
msgbox('MOTL loaded');




% --------------------------------------------------------------------
function merge_motls_Callback(hObject, eventdata, handles)
% hObject    handle to merge_motls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[motlFile motlPath] = uigetfile({'*.em'},'Select motif list stored as EM file.');

motl = tom_emread([motlPath motlFile]);
motl = motl.Value;

motl(:,end+1:end+size(handles.motl,2)) = handles.motl;

handles.motl = motl;

guidata(hObject,handles);





% --------------------------------------------------------------------
function save_motl_Callback(~, eventdata, handles)
% hObject    handle to save_motl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
directoryname = uigetdir(pwd,'Select directory for MOTL.');
filename = inputdlg('Select filename for current');

tom_emwrite([directoryname '/' filename{1}],handles.motl);


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function filter_volume_Callback(hObject, eventdata, handles)
% hObject    handle to filter_volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

kernelSize = str2double(inputdlg('Kernel size'));

if ~isfield(handles,{'originalVolume'})
    handles.originalVolume = handles.data_cube; 
end;

data_cube = tom_filter(handles.data_cube,kernelSize);

handles.data_cube = data_cube;

guidata(hObject,handles);

msgbox('Filtering done');



% --------------------------------------------------------------------
function reset_volume_Callback(hObject, eventdata, handles)
% hObject    handle to reset_volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data_cube = handles.originalVolume;

guidata(hObject,handles);


% --------------------------------------------------------------------
function rotate_volume_Callback(hObject, eventdata, handles)
% hObject    handle to rotate_volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
angles = str2mat(inputdlg('Define rotation as [phi psi theta]'));

if ~isfield(handles,{'originalVolume'})
    handles.originalVolume = handles.data_cube; 
end;

handles.data_cube = tom_rotate(handles.data_cube,angles);

guidata(hObject,handles);

msgbox('Rotation done');


% --------------------------------------------------------------------
function bin_volume_Callback(hObject, eventdata, handles)
% hObject    handle to bin_volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

binningFactor = str2double(inputdlg('How much binning do you want to apply?'));

if ~isfield(handles,{'originalVolume'})
    handles.originalVolume = handles.data_cube; 
end;

handles.binningFactor = binningFactor;
handles.data_cube = tom_bin(handles.data_cube,binningFactor);

guidata(hObject,handles);
msgbox('Binning done');

