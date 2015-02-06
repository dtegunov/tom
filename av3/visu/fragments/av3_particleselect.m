function varargout = av3_particleselect(varargin)
% AV3_PARTICLESELECT M-file for av3_particleselect.fig
%   motl       : motivelist
%   volume     : volume (tomogram!) containig the particles
%   nparticles : number of particles to be examined
%   iclass     : number assigned to class of "particles"
%   
%      AV3_PARTICLESELECT, by itself, creates a new AV3_PARTICLESELECT or raises the existing
%      singleton*.
%
%      H = AV3_PARTICLESELECT returns the handle to a new AV3_PARTICLESELECT or the handle to
%      the existing singleton*.
%
%      AV3_PARTICLESELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AV3_PARTICLESELECT.M with the given input arguments.
%
%      AV3_PARTICLESELECT('Property','Value',...) creates a new AV3_PARTICLESELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before av3_particleselect_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to av3_particleselect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help av3_particleselect

% Last Modified by GUIDE v2.5 04-Oct-2002 17:25:37

 % Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @av3_particleselect_OpeningFcn, ...
                   'gui_OutputFcn',  @av3_particleselect_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


%   initialization of variables
motl=varargin{1};
volume=varargin{2};
nparticles=varargin{3};
iclass=varargin{4};
if nparticles > size(motl, 2)
    nparticles = size(motl, 2)
end;
i = 1;
while i <= nparticles
    colormap gray;
    x = motl(8,i); y = motl(9,i);z = motl(10,i);
    %display cross
    Radius = 5; % size of cross
    uu = [x x x x-Radius x+Radius];vv = [y-Radius y+Radius y y y] ;      
    %line(uu,vv,'LineWidth',1,'Color',[0 1 1])%light blue 
    imagesc(volume(:,:,z));
    line(uu,vv,'LineWidth',1,'Color',[1 0 0])%light blue 
    drawnow;



% --- Executes just before av3_particleselect is made visible.
function av3_particleselect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to av3_particleselect (see VARARGIN)

% Choose default command line output for av3_particleselect
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes av3_particleselect wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = av3_particleselect_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
%function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function pushbutton1_Callback(hObject, eventdata, handles)
motl(20,i) = iclass;
    x = motl(8,i); y = motl(9,i);z = motl(10,i);
    %display cross
    uu = [x x x x-Radius x+Radius];vv = [y-Radius y+Radius y y y] ;      
    %line(uu,vv,'LineWidth',1,'Color',[0 1 1])%light blue 
    imagesc(volume(:,:,z));
    line(uu,vv,'LineWidth',1,'Color',[1 0 0])%light blue 
    drawnow;
    i = i+1;
    varargout{1} = motl;
    return;

% --- Executes on button press in pushbutton2.
%function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function pushbutton2_Callback(hObject, eventdata, handles)
if i <= nparticles
    x = motl(8,i); y = motl(9,i);z = motl(10,i);
    %display cross
    Radius = 5 % size of cross
    uu = [x x x x-Radius x+Radius];vv = [y-Radius y+Radius y y y] ;      
    %line(uu,vv,'LineWidth',1,'Color',[0 1 1])%light blue 
    imagesc(volume(:,:,z));
    line(uu,vv,'LineWidth',1,'Color',[1 0 0])%light blue 
    drawnow;  
    i = i+1;
else
    varargout{1} = motl;
    return;
end;


