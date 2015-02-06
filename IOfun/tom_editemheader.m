function varargout = tom_editemheader(varargin)
%TOM_EDITEMHEADER fills manually the Header of a EM Format image
%
%   varargout = tom_editemheader(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout           ...
%
%SYNTAX
%   tom_editemheader
%   tom_editemheader(Data)
%      
%DESCRIPTION
%   This is an interactive tool to fill manually the Header of an EM Format
%   image (Data)
%
%EXAMPLE
%   tom_editemheader
%   a=tom_emread('test.em');
%   tom_editemheader(a)
%
%REFERENCES
%
%SEE ALSO
%   TOM_EMREAD, TOM_EMWRITE
%
%   created by WDN 24/06/05
%   updated by WDN 10/02/06
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

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_editemheader_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_editemheader_OutputFcn, ...
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


% --- Executes just before tom_editemheader is made visible.
function tom_editemheader_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_editemheader (see VARARGIN)

% Choose default command line output for tom_editemheader
handles.output = hObject;
a=['byte    ';'short   ';'        ';'long int';'float   ';,...
    '        ';'        ';'complex '; 'double  '];
handles.DataType=cellstr(a);
if nargin==4
    handles.Data_org=varargin{1};
    handles.Data_mod=handles.Data_org;
    i_magic=handles.Data_org.Header.Magic(4);
    set(handles.size_type,'String',handles.DataType(i_magic));
    DisplayHeader(handles);
end
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = tom_editemheader_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%*************************************************************************
%------------------
% --- Menu LOAD ---
%------------------
function load_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile({'*.em';'*.*'}, 'Pick up an EM-file');
if isequal(filename,0) | isequal(pathname,0)
    disp('No data loaded.'); return;
end;
EM_name=[pathname filename];
handles.Data_org=tom_emread(EM_name);
handles.Data_mod=handles.Data_org;
i_magic=handles.Data_org.Header.Magic(4);
set(handles.size_type,'String',handles.DataType(i_magic));
DisplayHeader(handles);
guidata(hObject, handles);
%---------------------
% --- SET FILENAME ---
%---------------------
function filename_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Filename=get(hObject,'String');
guidata(hObject, handles);
%---------------------
% --- SET PATHNAME ---
%---------------------
function pathname_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Pathname=get(hObject,'String');
guidata(hObject, handles);
%-------------------
% --- SET SIZE X ---
%-------------------
function sizex_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Size(1)=str2num(get(hObject,'String'));
guidata(hObject, handles);
%-------------------
% --- SET SIZE Y ---
%-------------------
function sizey_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Size(2)=str2num(get(hObject,'String'));
guidata(hObject, handles);
%-------------------
% --- SET SIZE Z ---
%-------------------
function sizez_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Size(3)=str2num(get(hObject,'String'));
guidata(hObject, handles);
%--------------------
% --- SET VOLTAGE ---
%--------------------
function voltage_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Voltage=str2num(get(hObject,'String'));
guidata(hObject, handles);
%---------------
% --- SET CS ---
%---------------
function cs_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Cs=str2num(get(hObject,'String'));
guidata(hObject, handles);
%---------------------
% --- SET APERTURE ---
%---------------------
function aperture_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Aperture=str2num(get(hObject,'String'));
guidata(hObject, handles);
%--------------------------
% --- SET MAGNIFICATION ---
%--------------------------
function mag_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Magnification=str2num(get(hObject,'String'));
guidata(hObject, handles);
%------------------------------
% --- SET POSTMAGNIFICATION ---
%------------------------------
function postmag_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Postmagnification=str2num(get(hObject,'String'));
guidata(hObject, handles);
%--------------------------
% --- SET EXPOSURE TIME ---
%--------------------------
function Expo_time_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Exposuretime=str2num(get(hObject,'String'));
guidata(hObject, handles);
%-----------------------------
% --- SET OBJECT PIXELSIZE ---
%-----------------------------
function Obj_pix_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Objectpixelsize=str2num(get(hObject,'String'));
guidata(hObject, handles);
%-----------------------
% --- SET MICROSCOPE ---
%-----------------------
function microscope_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Microscope=get(hObject,'String');
guidata(hObject, handles);
%----------------------
% --- SET PIXELSIZE ---
%----------------------
function pixsize_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Pixelsize=str2num(get(hObject,'String'));
guidata(hObject, handles);
%--------------------
% --- SET CCDAREA ---
%--------------------
function ccdarea_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.CCDArea=str2num(get(hObject,'String'));
guidata(hObject, handles);
%--------------------
% --- SET DEFOCUS ---
%--------------------
function defocus_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Defocus=str2num(get(hObject,'String'));%save in Angstroem
guidata(hObject, handles);
%------------------------
% --- SET ASTIGMATISM ---
%------------------------
function astigmatism_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Astigmatism=str2num(get(hObject,'String'));
guidata(hObject, handles);
%-----------------------------
% --- SET ASTIGMATISM ANGLE---
%-----------------------------
function astigmatism_angle_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.AstigmatismAngle=str2num(get(hObject,'String'));
guidata(hObject, handles);
%---------------------------
% --- SET FOCUS INCREMENT---
%---------------------------
function focusincrement_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.FocusIncrement=str2num(get(hObject,'String'));
guidata(hObject, handles);
%-------------------------------
% --- SET COUNT PER ELECTRON ---
%-------------------------------
function CountsPerElectron_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.CountsPerElectron=str2num(get(hObject,'String'));
guidata(hObject, handles);
%----------------------
% --- SET INTENSITY ---
%----------------------
function intensity_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Intensity=str2num(get(hObject,'String'));
guidata(hObject, handles);
%------------------------------
% --- SET ENERGY SLIT WIDTH ---
%------------------------------
function EnergySlitwidth_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.EnergySlitwidth=str2num(get(hObject,'String'));
guidata(hObject, handles);
%--------------------------
% --- SET ENERGY OFFSET ---
%--------------------------
function EnergyOffset_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.EnergyOffset=str2num(get(hObject,'String'));
guidata(hObject, handles);
%-----------------------
% --- SET TILT ANGLE ---
%-----------------------
function tiltangle_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Tiltangle=str2num(get(hObject,'String'));
guidata(hObject, handles);
%----------------------
% --- SET TILT AXIS ---
%----------------------
function tiltaxis_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Tiltaxis=str2num(get(hObject,'String'));
guidata(hObject, handles);
%---------------------
% --- SET MARKER X ---
%---------------------
function Marker_X_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Marker_X=str2num(get(hObject,'String'));
guidata(hObject, handles);
%---------------------
% --- SET MARKER Y ---
%---------------------
function Marker_Y_Callback(hObject, eventdata, handles)
handles.Data_mod.Header.Marker_Y=str2num(get(hObject,'String'));
guidata(hObject, handles);
%--------------------
% --- SET COMMENT ---
%--------------------
function comment_Callback(hObject, eventdata, handles)
a=get(hObject,'String');c='';
if size(a,1)>1
    for i=1:size(a,1)
        c=[c deblank(a(i,:))];
    end
    co=c;
else
    co=a;
end
handles.Data_mod.Header.Comment=co;
guidata(hObject, handles);
%---------------------
% --- BUTTON APPLY ---
%---------------------
function apply_Callback(hObject, eventdata, handles)
path_org=pwd;
cd(handles.Data_mod.Header.Pathname)
d=dir;dd=struct2cell(d);
listfile=dd(1,3:size(dd,2));
for i=1:size(listfile,2)
    if strcmp(handles.Data_mod.Header.Filename,listfile(i))==1
        a=fullfile(handles.Data_mod.Header.Pathname,listfile{i});
        message=['The file ' a ' exists already. Do you want to overwrite it anyway?'];
        Question=questdlg(message,'Modify the header','Yes','No','No');
        if strcmp(Question,'Yes')
            break;%do nothing
        elseif strcmp(Question,'No')
            [myname, mypathname] = uiputfile('*.em', 'SAVE FILE AS');
            if isequal(myname,0)|isequal(mypathname,0)
                cd(path_org);
                return;%nothing because Cancel is ckicked
            else
                if isempty(findstr(myname,'.em'))
                    myname=strcat(myname,'.em');
                end
                handles.Data_mod.Header.Pathname=mypathname;
                handles.Data_mod.Header.Filename=myname;
                handles.Data_org.Header=handles.Data_mod.Header;
                set(handles.pathname,'String',mypathname);
                set(handles.filename,'String',myname);
                break;
            end
        end
    end
end
a=fullfile(handles.Data_mod.Header.Pathname,handles.Data_mod.Header.Filename);
tom_emwrite(a,handles.Data_mod);
message=([a ' saved']);
msgbox(message,'Data written');
cd(path_org);
guidata(hObject, handles);
%----------------------
% --- BUTTON RESET ----
%----------------------
function Reset_Callback(hObject, eventdata, handles)
DisplayHeader(handles)
guidata(hObject, handles);
%----------------------
% --- BUTTON QUIT ---
%----------------------
function Quit_Callback(hObject, eventdata, handles)
delete(handles.MEH)

% --------------------------------------------------------------------
% --------------------------------------------------------------------
% ----------- OTHER FUNCTIONS ----------------------------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------

function DisplayHeader(handles)
warning off;
set(handles.filename,'String',handles.Data_org.Header.Filename);
set(handles.pathname,'String',handles.Data_org.Header.Pathname);
if isempty(handles.Data_org.Header.Size)
    set(handles.sizex,'String','');
    set(handles.sizey,'String','');
    set(handles.sizez,'String','');
else
    set(handles.sizex,'String',handles.Data_org.Header.Size(1));
    set(handles.sizey,'String',handles.Data_org.Header.Size(2));
    set(handles.sizez,'String',handles.Data_org.Header.Size(3));
end
set(handles.voltage,'String',handles.Data_org.Header.Voltage);
set(handles.cs,'String',handles.Data_org.Header.Cs);
set(handles.aperture,'String',handles.Data_org.Header.Aperture);
set(handles.mag,'String',handles.Data_org.Header.Magnification);
set(handles.postmag,'String',handles.Data_org.Header.Postmagnification);
set(handles.Expo_time,'String',handles.Data_org.Header.Exposuretime);
set(handles.Obj_pix,'String',handles.Data_org.Header.Objectpixelsize);
set(handles.microscope,'String',handles.Data_org.Header.Microscope);
set(handles.pixsize,'String',handles.Data_org.Header.Pixelsize);
set(handles.ccdarea,'String',handles.Data_org.Header.CCDArea);
set(handles.defocus,'String',handles.Data_org.Header.Defocus);
set(handles.astigmatism,'String',handles.Data_org.Header.Astigmatism);
set(handles.astigmatism_angle,'String',handles.Data_org.Header.AstigmatismAngle);
set(handles.focusincrement,'String',handles.Data_org.Header.FocusIncrement);
set(handles.CountsPerElectron,'String',handles.Data_org.Header.CountsPerElectron);
set(handles.intensity,'String',handles.Data_org.Header.Intensity);
set(handles.EnergySlitwidth,'String',handles.Data_org.Header.EnergySlitwidth);
set(handles.EnergyOffset,'String',handles.Data_org.Header.EnergyOffset);
set(handles.tiltangle,'String',handles.Data_org.Header.Tiltangle);
set(handles.tiltaxis,'String',handles.Data_org.Header.Tiltaxis);
set(handles.Marker_X,'String',handles.Data_org.Header.Marker_X);
set(handles.Marker_Y,'String',handles.Data_org.Header.Marker_Y);
set(handles.comment,'String',char(handles.Data_org.Header.Comment'));
if isempty(handles.Data_org.Header.Magic)
    set(handles.magic1,'String','');
    set(handles.magic2,'String','');
    set(handles.magic3,'String','');
    set(handles.magic4,'String','');
    set(handles.size_type,'String','');
else
    set(handles.magic1,'String',handles.Data_org.Header.Magic(1));
    set(handles.magic2,'String',handles.Data_org.Header.Magic(2));
    set(handles.magic3,'String',handles.Data_org.Header.Magic(3));
    set(handles.magic4,'String',handles.Data_org.Header.Magic(4));
end
warning on;






