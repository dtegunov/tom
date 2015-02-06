function varargout = tom_lamboot2(varargin)
%TOM_lamboot2 sends the lamboot command to a list of computer
%
%   varargout = tom_lamboot2(varargin)
%
%Test the network connection of the Rechenzentrum and Baumeister computer.
%The user can select the ones connected to the network and send a lamboot
%command. 
%This GUI just work under linux computer. Do not run it under windows.
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
%   ... = tom_lamboot2(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   UNIX, SYSTEM
%
%   created by WDN 03/07/05
%   updated by WDN 04/07/05
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
                   'gui_OpeningFcn', @tom_lamboot2_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_lamboot2_OutputFcn, ...
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


% --- Executes just before tom_lamboot2 is made visible.
function tom_lamboot2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_lamboot2 (see VARARGIN)

% Choose default command line output for tom_lamboot2
handles.output = hObject;
handles.NrOpteron = 5;
handles.NrCluster = 8;
%handles.NrRz=handles.NrSally+ handles.NrCathy+handles.NrOpteron;
a=['cairo   ';'dublin  ';'haifa   ';'jakarta '; 'kyoto   ';,...
   'lima    ';'montreal';'nairobi ';'tucson  ';'xian    ';'berlin  '];
b=['cluster01';'cluster02';'cluster03';'cluster04';'cluster05';,...
   'cluster06';'cluster07';'cluster08'];
c=['opteron  ';'opteron01';'opteron02';'opteron03';'fire5    '];
d=['cluster01';'cluster02';'cluster03';'cluster04';'cluster05';,...
   'cluster06';'cluster07';'cluster08';,...
   'opteron  ';'opteron01';'opteron02';'opteron03';'fire5    '];
handles.BaumeisterComputer=cellstr(a);
handles.RzCluster=cellstr(b);
handles.RzOpteron=cellstr(c);
handles.RzComputer=cellstr(d);
handles.Gray=get(hObject,'Color');%default gray
handles.Green=[0 1 0];%green
handles.Orange=[1 0.75 0];%orange
handles.Red=[1 0 0];%red
handles.Blue=[0 0.5 1];%blue
set(handles.status_1,'BackgroundColor',handles.Green);%green
set(handles.status_2,'BackgroundColor',handles.Orange);%orange
set(handles.status_3,'BackgroundColor',handles.Blue);%blue
set(handles.status_4,'BackgroundColor',handles.Red);%red
HowManyCpu_Callback(hObject, eventdata, handles);
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = tom_lamboot2_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%-----------------------------------------
% --- CHECKBOX ALL BAUMEISTER COMPUTER ---
%-----------------------------------------
function baumeister_all_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    for i=1:size(handles.BaumeisterComputer,1)
        set(handles.([handles.BaumeisterComputer{i}]),...
            'Value',1,...
            'BackgroundColor',handles.Gray);%default gray
    end
else
    for i=1:size(handles.BaumeisterComputer,1)
        set(handles.([handles.BaumeisterComputer{i}]),...
            'Value',0,...
            'BackgroundColor',handles.Gray);%default gray;
    end
end
HowManyCpu_Callback(hObject, eventdata, handles);
%----------------------------------------
% --- CHECKBOX BAUMEISTER ALL XOSVIEW ---
%----------------------------------------
function bau_xos_all_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    for i=1:size(handles.BaumeisterComputer,1)
        set(handles.(['xos_' handles.BaumeisterComputer{i}]),'Value',1);
    end
    unix('ssh -X cairo /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+0 &');pause(0.5);
    unix('ssh -X dublin /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+95 &');pause(0.5);
    unix('ssh -X haifa /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+164 &');pause(0.5);
    unix('ssh -X jakarta /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+233 &');pause(0.5);
    unix('ssh -X kyoto /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+302 &');pause(0.5);
    unix('ssh -X lima /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+371 &');pause(0.5);
    unix('ssh -X montreal /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+440 &');pause(0.5);
    unix('ssh -X nairobi /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+509 &');pause(0.5);
    unix('ssh -X tucson /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+578 &');pause(0.5);
    unix('ssh -X xian /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+647 &');pause(0.5);
    unix('ssh -X berlin /usr/bin/xosview -load -mem -swap -disk -page -int -net -geometry 220x482+230+0 &');
else
    for i=1:size(handles.BaumeisterComputer,1)
        set(handles.(['xos_' handles.BaumeisterComputer{i}]),'Value',0);
        unix(['ssh ' handles.BaumeisterComputer{i} ' killall xosview &']);
    end
end
%--------------------------------
% --- CHECKBOX XOSVIEW BERLIN ---
%--------------------------------
function xos_berlin_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X berlin /usr/bin/xosview -load -mem -swap -disk -page -int -net -geometry 220x482+230+20 &');
else
    unix('ssh berlin killall xosview &');
end

%-------------------------------
% --- CHECKBOX XOSVIEW CAIRO ---
%-------------------------------
function xos_cairo_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X cairo /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+0 &');
else
    unix('ssh cairo killall xosview &');
end

%--------------------------------
% --- CHECKBOX XOSVIEW DUBLIN ---
%--------------------------------
function xos_dublin_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X dublin /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+95 &');
else
    unix('ssh dublin killall xosview &');
end

%-------------------------------
% --- CHECKBOX XOSVIEW HAIFA ---
%-------------------------------
function xos_haifa_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X haifa /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+164 &');
else
    unix('ssh haifa killall xosview &');
end

%---------------------------------
% --- CHECKBOX XOSVIEW JAKARTA ---
%---------------------------------
function xos_jakarta_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X jakarta /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+233 &');
else
    unix('ssh jakarta killall xosview &');
end
%-------------------------------
% --- CHECKBOX XOSVIEW KYOTO ---
%-------------------------------
function xos_kyoto_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X kyoto /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+302 &');
else
    unix('ssh kyoto killall xosview &');
end
%------------------------------
% --- CHECKBOX XOSVIEW LIMA ---
%------------------------------
function xos_lima_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X lima /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+371 &');
else
    unix('ssh lima killall xosview &');
end

%----------------------------------
% --- CHECKBOX XOSVIEW MONTREAL ---
%----------------------------------
function xos_montreal_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X montreal /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+440 &');
else
    unix('ssh montreal killall xosview &');
end
%---------------------------------
% --- CHECKBOX XOSVIEW NAIROBI ---
%---------------------------------
function xos_nairobi_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X nairobi /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+509 &');
else
    unix('ssh nairobi killall xosview &');
end

%--------------------------------
% --- CHECKBOX XOSVIEW TUCSON ---
%--------------------------------
function xos_tucson_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X tucson /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+578 &');
else
    unix('ssh tucson killall xosview &');
end

%------------------------------
% --- CHECKBOX XOSVIEW XIAN ---
%------------------------------
function xos_xian_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X xian /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x42+0+647 &');
else
    unix('ssh xian killall xosview &');
end
%-----------------------------
% --- CHECKBOX ALL CLUSTER ---
%-----------------------------
function cluster_all_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.cluster01,'Value',1,'BackgroundColor',handles.Gray);%default gray
    set(handles.cluster02,'Value',1,'BackgroundColor',handles.Gray);%default gray
    set(handles.cluster03,'Value',1,'BackgroundColor',handles.Gray);%default gray
    set(handles.cluster04,'Value',1,'BackgroundColor',handles.Gray);%default gray
    set(handles.cluster05,'Value',1,'BackgroundColor',handles.Gray);%default gray
    set(handles.cluster06,'Value',1,'BackgroundColor',handles.Gray);%default gray
    set(handles.cluster07,'Value',1,'BackgroundColor',handles.Gray);%default gray
    set(handles.cluster08,'Value',1,'BackgroundColor',handles.Gray);%default gray
else
	set(handles.cluster01,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.cluster02,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.cluster03,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.cluster04,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.cluster05,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.cluster06,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.cluster07,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.cluster08,'Value',0,'BackgroundColor',handles.Gray);
end
HowManyCpu_Callback(hObject, eventdata, handles);
%-------------------------------------
% --- CHECKBOX CLUSTER ALL XOSVIEW ---
%-------------------------------------
function cluster_xos_all_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    for i=1:size(handles.RzCluster,1)
        set(handles.(['xos_' handles.RzCluster{i}]),'Value',1);
    end
    unix('ssh -X cluster01 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+460+0 &');pause(0.5);
    unix('ssh -X cluster02 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+460+375 &');pause(0.5);
    unix('ssh -X cluster03 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+690+0 &');pause(0.5);
    unix('ssh -X cluster04 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+690+375 &');pause(0.5);
    unix('ssh -X cluster05 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+920+0 &');pause(0.5);
    unix('ssh -X cluster06 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+920+375 &');pause(0.5);
    unix('ssh -X cluster07 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+1150+0 &');pause(0.5);
    unix('ssh -X cluster08 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+1150+375 &');pause(0.5);
else
    for i=1:size(handles.RzCluster,1)
        set(handles.(['xos_' handles.RzCluster{i}]),'Value',0);
        unix(['ssh ' handles.RzCluster{i} ' pkill xosview &']);
    end
end
%-----------------------------------
% --- CHECKBOX XOSVIEW cluster01 ---
%-----------------------------------
function xos_clu01_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X cluster01 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+460+0 &');
else
    unix('ssh cluster01 pkill xosview &');
end
%-----------------------------------
% --- CHECKBOX XOSVIEW cluster02 ---
%-----------------------------------
function xos_clu02_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X cluster02 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+460+375 &');
else
    unix('ssh cluster02 pkill xosview &');
end
%-----------------------------------
% --- CHECKBOX XOSVIEW cluster03 ---
%-----------------------------------
function xos_cluster03_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X cluster03 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+690+0 &');
else
    unix('ssh cluster03 pkill xosview &');
end
%-----------------------------------
% --- CHECKBOX XOSVIEW cluster04 ---
%-----------------------------------
function xos_clu04_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X cluster04 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+690+375 &');
else
    unix('ssh cluster04 pkill xosview &');
end
%-----------------------------------
% --- CHECKBOX XOSVIEW cluster05 ---
%-----------------------------------
function xos_cluster05_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X cluster05 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+920+0 &');
else
    unix('ssh cluster05 pkill xosview &');
end
%-----------------------------------
% --- CHECKBOX XOSVIEW cluster06 ---
%-----------------------------------
function xos_clu06_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X cluster06 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+920+375 &');
else
    unix('ssh cluster06 pkill xosview &');
end
%-----------------------------------
% --- CHECKBOX XOSVIEW cluster07 ---
%-----------------------------------
function xos_cluster07_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh -X cluster07 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+1150+0 &');
else
    unix('ssh cluster07 pkill xosview &');
end
%-----------------------------------
% --- CHECKBOX XOSVIEW cluster08 ---
%-----------------------------------
function xos_clu08_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh cluster08 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 220x322+1150+375 &');
else
    unix('ssh cluster08 pkill xosview &');
end
%-----------------------------
% --- CHECKBOX ALL OPTERON ---
%-----------------------------
function opteron_all_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.opteron,'Value',1,'BackgroundColor',handles.Gray);%default gray
    set(handles.opteron01,'Value',1,'BackgroundColor',handles.Gray);%default gray
    set(handles.opteron02,'Value',1,'BackgroundColor',handles.Gray);%default gray
    set(handles.opteron03,'Value',1,'BackgroundColor',handles.Gray);%default gray
    set(handles.fire5,'Value',1,'BackgroundColor',handles.Gray);%default gray
else
	set(handles.opteron,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.opteron01,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.opteron02,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.opteron03,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.fire5,'Value',0,'BackgroundColor',handles.Gray);
end
HowManyCpu_Callback(hObject, eventdata, handles);
%-------------------------------------
% --- CHECKBOX OPTERON ALL XOSVIEW ---
%-------------------------------------
function opteron_xos_all_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    for i=1:size(handles.RzOpteron,1)
        set(handles.(['xos_' handles.RzOpteron{i}]),'Value',1);
    end
    unix('ssh -X opteron /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 240x42+1380+0 &');pause(0.5);
    unix('ssh -X opteron01 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 240x42+1380+95 &');pause(0.5);
    unix('ssh -X opteron02 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 240x42+1380+164 &');pause(0.5);
    unix('ssh -X opteron03 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 240x42+1380+233 &');pause(0.5);
    unix('ssh -X fire5 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 240x42+1380+302 &');pause(0.5);
else
    for i=1:size(handles.RzOpteron,1)
        set(handles.(['xos_' handles.RzOpteron{i}]),'Value',0);
        unix(['ssh ' handles.RzOpteron{i} ' pkill xosview &']);
    end
end
%---------------------------------
% --- CHECKBOX XOSVIEW OPTERON ---
%---------------------------------
function xos_opt_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh opteron /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 240x42+1380+0 &');
else
    unix('ssh opteron pkill xosview &');
end
%-----------------------------------
% --- CHECKBOX XOSVIEW OPTERON01 ---
%-----------------------------------
function xos_opt01_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh opteron01 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 240x42+1380+95 &');
else
    unix('ssh opteron01 pkill xosview &');
end
%-----------------------------------
% --- CHECKBOX XOSVIEW OPTERON02 ---
%-----------------------------------
function xos_opt02_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh opteron02 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 240x42+1380+164 &');
else
    unix('ssh opteron02 pkill xosview &');
end
%-----------------------------------
% --- CHECKBOX XOSVIEW OPTERON03 ---
%-----------------------------------
function xos_opt03_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh opteron03 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 240x42+1380+233 &');
else
    unix('ssh opteron03 pkill xosview &');
end
%-------------------------------
% --- CHECKBOX XOSVIEW FIRE5 ---
%-------------------------------
function xos_fire5_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    unix('ssh fire5 /usr/bin/X11/xosview -load -mem -swap -disk -page -int -net -geometry 240x42+1380+302 &');
else
    unix('ssh fire5 pkill xosview &');
end
%----------------------------------------------------------------
% --- BUTTON TEST NETWORK CONNECTION RZ & BAUMEISTER COMPUTER ---
%----------------------------------------------------------------
function test_all_Callback(hObject, eventdata, handles)
set(handles.message,'String','Wait a moment. Processing ....');
for i=1:size(handles.RzComputer,1)
    set(handles.(handles.RzComputer{i}),'BackgroundColor',handles.Gray);%default gray
end
for i=1:size(handles.BaumeisterComputer,1)
    set(handles.([handles.BaumeisterComputer{i}]),'BackgroundColor',handles.Gray);%default gray
end
%%%%%%%%%%  test baumeister computer  %%%%%%%%%%
for i=1:size(handles.BaumeisterComputer,1)
    if get(handles.([handles.BaumeisterComputer{i}]),'Value')
        [status,message]=unix(['ping -w 1 -c 1 ' handles.BaumeisterComputer{i}]);
        if status~=1
            [status,message]=unix(['ssh ' handles.BaumeisterComputer{i} ' uname -a']);
            if status==0
                set(handles.([handles.BaumeisterComputer{i}]),'BackgroundColor',handles.Green);%green
            else
                set(handles.([handles.BaumeisterComputer{i}]),'Value',0,...
                    'BackgroundColor',handles.Red);%red
            end
        else
            set(handles.([handles.BaumeisterComputer{i}]),...
                'Value',0,...
                'BackgroundColor',handles.Red);%red
            set(handles.baumeister_all,'Value',0);
        end
        drawnow;
    end
end
%%%%%%%%%%  test RZ computer %%%%%%%%%%
for i=1:size(handles.RzComputer,1)
    if get(handles.(handles.RzComputer{i}),'Value')
        [status,message]=unix(['ping -w 1 -c 1 ' handles.RzComputer{i}]);
        if status~=1
            [status,message]=unix(['ssh ' handles.RzComputer{i} ' uname -a']);
            if status==0
                set(handles.(handles.RzComputer{i}),'BackgroundColor',handles.Green);%green
            else
                set(handles.(handles.RzComputer{i}),'Value',0,...
                    'BackgroundColor',handles.Red);%red
            end
        else
            set(handles.(handles.RzComputer{i}),...
                'Value',0,...
                'BackgroundColor',handles.Red);%red
        end
        drawnow;
    end
end
HowManyCpu_Callback(hObject, eventdata, handles);
set(handles.message,'String','');
%--------------------------------------------------
% --- BUTTON RESET ALL RZ & BAUMEISTER COMPUTER ---
%--------------------------------------------------
function reset_all_Callback(hObject, eventdata, handles)
%%%%%%%%%%  reset_all RZ  %%%%%%%%%%
set(handles.opteron_all,'Value',0);
set(handles.cluster_all,'Value',0);
for i=1:size(handles.RzComputer,1)
    set(handles.([handles.RzComputer{i}]),...
        'Value',0,...
        'Enable','on',...
        'BackgroundColor',handles.Gray);%default gray
end
set(handles.cluster_xos_all,'Value',0);
cluster_xos_all_Callback(handles.cluster_xos_all, eventdata, handles);
set(handles.opteron_xos_all,'Value',0);
opteron_xos_all_Callback(handles.opteron_xos_all, eventdata, handles);
%%%%%%%%%%  reset_all baumeister  %%%%%%%%%%
set(handles.baumeister_all,'Value',0);
for i=1:size(handles.BaumeisterComputer,1)
    set(handles.([handles.BaumeisterComputer{i}]),...
        'Value',0,...
        'Enable','on',...
        'BackgroundColor',handles.Gray);%default gray
end
set(handles.cpu_tot,'String','');
set(handles.bau_xos_all,'Value',0);
bau_xos_all_Callback(handles.bau_xos_all, eventdata, handles);
%-----------------------------
% --- BUTTON START LAMBOOT ---
%-----------------------------
function Start_lamb_Callback(hObject, eventdata, handles)
[s,w]=unix('whoami');
pathname=['/home/' fliplr(deblank(fliplr(deblank(w))))];
[s,w]=unix('hostname');computername=fliplr(deblank(fliplr(deblank(w))));
fid=fopen([pathname  '/lamhosts'],'w');
%%%%%%%%%%  RZ computer  %%%%%%%%%%
cn=0;
for i=1:size(handles.RzComputer,1)
    if get(handles.([handles.RzComputer{i}]),'Value')
        if strcmp([handles.RzComputer{i}],computername)
            cn=1;
        end
        fprintf(fid,'%s\n',[handles.RzComputer{i} ' cpu=' get(handles.(['cpu_' handles.RzComputer{i}]),'String')]);
    end
end
%%%%%%%%%%  baumeister computer  %%%%%%%%%%
for i=1:size(handles.BaumeisterComputer,1)
    if get(handles.([handles.BaumeisterComputer{i}]),'Value')
        if strcmp([handles.BaumeisterComputer{i}],computername)
            cn=1;
        end
        fprintf(fid,'%s\n',[handles.BaumeisterComputer{i} ' cpu=' get(handles.(['cpu_' handles.BaumeisterComputer{i}]),'String')]);
    end
end
fclose(fid);
if cn==0
    message=['Error. The command lamboot requires to have the computer ' computername ' selected. So please select it and press again the button Start Lamboot.'];
    msgbox(message,'Lamboot error','error');
    return;
end
path_org=pwd;
cd(pathname);
[s,w]=unix('chmod 755 lamhosts');
[s,w]=unix('ls -l lamhosts');
disp('');disp('********** lamhosts **********');
[s,w]=unix('cat lamhosts');disp(w);
[sta,out]=unix('lamclean');
[sta,out]=unix('lamboot -v lamhosts');disp(w);
disp('********** lamboot ready **********');
cd(path_org);


%-----------------------
% --- BUTTON REFRESH ---
%-----------------------
function Refresh_Callback(hObject, eventdata, handles)
%%%%%%%%%%  refresh RZ computer  %%%%%%%%%%
set(handles.message,'String','Wait a moment. Processing ....');
drawnow;
list=[ ];j=1;
for i=1:size(handles.RzComputer,1);
    if get(handles.([handles.RzComputer{i}]),'Value')==1
        list{j}=handles.RzComputer{i};
        j=j+1;
    end
end
list=list'; j=j-1;       
out=tom_systeminfo(list);
for i=1:j
    if out(i).Accessible==0;
        set(handles.(out(i).Name),'BackgroundColor',handles.Red);%red
    elseif out(i).LoadAverage(2)<1;
        %cpu free
        set(handles.(out(i).Name),'BackgroundColor',handles.Green);%green
    elseif  out(i).LoadAverage(2)>1 & out(i).LoadAverage(2)<=2.5;
        %cpu occupied but not full
        set(handles.(out(i).Name),'BackgroundColor',handles.Orange);%orange
    else out(i).LoadAverage(2)>2.5;
        %cpu is full
        set(handles.(out(i).Name),'BackgroundColor',handles.Blue);%blue
    end
end

%%%%%%%%%%  refresh baumeister computer  %%%%%%%%%%
list=[ ];j=1;
for i=1:size(handles.BaumeisterComputer,1);
    if get(handles.([handles.BaumeisterComputer{i}]),'Value')==1
        list{j}=handles.BaumeisterComputer{i};
        j=j+1;
    end
end
list=list'; j=j-1;       
out=tom_systeminfo(list);
for i=1:j
    if out(i).Accessible==0;
        set(handles.(out(i).Name),'BackgroundColor',handles.Red);%red
    elseif out(i).LoadAverage(2)<1;
        %cpu free
        set(handles.(out(i).Name),'BackgroundColor',handles.Green);%green
    elseif  out(i).LoadAverage(2)>1 & out(i).LoadAverage(2)<=2.5;
        %cpu occupied but not full
        set(handles.(out(i).Name),'BackgroundColor',handles.Orange);%orange
    else out(i).LoadAverage(2)>2.5;
        %cpu is full
        set(handles.(out(i).Name),'BackgroundColor',handles.Blue);%blue
    end
end
set(handles.message,'String','');
%if out.Accessible==0;
%    set(handles.(handles.BaumeisterComputer{i}),'BackgroundColor',handles.Red);%red
%elseif out.LoadAverage(1)<1; 
    %cpu free
%    set(handles.(handles.BaumeisterComputer{i}),'BackgroundColor',handles.Green);%green
%elseif  out.LoadAverage(1)>1 & out.LoadAverage(1)<=2.5;
    %cpu occupied but not full
%    set(handles.(handles.BaumeisterComputer{i}),'BackgroundColor',handles.Orange);%orange
%else out.LoadAverage(1)>2.5;
    %cpu is full
%    set(handles.(handles.BaumeisterComputer{i}),'BackgroundColor',handles.Blue);%Blue
%end


%----------------------------
% --- FUNCTION HOWMANYCPU ---
%----------------------------
function HowManyCpu_Callback(hObject, eventdata, handles)
Nbcpu=0;
%%%%%%%%%%  RZ computer  %%%%%%%%%%
for i=1:size(handles.RzComputer,1)
    if get(handles.([handles.RzComputer{i}]),'Value')
        Nbcpu=Nbcpu + str2num(get(handles.(['cpu_' handles.RzComputer{1}]),'String'));
    end
end
%%%%%%%%%%  baumeister computer  %%%%%%%%%%
for i=1:size(handles.BaumeisterComputer,1)
    if get(handles.([handles.BaumeisterComputer{i}]),'Value')
        Nbcpu=Nbcpu + str2num(get(handles.(['cpu_' handles.BaumeisterComputer{i}]),'String'));
    end
end

set(handles.cpu_tot,'String',num2str(Nbcpu));

