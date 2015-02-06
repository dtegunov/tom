function varargout = tom_lamboot(varargin)
%TOM_LAMBOOT sends the lamboot command to a list of computer
%
%   varargout = tom_lamboot(varargin)
%
%Test the network connection of the Rechenzentrum and Baumeister computer.
%The user can select the ones connected to the network and send a lamboot
%command. 
%This GUI just work under linux computer. Do not run it under windows.
%
%PARAMETERS
%
%  INPUT
%   filelabel           ...
%   threshold           ...
%   label               ...
%   color               ...
%   transformmatrix     ...
%   iconposition        ...
%   host                ...
%  
%  OUTPUT
%   data		...
%
%EXAMPLE
%   ... = tom_lamboot(...);
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
                   'gui_OpeningFcn', @tom_lamboot_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_lamboot_OutputFcn, ...
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


% --- Executes just before tom_lamboot is made visible.
function tom_lamboot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_lamboot (see VARARGIN)

% Choose default command line output for tom_lamboot
handles.output = hObject;
%handles.NrSally = 22;
%handles.NrCathy = 9;
%handles.NrOpteron = 5;
handles.NrCluster = 4;
%handles.NrRz=handles.NrSally+ handles.NrCathy+handles.NrOpteron;
a=['berlin  ';'cairo   ';'dublin  ';'haifa   ';'jakarta '; 'kyoto   ';,...
   'lima    ';'montreal';'nairobi ';'tucson  ';'xian    '];
% b=['sally1   ';'sally2   ';'sally3   ';'sally4   ';'sally5   ';'sally6   ';'sally7   ';,...
%    'sally8   ';'sally9   ';'sally10  ';'sally11  ';'sally12  ';'sally13  ';'sally14  ';,...
%    'sally15  ';'sally16  ';'sally17  ';'sally18  ';'sally19  ';'sally20  ';'sally21  ';,...
%    'sally22  ';'cathy1   ';'cathy2   ';'cathy3   ';'cathy4   ';'cathy5   ';'cathy6   ';,...
%    'cathy7   ';'cathy8   ';'cathy9   ';'opteron  ';'opteron01';'opteron02';,...
%    'opteron03';'fire5    ';'cluster01';'cluster02';'cluster03';'cluster04';'cluster05'];
b=['cluster09';'cluster10';'cluster11';'cluster12'];
handles.BaumeisterComputer=cellstr(a);
handles.RzComputer=cellstr(b);
handles.Gray=get(hObject,'Color');%default gray
handles.Green=[0 1 0];%green
handles.Orange=[1 0.75 0];%orange
handles.Red=[1 0 0];%red
set(handles.status_1,'BackgroundColor',handles.Green);%green
set(handles.status_2,'BackgroundColor',handles.Orange);%orange
set(handles.status_3,'BackgroundColor',handles.Red);%red
HowManyCpu_Callback(hObject, eventdata, handles);
% Update handles structure
guidata(hObject, handles);


% UIWAIT makes tom_lamboot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_lamboot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%---------------------------
% --- CHECKBOX ALL SALLY ---
%---------------------------
function sally_all_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    for i=1:handles.NrSally
        set(handles.(['sally' num2str(i)]),...
            'Value',1,...
            'BackgroundColor',handles.Gray);%default gray
    end
else
    for i=1:handles.NrSally
        set(handles.(['sally' num2str(i)]),...
            'Value',0,...
            'BackgroundColor',handles.Gray);%default gray;
    end    
end
HowManyCpu_Callback(hObject, eventdata, handles);
%---------------------------
% --- CHECKBOX ALL CATHY ---
%---------------------------
function cathy_all_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    for i=1:handles.NrCathy
        set(handles.(['cathy' num2str(i)]),...
            'Value',1,...
            'BackgroundColor',handles.Gray);%default gray
    end
else
    for i=1:handles.NrCathy
        set(handles.(['cathy' num2str(i)]),...
            'Value',0,...
            'BackgroundColor',handles.Gray);%default gray
    end    
end
HowManyCpu_Callback(hObject, eventdata, handles);
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
else
	set(handles.cluster01,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.cluster02,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.cluster03,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.cluster04,'Value',0,'BackgroundColor',handles.Gray);
    set(handles.cluster05,'Value',0,'BackgroundColor',handles.Gray);
end
HowManyCpu_Callback(hObject, eventdata, handles);

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

%---------------------------------------------
% --- BUTTON TEST RZ & BAUMEISTER COMPUTER ---
%---------------------------------------------
function test_all_Callback(hObject, eventdata, handles)
%%%%%%%%%%  test RZ computer %%%%%%%%%%
for i=1:size(handles.RzComputer,1)
    set(handles.(handles.RzComputer{i}),'BackgroundColor',handles.Gray);%default gray
end
for i=1:size(handles.RzComputer,1)
    if get(handles.(handles.RzComputer{i}),'Value')
        [status,message]=unix(['ping -w 1 -c 1 ' handles.RzComputer{i}]);
        if status~=1
            set(handles.(handles.RzComputer{i}),'BackgroundColor',handles.Green);%green
        else
            set(handles.(handles.RzComputer{i}),...
                'Value',0,...
                'BackgroundColor',handles.Red);%red
            if i>=1 & i<=handles.NrSally
                set(handles.sally_all,'Value',0);
            elseif i>handles.NrSally & i<=(handles.NrSally+handles.NrCathy)
                set(handles.cathy_all,'Value',0);
            elseif i>(handles.NrSally+handles.NrCathy) & i<=(handles.NrSally+handles.NrCathy+handles.NrOpteron+handles.NrCluster)
                set(handles.opteron_all,'Value',0);
            end
        end
        drawnow;
    end
end
%%%%%%%%%%  test baumeister computer  %%%%%%%%%%
for i=1:size(handles.BaumeisterComputer,1)
    set(handles.([handles.BaumeisterComputer{i}]),'BackgroundColor',handles.Gray);%default gray
end
for i=1:size(handles.BaumeisterComputer,1)
    if get(handles.([handles.BaumeisterComputer{i}]),'Value')
        [status,message]=unix(['ping -w 1 -c 1 ' handles.BaumeisterComputer{i}]);
        if status~=1
            set(handles.([handles.BaumeisterComputer{i}]),'BackgroundColor',handles.Green);%green
        else
            set(handles.([handles.BaumeisterComputer{i}]),...
                'Value',0,...
                'BackgroundColor',handles.Red);%red
            set(handles.baumeister_all,'Value',0);
        end
        drawnow;
    end
end
HowManyCpu_Callback(hObject, eventdata, handles);

%--------------------------------------------------
% --- BUTTON RESET ALL RZ & BAUMEISTER COMPUTER ---
%--------------------------------------------------
function reset_all_Callback(hObject, eventdata, handles)
%%%%%%%%%%  reset_all RZ  %%%%%%%%%%
set(handles.sally_all,'Value',0);
set(handles.cathy_all,'Value',0);
set(handles.opteron_all,'Value',0);
for i=1:size(handles.RzComputer,1)
    set(handles.([handles.RzComputer{i}]),...
        'Value',0,...
        'Enable','on',...
        'BackgroundColor',handles.Gray);%default gray
end
%%%%%%%%%%  reset_all baumeister  %%%%%%%%%%
set(handles.baumeister_all,'Value',0);
for i=1:size(handles.BaumeisterComputer,1)
    set(handles.([handles.BaumeisterComputer{i}]),...
        'Value',0,...
        'Enable','on',...
        'BackgroundColor',handles.Gray);%default gray
end
set(handles.cpu_tot,'String','');

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







