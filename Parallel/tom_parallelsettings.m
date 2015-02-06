function varargout = tom_parallelsettings(varargin)
%TOM_PARALLELSETTINGS is a GUI for creating a parallel settings structure
%
%   varargout = tom_parallelsettings(varargin)
%
%TOM_PARALLELSETTINGS is a GUI for creating a parallel settings structure used
%by various functions
%
%PARAMETERS
%
%  INPUT
%   in_struct           input parallel settings structure (optional)
%  
%  OUTPUT
%   out_struct          output parallel structure
%               
%                       outstruct.jobmanager:     Name of the jobmanager
%                       outstruct.packageloss:    maximum allowed package loss [0..1]
%                       outstruct.number_of_tasks: number of tasks in which the job will be split
%                       outstruct.workers.min:     minimum number of workers to use
%                       outstruct.workers.max:     maximum number of workers to use
%                       outstruct.timeout:         timeout value for each task in seconds
%                       outstruct.restart_workers: 1 = restart workers before job
%                       run, 0 = do not restart workers
%
%EXAMPLE
%   ... = tom_parallelsettings(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 01/25/06
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


%GUI for creating a parallel settings structure used by various functions
%
%SYNTAX
%out_struct = tom_parallelsettings(in_struct)
%
% INPUT
% in_struct:    input parallel settings structure (optional)
%
% OUTPUT
% out_struct:   output parallel structure
%               
%               outstruct.jobmanager:     Name of the jobmanager
%               outstruct.packageloss:    maximum allowed package loss [0..1]
%               outstruct.number_of_tasks: number of tasks in which the job will be split
%               outstruct.workers.min:     minimum number of workers to use
%               outstruct.workers.max:     maximum number of workers to use
%               outstruct.timeout:         timeout value for each task in seconds
%               outstruct.restart_workers: 1 = restart workers before job run, 0 = do not restart workers
%
%Copyright (c) 2006
%TOM toolbox for Electron Tomography
%Max-Planck-Institute for Biochemistry
%Dept. Molecular Structural Biology
%82152 Martinsried, Germany
%http://www.biochem.mpg.de/tom
%
%Created: 25/01/06 AK

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_parallelsettings_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_parallelsettings_OutputFcn, ...
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Opening Function                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tom_parallelsettings_OpeningFcn(hObject, eventdata, handles, varargin)

if nargin > 3
    handles.parstruct = varargin{1};
else
    handles.parstruct = struct();
    handles.parstruct.jobmanager = 'default_jobmanager';
    handles.parstruct.packageloss = 0;
    handles.parstruct.number_of_tasks = 1;
    handles.parstruct.workers.min = 1;
    handles.parstruct.workers.max = 16;
    handles.parstruct.timeout = 3600;
    handles.parstruct.restart_workers = 0;
end

set(handles.par_packageloss,'String',num2str(handles.parstruct.packageloss));
set(handles.par_numberoftasks,'String',num2str(handles.parstruct.number_of_tasks));
set(handles.par_numworkers_min,'String',num2str(handles.parstruct.workers.min));
set(handles.par_numworkers_max,'String',num2str(handles.parstruct.workers.max));
set(handles.par_timeout,'String',num2str(handles.parstruct.timeout));
set(handles.par_restartworkers,'Value',handles.parstruct.restart_workers);

%get information about all jobmanagers
tmp_managers{1} = findResource('scheduler','type','jobmanager','lookupurl','cluster09');
all_managers = {};
if isempty(tmp_managers{1})
    try
        all_managers{1} = findResource('scheduler','type','jobmanager','lookupurl','muc');
    catch
    end;
    
else
     all_managers{1} = tmp_managers{1};
    try
        all_managers{2} = findResource('scheduler','type','jobmanager','lookupurl','muc');
    catch
    end;
end;

%    try
%    all_managers = {};
    %all_managers{1} = findResource('scheduler','type','jobmanager','lookupurl','cluster02');
%    all_managers{1} = findResource('scheduler','type','jobmanager','lookupurl','cluster09');    
 %   all_managers{2} = findResource('scheduler','type','jobmanager','lookupurl','vancouver');
%    all_managers{2} = findResource('scheduler','type','jobmanager','lookupurl','muc');
%catch
%    all_managers{1} = findResource('scheduler','type','jobmanager','lookupurl','vancouver');
%    all_managers{1} = findResource('scheduler','type','jobmanager','lookupurl','muc');

    %      err = lasterror;
%      if strcmp(err.identifier,'distcomp:findResource:LicenseUnavailable')
%          error('All distributed computing licenses are in use, please try again later.');
%          return;
%      else
%          error('Unknown error');
%          return;
%      end
%end

if isempty(all_managers)
    error('No jobmanager found');
    return;
end

jmstring = '';
for i=1:size(all_managers,2)
% for i=1:1
    if i == 1
        jmname = all_managers{1}.HostName;
    end
    handles.jmstruct(i) = get(all_managers{i});
    jmstring = strvcat(jmstring,all_managers{i}.HostName);
end

set(handles.par_jobmanager,'String',jmstring);
handles = update_jm_infobox(jmname,handles);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_parallelsettings wait for user response (see UIRESUME)
uiwait(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_parallelsettings_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure

varargout{1} = handles.output;
delete(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Jobmanager dropdown box                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par_jobmanager_Callback(hObject, eventdata, handles)

contents = get(hObject,'String');
if ischar(contents) && size(contents,1) == 1
    tmp = contents;
else
    tmp = deblank(contents(get(hObject,'Value'),:));
end
handles = update_jm_infobox(tmp,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Update jobmanager infobox                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = update_jm_infobox(jmname,handles);

jm = findResource('scheduler','type','jobmanager','lookupurl',jmname);

ip_help=get(jm, 'Hostaddress');

if iscell(ip_help)
    ip_help = ip_help{end};
end
%x =get(jm, 'Hostname');
infostring = ['Host:' get(jm, 'Hostname')];
infostring = strvcat(infostring,['IP Address: ',ip_help]);
infostring = strvcat(infostring,['state: ',get(jm, 'State')]);
infostring = strvcat(infostring,['Idle workers: ',num2str(get(jm, 'NumberOfIdleWorkers'))]);
infostring = strvcat(infostring,['Busy workers: ',num2str(get(jm, 'NumberOfBusyWorkers'))]);

set(handles.par_jobmanagerinfo,'String',infostring);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  number of tasks                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par_packageloss_Callback(hObject, eventdata, handles)

string = str2num(get(hObject,'String'));
if ~isnumeric(string) | string > 100 | string < 0
    errordlg('package loss must be a number between 0 and 100!');
    set(hObject,'String',num2str(handles.parstruct.packageloss));
    return;
end
set(hObject,'String',num2str(round(abs(string))));
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  number of tasks                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par_numberoftasks_Callback(hObject, eventdata, handles)

string = str2num(get(hObject,'String'));
if ~isnumeric(string)
    errordlg('number of tasks must be a number!');
    set(hObject,'String',num2str(handles.parstruct.number_of_tasks));
    return;
end
set(hObject,'String',num2str(round(abs(string))));
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  min workers                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par_numworkers_min_Callback(hObject, eventdata, handles)

string = str2num(get(hObject,'String'));
if ~isnumeric(string)
    errordlg('min workers value must be a number!');
    set(hObject,'String',num2str(handles.parstruct.workers.min));
    return;
end
set(hObject,'String',num2str(round(abs(string))));
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  max workers                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par_numworkers_max_Callback(hObject, eventdata, handles)

string = str2num(get(hObject,'String'));
if ~isnumeric(string)
    errordlg('max workers value must be a number!');
    set(hObject,'String',num2str(handles.parstruct.workers.max));
    return;
end
set(hObject,'String',num2str(round(abs(string))));
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Timeout                                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par_timeout_Callback(hObject, eventdata, handles)

string = str2num(get(hObject,'String'));
if ~isnumeric(string)
    errordlg('timeout value must be a number!');
    set(hObject,'String',num2str(handles.parstruct.timeout));
    return;
end
set(hObject,'String',num2str(round(abs(string))));
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Button OK                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function part_buttonok_Callback(hObject, eventdata, handles)

handles.parstruct.packageloss = str2num(get(handles.par_packageloss,'String'));
handles.parstruct.number_of_tasks = str2num(get(handles.par_numberoftasks,'String'));
handles.parstruct.workers.min = str2num(get(handles.par_numworkers_min,'String'));
handles.parstruct.workers.max = str2num(get(handles.par_numworkers_max,'String'));
handles.parstruct.timeout = str2num(get(handles.par_timeout,'String'));
handles.parstruct.restart_workers = get(handles.par_restartworkers,'Value');

contents = get(handles.par_jobmanager,'String');

%tmp = contents{get(handles.par_jobmanager,'Value')};
tmp = deblank(contents(get(handles.par_jobmanager,'Value'),:));
handles.parstruct.jobmanager = tmp;

handles.output = handles.parstruct;
guidata(hObject, handles);
uiresume(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Button Cancel                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function part_buttoncancel_Callback(hObject, eventdata, handles)

handles.output = '';
guidata(hObject, handles);
uiresume(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par_timeout_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function par_numworkers_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function par_numworkers_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function par_numberoftasks_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function par_jobmanager_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function par_packageloss_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function par_restartworkers_Callback(hObject, eventdata, handles)