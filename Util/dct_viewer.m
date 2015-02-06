function varargout = dct_viewer(varargin)
% DCT_VIEWER M-file for dct_viewer.fig
%      DCT_VIEWER, by itself, creates a new DCT_VIEWER or raises the existing
%      singleton.
%
%      H = DCT_VIEWER returns the handle to a new DCT_VIEWER or the handle to
%      the existing singleton.
%
%      DCT_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DCT_VIEWER.M with the given input
%      arguments.
%
%      DCT_VIEWER('Property','Value',...) creates a new DCT_VIEWER or raises the
%      existing singleton.  Starting from the left, property value pairs are
%      applied to the GUI before dct_viewer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dct_viewer_OpeningFcn via varargin.
%
%  Tim Farajian, Andy Jardine, 2005


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @dct_viewer_OpeningFcn, ...
    'gui_OutputFcn',  @dct_viewer_OutputFcn, ...
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




% --- GUI Output function ------------------------------------------
function varargout = dct_viewer_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;




% --- Close Figure -------------------------------------------------
function Figure_CloseRequestFcn(hObject, eventdata, handles)

try
    if strcmp(get(handles.updatetimer,'Running'),'on')
        stop(handles.updatetimer)
    end
    delete([handles.updatetimer handles.starttimer])
end
delete(hObject);




% --- Open Figure ---------------------------------------------------
function dct_viewer_OpeningFcn(hObject, eventdata, handles)

% Choose default command line output for dct_viewer
handles.output = hObject;

% define arrays in which to keep references to DCT elements on the network
handles.jm = [];
%handles.Workers = [];
handles.Jobs = [];
handles.Tasks = [];

% define arrays to keep reference to current elements
handles.cjm = [];
%handles.cWorker = [];
handles.cJob = [];
handles.cTask = [];

% Create update timer:
%   Specify timer period - Want smaller values on right of slider
%   Create timer, but dont start it - can be started when checkbox
%   ticked...

%mod = get(handles.RateSlider,'Max') + get(handles.RateSlider,'Min');
%SecPerUpdate =  mod - get(handles.RateSlider,'Value');
%set(handles.RateText,'String',sprintf('%0.1f',SecPerUpdate))
%oldwarn = warning('off','MATLAB:TIMER:RATEPRECISION');

SecPerUpdate = 1;
handles.updatetimer = timer('busymode','drop',...
    'ExecutionMode','fixedSpacing',...
    'Period',SecPerUpdate,'startdelay',0,...
    'TimerFcn',@(a,b)Update(a,b,hObject),...
    'TasksToExecute',1e6);

% % Create timer to start populating listboxes after figure is shown
% handles.starttimer = timer('startdelay',0.1,...
%     'TimerFcn',@PopulateJobManagers,'UserData',handles.Figure);

%warning(oldwarn); % Reset warning state


% LOAD PLUGINS:
%   Any M-file in the 'plugins' subdirectory will be listed in Plugins menu.
%   The plugin function will be passed a handle to the current job manager as
%   the only input.
%   Keep array of new Plugin menu items in handles.plugin
if exist('plugin','dir')
    handles.plugin = []; % Initialize demo handles
    olddir = cd('plugin'); % Change to demo directory
    files = dir('*.m'); % Read directory
    for n = 1:length(files)
        [ignore PluginName ext]=fileparts(files(n).name);
        handles.plugin(end+1) = uimenu(handles.MenuBarDemo,...
            'Label', PluginName, 'callback', @PluginCallback,...
            'UserData',{handles.Figure,str2func(PluginName)});
    end
    cd(olddir)
end

% put a flag in, showing the worker plot is not initialised...
handles.numWorkersPlotted = -1;
handles.selectedworker = -1;
handles.masterUpdate = 0;

% Record some resize information - this will be used to scale the GUI
% if we resize the window later
d = get(handles.Figure,'position');
handles.form_initwidth = d(3);
handles.form_initheight = d(4);
handles.pos_panel_jms = get(handles.panel_jms,'position');
handles.pos_panel_jobs = get(handles.panel_jobs,'position');
handles.pos_panel_tasks = get(handles.panel_tasks,'position');
handles.pos_panel_workers = get(handles.panel_workers,'position');

handles.pos_panel_showworkers = get(handles.panel_showworkers,'position');
handles.pos_panel_display = get(handles.panel_display,'position');
handles.pos_panel_worker = get(handles.panel_worker,'position');

handles.pos_ax_workers = get(handles.ax_workers,'position');
handles.pos_panel_axbox = get(handles.panel_axbox,'position');
handles.pos_sld_Workers = get(handles.sld_Workers,'position');

handles.pos_MasterUpdate = get(handles.MasterUpdate,'position');
handles.pos_but_updateNow = get(handles.but_updateNow,'position');
handles.pos_edit_status = get(handles.edit_status,'position');


% Save handles
guidata(hObject, handles);

% Find Jobmanagers and work from there...
PopulateJobManagers(hObject)







% --- Populate manager list with all managers on network -----------------
function PopulateJobManagers(hObject)

% get handles
handles = guidata(hObject);

% Set GUI update indicators and clear list boxes
set([handles.jobbox handles.taskbox],'String',{},'enable','inactive')
set(handles.jmbox,'String',{'Finding job managers...'},'value',1,'enable','inactive')
set(handles.edit_status,'string','Finding job managers...');
drawnow

% Clear any references to current jobs, and tasks
handles.cJob = [];
handles.cTask = [];
handles.Jobs = [];
handles.Tasks = [];

% Get all job managers and put names in listbox
handles.jm = findResource('jobmanager');


% sort jobmanagers into alphabetical order
jmstr = get(handles.jm,{'Name'});
[jmstr ix] = sort(jmstr);
handles.jm = handles.jm(ix);

jmstr = flipud(jmstr);
handles.jm = flipud(handles.jm);


% put names into jobmanagers listbox
set(handles.jmbox,{'value','String'},{1,jmstr});

% Get handles to manager menu items
hMenu = get(get(handles.jmbox,'uiContextMenu'),'Children');
if isempty(handles.jm)
    handles.cjm = []; % Indicate there is no current job manager
    set(hMenu,'Enable','off');
    set(handles.MenuBarJM,'Enable','off')
else
    set(hMenu,'Enable','on');
    set(handles.MenuBarJM,'Enable','on')
    handles.cjm = handles.jm(1); % Set current job manager
end
guidata(hObject,handles)

% TODO: add job manager pull down menu functionality

% Set GUI to normal
set(handles.jmbox,'enable','on')
drawnow

% Populate the list boxes and update the workers for the selected jm
PopulateJobsTasksWorkers(hObject)

% % If update timer is not already started -> start it
% if strcmp(get(handles.updatetimer,'Running'),'off')
%    start(handles.updatetimer)
% end






% --- Populates all listboxes for the first time -------------------------
function PopulateJobsTasksWorkers(hObject)

if evalin('base','exist(''debugflag'')'), disp('PopulateJobsTasksWorkers'), end

handles = guidata(hObject);

% Populate the Jobs list box
set(handles.edit_status,'string','Finding Job information...');
PopulateJobs([],[],hObject);

% Populate the Tasks list box
set(handles.edit_status,'string','Finding Task information...');
PopulateTasks([],[],hObject);

% Update workers using graphical display
set(handles.edit_status,'string','Finding Worker information...');
workers = findSelectedWorkers(handles);
plotworkermatrix(workers,handles.ax_workers)
update_currentworker(workers,handles)

% Update everything
Update([],[],hObject)

% Update status bar to say we've finished
set(handles.edit_status,'string','Idle');





% --- Populate job list --------------------------------------------------
function PopulateJobs(ignore, ignore2,hObject)

% Return if figure has already closed
if ~ishandle(hObject)
    return
end
handles = guidata(hObject);
set(handles.jobbox,'String',{'Finding jobs...'},'value',1,'enable','inactive')
drawnow

% Things may have changed from the drawnow
if ~ishandle(hObject)
    % Figure already closed
    return
end
handles = guidata(hObject);

% Initialize variables
jobstr = {};
jobloc = 1;

% Get all jobs on current job manager
handles.Jobs = get(handles.cjm,'Jobs');

% Get handles to all job menu items
hMenu = get(get(handles.jobbox,'uiContextMenu'),'Children');
if isempty(handles.Jobs)
    % No jobs -> disable menu
    handles.cJob = [];
    set(hMenu,'Enable','off');
    set(handles.MenuBarJob,'Enable','off')
else
    set(hMenu,'Enable','on');
    set(handles.MenuBarJob,'Enable','on')
    handles.cJob = handles.Jobs(jobloc);
end
handles.JobData = get(handles.Jobs,{'Name','UserName'});
state = get(handles.Jobs,{'State'});
guidata(hObject,handles)

% create a cell array of job name strings, and put in listbox
for n = 1:length(handles.Jobs)
    jobstr{n} = sprintf('%s - %s - %s',handles.JobData{n,:},state{n});
end
set(handles.jobbox,{'Value','String'},{jobloc,jobstr});
set(handles.jobbox,'enable','on')

% set context menu items to allow job to be submitted...
if isempty(state) || ~strcmp(state{jobloc},'pending')
    set(handles.MenuSubmit,'Enable','on');
else
    set(handles.MenuSubmit,'Enable','off');
end





% --- Populate task list box -----------------------------
function PopulateTasks(ignore, ignore2,hObject)

% Return if figure already closed
if ~ishandle(hObject)
    return
end

% get hold of GUI information and set update indicators
handles = guidata(hObject);
set(handles.taskbox,'String',{'Finding tasks...'},'value',1,'enable','inactive')
drawnow

% Things may have changed from the drawnow
if ~ishandle(hObject)
    % Figure already closed
    return
end
handles = guidata(hObject);

% Initialize variables
taskstr = {};
taskloc = 1;
handles.Tasks = get(handles.cJob,'Tasks'); % Get handles to all tasks

% Get handles to all task menu items
hMenu = get(get(handles.taskbox,'uiContextMenu'),'Children');
if isempty(handles.Tasks)
    % No tasks -> disable menu
    handles.cTask = [];
    set(hMenu,'Enable','off');
    set(handles.MenuBarTask,'Enable','off')
else
    set(hMenu,'Enable','on');
    set(handles.MenuBarTask,'Enable','on')
    handles.cTask = handles.Tasks(taskloc); % Specify current task
end

% Query state of the tasks
dat = get(handles.Tasks,{'ID','State'});
handles.TaskData = dat(:,1);
state = dat(:,2);

% create a cell array with the current tasks in, then put in listbox
guidata(hObject,handles)
for n = 1:length(handles.Tasks)
    taskstr{n} = sprintf('Task%d - %s',handles.TaskData{n,:},state{n});
end
set(handles.taskbox,{'value','String'},{taskloc,taskstr});

% restore normal state of GUI
set(handles.taskbox,'enable','on')








% --- Job manager changed by user ------------------------------------
function jmbox_Callback(hObject, ignore, handles)

% Bail if no job manager in listbox
if isempty(handles.jm)
    return
end

% Get currently selected manager and check if it has changed
jm = handles.jm(get(handles.jmbox,'Value'));
if ~isequal(handles.cjm,jm)

    % Yes, job manager has changed - so reassign handles.cjm
    handles.cjm = jm;

    % Stop update timer
    restarttimer = 0;
    if strcmp(handles.updatetimer.Running,'on')
        stop(handles.updatetimer)
        restarttimer = 1;
    end

    % Clear list boxes
    set([handles.jobbox handles.taskbox],'String',{})

    % Reset current worker, job, and task
    handles.cWorker = [];
    handles.cJob = [];
    handles.cTask = [];
    handles.Workers = [];
    handles.Jobs = [];
    handles.Tasks = [];

    % Set to show workers by jobmanager
    set(handles.rad_jmworkers,'value',1)

    guidata(hObject,handles);
    % Populate list boxes
    PopulateJobsTasksWorkers(hObject)

    % Restart the update timer
    if strcmp(handles.updatetimer.Running,'off')
        if restarttimer
            start(handles.updatetimer)
        end
    end
end






%--- Job changed by user -----------------------------------
function jobbox_Callback(hObject, ignore, handles)

% Bail if there are no jobs in listbox
if isempty(handles.Jobs)
    return
end

% check whether job has been changed or not
job = handles.Jobs(get(handles.jobbox,'Value'));
if ~isequal(handles.cJob,job)
    % Yes, job changed, reassign
    handles.cJob = job;
    set(handles.taskbox,'String',{})
    handles.cTask = [];
    guidata(hObject,handles);
    PopulateTasks([],[],hObject);
end

% Submit menu item is disabled if current job is not pending
state = job.state;
if strcmp(state,'pending')
    set(handles.MenuSubmit,'Enable','on');
else
    set(handles.MenuSubmit,'Enable','off')
end

% JobReport menu item is disabled if current job is pending or queued
if strcmp(state,'pending') || strcmp(state,'queued')
    set(handles.MenuSubmit,'Enable','off');
else
    set(handles.MenuSubmit,'Enable','on')
end

% Make sure everything is up to date
Update([],[],hObject)



%--- Task changed by user ------------------------------------------------
function taskbox_Callback(hObject, eventdata, handles)

% Bail if there are no tasks
if isempty(handles.Tasks)
    return
end

% Reassign current task
handles.cTask = handles.Tasks(get(handles.taskbox,'Value')); % Selected task
guidata(hObject,handles);

% Make sure everything is up to date
Update([],[],hObject)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Context Menu Item Callbacks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Job Manager menu callbacks %%%

% --- Export manager to base workspace -----------------------------------
function MenuExportJm_Callback(hObject, eventdata, handles)

jm = inputdlg('Export handle to job manager as name:','Export Job Manager');
if ~isempty(jm)
    assignin('base',jm{1},handles.cjm);
    disp(['Variable ' jm{1} ' created in the workspace.']);
end


% --- Display property inspector -----------------------------------------
function MenuPropsJM_Callback(hObject, eventdata, handles)
inspect(handles.jm(get(handles.jmbox,'Value')))




%%% Job menu callbacks %%%

% --- Export current job to base workspace -------------------------------
function MenuExportJob_Callback(hObject, eventdata, handles)

name = inputdlg('Export handle to job as name:','Export Job');
if ~isempty(name)
    assignin('base',name{1},handles.cJob);
    disp(['Variable ' name{1} ' created in the workspace.']);
end


% --- Submit current job -------------------------------------------------
function MenuSubmit_Callback(hObject, eventdata, handles)
job = handles.Jobs(get(handles.jobbox,'Value'));
submit(job)


% --- Destroy current job ------------------------------------------------
function MenuDestroy_Callback(hObject, eventdata, handles)

job = handles.Jobs(get(handles.jobbox,'Value'));
if ishandle(job)
    destroy(job)
end
Update([],[],hObject)


% --- Current Job report -------------------------------------------------
function MenuReport_Callback(hObject, eventdata, handles)

%JobReport(handles.cJob)
plottasks(handles.cjm.name,handles.cJob.name)


% --- Properties on Current Job ------------------------------------------
function MenuPropsJob_Callback(hObject, eventdata, handles)
inspect(handles.Jobs(get(handles.jobbox,'Value')))




%%% Task menu callbacks %%%

% --- Export current task to workspace -----------------------------------
function MenuExportTask_Callback(hObject, eventdata, handles)

name = inputdlg('Export handle to task as name:','Export Task');
if ~isempty(name)
    assignin('base',name{1},handles.cTask);
    disp(['Variable ' name{1} ' created in the workspace.']);
end


% --- Cancel current task ------------------------------------------------
function MenuCancel_Callback(hObject, eventdata, handles)
if ishandle(handles.cTask)
    cancel(handles.cTask)
end


% --- Destroy current task -----------------------------------------------
function MenuDestroyTask_Callback(hObject, eventdata, handles)

if ishandle(handles.cTask)
    destroy(handles.cTask)
end
Update([],[],hObject)


% --- Inspect props of current task --------------------------------------
function MenuPropsTask_Callback(hObject, eventdata, handles)
inspect(handles.Tasks(get(handles.taskbox,'Value')))





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional Menu Bar Callbacks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Close the current figure -------------------------------------------
function MenuClose_Callback(hObject, eventdata, handles)
close(handles.Figure)


% --- Run one of the plugins ---------------------------------------------
function PluginCallback(hObject,ignore)

data = get(hObject,'UserData');
handles = guidata(data{1});
fcn = data{2};
fcn(handles.cjm);


% --- Refresh the list of Job managers -----------------------------------
function MenuBarRefresh_Callback(hObject, eventdata, handles)

PopulateJobManagers(hObject,handles)






%%%%%%%%%%%%%%%%%%%%
% Update listboxes %
%%%%%%%%%%%%%%%%%%%%

% These keep the listboxes updated.  Different operations are performed to
% initially populate the listboxes using the Populate* functions.

% --- Update everything in turn ------------------------------------------
function Update(ignore,ignore2,hObject)
% Timer callback - Updates all listboxes with current info from job manager
% Each function is similar.  Each is organized to minimize the number of
% times the job manager is queried for information, which is expensive.

% This function is called either from the timer, when the checkbox is
% ticked or from the manual update button

% Exit if the live update tickbox is not ticked
handles = guidata(hObject);

if evalin('base','exist(''debugflag'')'), disp('Update'), end

set(handles.edit_status,'string','Updating Worker information...');
drawnow
UpdateWorkers(hObject);

set(handles.edit_status,'string','Updating Job information...');
drawnow
UpdateJobs(hObject);

set(handles.edit_status,'string','Updating Task information...');
drawnow
UpdateTasks(hObject);

set(handles.edit_status,'string','Idle');






% --- UPDATE WORKERS -----------------------------------------------------
function UpdateWorkers(hObject)

% Figure has been closed
if ~ishandle(hObject)
    return
end
handles = guidata(hObject);

% if get(handles.rad_jmworkers,'value') == 1
%     % Find manager's associated workers
%     workers = get(handles.cjm,{'IdleWorkers','BusyWorkers'});
%     workers = [workers{:,1}; workers{:,2}];
% else
%     % Find all workers on the network
%     workers = findResource('worker');
% end

workers = findSelectedWorkers(handles);

% call the graphical display of the workers...
plotworkermatrix(workers,handles.ax_workers)
update_currentworker(workers,handles)

% Store GUI data again
guidata(hObject,handles)





%%% --- Updates current job info from the job manager --------------------
function UpdateJobs(hObject)

% Check if figure has closed
if ~ishandle(hObject)
    return
end
handles = guidata(hObject);

% Read all jobs in current job manager
jobs = get(handles.cjm,'Jobs');

% Remove jobs that are no longer there
stillthere = ismember(handles.Jobs, jobs);
handles.Jobs(~stillthere) = [];
handles.JobData(~stillthere,:) = [];

% Ensure empty is 0-by-n rather than n-by-0
if size(handles.Jobs,2)==0
    handles.Jobs = handles.Jobs.';
end

% Ensure empty is 0-by-n rather than n-by-0
if size(handles.JobData,2)==0
    handles.JobData = handles.JobData.';
end

% Add new jobs
notnew = ismember(jobs, handles.Jobs);
handles.Jobs = [handles.Jobs; jobs(~notnew)];

% Get data about new jobs
r = get(jobs(~notnew),{'Name','UserName'}); % Only need to obtain when new
state = get(handles.Jobs,{'State'}); % Need to obtain every time
handles.JobData = [handles.JobData; r]; % Include new data in list

nJobs = length(handles.Jobs);
jobstr = cell(nJobs,1);
for n = 1:nJobs
    jobstr{n} = sprintf('%s - %s - %s',handles.JobData{n,:},state{n});
end

% Get handles to all context menu items
hMenu = get(get(handles.jobbox,'uiContextMenu'),'Children');
if nJobs == 0
    % No jobs in list
    jobloc = 1; % Listbox value
    handles.cJob = []; % Assign current job
    % Disable context and pulldown menu items
    set(hMenu,'Enable','off');
    set(handles.MenuBarJob,'Enable','off')
else
    jobloc = min(nJobs,get(handles.jobbox,'Value')); % Listbox value
    handles.cJob = handles.Jobs(jobloc); % Assign current job

    % Enable context and pull down menu items
    set(hMenu,'Enable','on');
    set(handles.MenuBarJob,'Enable','on')
end
set(handles.jobbox,{'Value','String'},{jobloc,jobstr});

% Disable Submit menu item if job is not pending
if isempty(state) || ~strcmp(state{jobloc},'pending')
    set(handles.MenuSubmit,'Enable','off');
    set(handles.MenuBarSubmit,'Enable','off');
else
    set(handles.MenuSubmit,'Enable','on');
    set(handles.MenuBarSubmit,'Enable','on');
end

% Store handles again
guidata(hObject,handles)




%--- Update current list of tasks ----------------------------------------
function UpdateTasks(hObject)

% Check if figure has closed
if ~ishandle(hObject)
    return
end
handles = guidata(hObject);

% Find all tasks associated with current job
tasks = get(handles.cJob,'Tasks');

% Remove tasks from list that are no longer there
stillthere = ismember(handles.Tasks, tasks);
handles.Tasks(~stillthere) = [];
handles.TaskData(~stillthere) = [];

% Ensure empty is 0-by-n rather than n-by-0
if size(handles.Tasks,2)==0
    handles.Tasks = handles.Tasks.';
    handles.TaskData = handles.TaskData.';
end

% Add new tasks to list
notnew = ismember(tasks, handles.Tasks);
handles.Tasks = [handles.Tasks; tasks(~notnew)];

ID = get(tasks(~notnew),{'ID'}); % Only need to obtain when new
state = get(handles.Tasks,{'State'}); % Need to obtain every time
handles.TaskData = [handles.TaskData; ID];
if length(handles.TaskData) == 1
    if ~iscell(handles.TaskData)
        handles.TaskData = {handles.TaskData};
    end
end
nTasks = length(handles.Tasks);
taskstr = {};
for n = 1:length(handles.Tasks)
    taskstr{n} = sprintf('Task%d - %s',handles.TaskData{n},state{n});
end
hMenu = get(get(handles.taskbox,'uiContextMenu'),'Children');
if nTasks == 0
    % No tasks in list
    taskloc = 1; % Listbox value
    handles.cTask = []; % Assign current task
    set(hMenu,'Enable','off');  % Disable context menu items
    set(handles.MenuBarTask,'Enable','off')
else
    taskloc = min(nTasks,get(handles.taskbox,'Value')); % Listbox value
    handles.cTask = handles.Tasks(taskloc); % Assign current task
    set(hMenu,'Enable','on'); % Enable context menu items
    set(handles.MenuBarTask,'Enable','on')
end
set(handles.taskbox,{'value','String'},{taskloc,taskstr});

% Store handles again
guidata(hObject,handles)






%%%%%%%%%%%%%%%%%%%%%%%
% Uicontrol Callbacks %
%%%%%%%%%%%%%%%%%%%%%%%

% --- Callback for when checkbox clicked ---------------------------------
function MasterUpdate_Callback(hObject, eventdata, handles)

% Disable the manual update button if automatic updating is enabled and
% vice versa.
v = get(hObject,'Value');
if v
    %turn on automatic updating
    start(handles.updatetimer)
    set(handles.but_updateNow,'Enable','off')
else
    %turn off automatic updating
    stop(handles.updatetimer)
    set(handles.but_updateNow,'Enable','on')
end




% --- Change the rate at which the interface is refreshed ----------------
function RateSlider_Callback(hObject, eventdata, handles)

% Stop timer if running
state = strcmp(handles.updatetimer.Running,'on');
if state
    stop(handles.updatetimer)
end

% Change timer period - Want smaller values on right of slider
mod = get(hObject,'Max') + get(hObject,'Min');
SecPerUpdate =  mod - get(hObject,'Value');
set(handles.RateText,'String',sprintf('%0.1f',SecPerUpdate))
oldwarn = warning('off','MATLAB:TIMER:RATEPRECISION');
set(handles.updatetimer,'Period', SecPerUpdate);
warning(oldwarn);
%drawnow

% start the timer again...
if state
    start(handles.updatetimer)
end





% --- Display some help about the program --------------------------------
function MenuBarHelp_Callback(hObject, eventdata, handles)

% read in help information from file
fid = fopen('dct_viewer_help.txt');
str = fread(fid,inf,'*char').';
fclose(fid);

% display help information in another figure
f = figure('menubar','none','NumberTitle','off','Name','About DCT Viewer');
h=uicontrol('style','edit','parent',f,'min',0,'max',2,'enable','inactive'...
    ,'units','normalized','position',[0 0 1 1],'string',str,...
    'horizontalalignment','left','fontsize',10);




% --- Executes during object creation, after setting all properties ------
function sld_Workers_CreateFcn(hObject, eventdata, handles)

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes during object creation, after setting all properties ------
function edit_status_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%%%%%%%%%%%%%%%%%%%%%%%
% Worker Display Code %
%%%%%%%%%%%%%%%%%%%%%%%

%TODO: sort the workers into some reliable order (eg. by name)


% --- function to plot a matrix of workers -------------------------------
function plotworkermatrix(workers,hAxis)

if evalin('base','exist(''debugflag'')'), disp('plotworkermatrix'), end

% get hold of persistent information about the number of workers plotted
hObject = get(hAxis,'parent');
handles = guidata(hObject);

% init the worker plot if the number of workers has changed, else update
n = size(workers,1);
n_store = handles.numWorkersPlotted;
if n ~= n_store

    % Draw patches representing workers on an axis
    cla(hAxis)
    if n == 0
        % there are no workers, just clear the axis and bail
        whandles = 0;
        set(handles.sld_Workers,'enable','off')

    else
        % there are workers...
        % work out the side of the worker drawing, set axes & worker colour
        c = [0.9176 0.9255 0.8353];
        worker_pix_zl1 = [190 140];
        worker_pix_zl2 = worker_pix_zl1*0.5;

        if get(handles.rad_zoom1,'value') == 1
            worker_pix_zoom = worker_pix_zl1;
        else
            worker_pix_zoom = worker_pix_zl2;
        end

        pos = get(hAxis,'position');
        nx = floor(pos(3)/worker_pix_zoom(1));
        ny = floor(pos(4)/worker_pix_zoom(2));

        % plot a bunch of patches, with text
        for p = 1:n

            % work out co-ordinates of worker in list
            j = mod(p-1,nx);
            i = floor((p-1)/nx);

            % draw patches to represent worker
            hPworker1(p) = patch([0.20 0.80 0.80 0.20]+j,[0.65 0.65 0.05 0.05]+i,c,'parent',hAxis,'ButtonDownFcn',{@worker_click_callback},'tag',['wp' num2str(p)]);
            hPworker2(p) = patch([0.05 0.95 0.95 0.05]+j,[0.95 0.95 0.70 0.70]+i,c,'parent',hAxis,'ButtonDownFcn',{@worker_click_callback},'tag',['wp' num2str(p)]);
            hPworker3(p) = patch([0.25 0.75 0.75 0.25]+j,[0.60 0.60 0.10 0.10]+i,'w','parent',hAxis,'ButtonDownFcn',{@worker_click_callback},'tag',['wp' num2str(p)]);

            % add text to the worker pictures if at biggest zoom level
            if get(handles.rad_zoom1,'value') == 1
                hTname(p) = text(0.5+j,0.825+i,'Name of Worker','horizontalalignment','center','parent',hAxis,'interpreter','none','fontsize',8,'ButtonDownFcn',{@worker_click_callback},'tag',['wp' num2str(p)]);
                hTstate(p) = text(0.5+j,0.28+i,'Idle','horizontalalignment','center','parent',hAxis,'interpreter','none','fontsize',10,'ButtonDownFcn',{@worker_click_callback},'tag',['wp' num2str(p)]);
                hTstate2(p) = text(0.5+j,0.46+i,'No Task','horizontalalignment','center','parent',hAxis,'interpreter','none','fontsize',8,'ButtonDownFcn',{@worker_click_callback},'tag',['wp' num2str(p)]);
            end
        end

        ylim = [0 pos(4)/pos(3)*nx*1.3];
        set(hAxis,'xlim',[0 nx],'ylim',ylim);
        handles.workerplot_ylim = ylim;

        if ylim(2)>(i+1)
            % all the workers fit, dont need scrollbar
            set(handles.sld_Workers,'enable','off');

        else
            % need to set up and turn on the scroll bar
            set(handles.sld_Workers,'enable','on')
            set(handles.sld_Workers,'min',ylim(2)-i-1)
            set(handles.sld_Workers,'sliderstep',[0.1 ylim(2)/(i+1)])
        end

        % return worker object handles
        if get(handles.rad_zoom1,'value') == 1
            whandles = {hPworker1 hPworker2 hPworker3 hTname hTstate hTstate2};
        else
            whandles = {hPworker1 hPworker2 hPworker3};
        end

    end

    handles.whandles = whandles;
    handles.numWorkersPlotted = n;
    handles.selectedworker = -1;

else
    % The number of workers hasn't changed
    whandles = handles.whandles;
end

% store data back to GUI again
guidata(hObject,handles);


% Update the information on the worker plot if necessary
if n > 0
    % collect information from workers, in one go
    workerData = get(workers,{'name','currentjob','currenttask','state'});

    for p = 1:n
        winfo(p).name = workerData(p,1);

        % get information about current job, if there is one
        if get(handles.rad_zoom1,'value') == 1
            if size(workerData{p,2}) == [0 0]
                winfo(p).jobname = 'Idle';
            else
                winfo(p).jobname = get(workerData{p,2},'name');
            end
        end

        % get information about the current task, if there is one
        if size(workerData{p,3}) == [0 0]
            winfo(p).taskname = 'No Task';
            winfo(p).facecolor = 'w';
        else
            if get(handles.rad_zoom1,'value') == 1
                winfo(p).taskname = ['Task ' num2str(get(workerData{p,3},'id'))];
            end
            winfo(p).facecolor = 'g';
        end
    end

    % update the plot, in one go
    for p = 1:n
        if get(handles.rad_zoom1,'value') == 1
            set(whandles{4}(p),'string',winfo(p).name)
            set(whandles{5}(p),'string',winfo(p).jobname);
            set(whandles{6}(p),'string',winfo(p).taskname);
        end
        set(whandles{3}(p),'facecolor',winfo(p).facecolor);
        if strcmp(workerData(p,4),'unavailable')
            set(whandles{1}(p),'facecolor','r')
        else
            set(whandles{1}(p),'facecolor',[0.9176 0.9255 0.8353])
        end
    end
end

%if n > 0
%    update_workerplot(workers,whandles,hObject)
%end


% store data back to GUI again
guidata(hObject,handles);







% --- Executes on mouse press over one of the 'worker patches' -----------
function worker_click_callback(wpatch, eventdata)

if evalin('base','exist(''debugflag'')'), disp('worker_click_callback'), end

% get hold of the GUI handles information
handles = guidata(wpatch);
whandles = handles.whandles;

% work out which worker has been selected
wtag = get(wpatch,'tag');
wnum = str2num(wtag(3:end));

% if we have just selected the same worker, de-select it, else highlight it
if handles.selectedworker == wnum
    %de-select worker
    set(whandles{1}(handles.selectedworker),'linewidth',0.5);
    set(whandles{2}(handles.selectedworker),'linewidth',0.5);
    set(whandles{3}(handles.selectedworker),'linewidth',0.5);
    handles.selectedworker = -1;
else
    % deselect the previous worker, providing one has been selected
    if handles.selectedworker ~= -1
        set(whandles{1}(handles.selectedworker),'linewidth',0.5);
        set(whandles{2}(handles.selectedworker),'linewidth',0.5);
        set(whandles{3}(handles.selectedworker),'linewidth',0.5);
    end

    % highlight the newly selected worker
    %set(whandles{1}(wnum),'linewidth',2);
    %set(whandles{2}(wnum),'linewidth',2);
    %set(whandles{3}(wnum),'linewidth',2);
    handles.selectedworker = wnum;

    % put the relevant information in the text display

end

% store data back to the GUI
guidata(wpatch,handles);

% call to update the selected worker details
workers = findSelectedWorkers(handles);
update_currentworker(workers,handles)








% --- Update the information in the current worker details panel ---------
function update_currentworker(workers,handles)

% check that there are workers to plot!
if size(workers,1) == 0
    h = [handles.txt_selworkername handles.txt_selworkerhost handles.txt_selworkerstatus handles.txt_selworkerjob handles.txt_selworkertask];
    set(h,'string','')
    return
end

if handles.selectedworker == -1
    h = [handles.txt_selworkername handles.txt_selworkerhost handles.txt_selworkerstatus handles.txt_selworkerjob handles.txt_selworkertask];
    set(h,'string','')
else
    if(handles.selectedworker <= length(workers))
        w = workers(handles.selectedworker);
        wstate = get(w,{'name','hostname','state','currentjob','currenttask'});
        set(handles.txt_selworkername,'string',wstate{1})
        set(handles.txt_selworkerhost,'string',wstate{2})

        if size(wstate{4}) == [0 0]
            set(handles.txt_selworkerstatus,'string',['Idle (' wstate{3} ')'])
            set([handles.txt_selworkerjob handles.txt_selworkertask],'string','')
        else
            set(handles.txt_selworkerstatus,'string',['Busy (' wstate{3} ')'])
            set(handles.txt_selworkerjob,'string',get(wstate{4},'name'))
            set(handles.txt_selworkertask,'string',num2str(get(wstate{5},'id')))
        end
    end

end







% --- return the workers that we are interested in -----------------------
function workers = findSelectedWorkers(handles)

% call to update the selected worker details
if get(handles.rad_jmworkers,'value') == 1
    % Find manager's associated workers
    workers = get(handles.cjm,{'IdleWorkers','BusyWorkers'});
    workers = [workers{:,1}; workers{:,2}];
else
    % Find all workers on the network
    workers = findResource('worker');
end

% sort the workers into order
names = get(workers,{'name'});
[sorted ix] = sort(names);
workers = workers(ix);






% --- Executes on button press in rad_jmworkers or rad_allworkers --------
function rad_workers_Callback(hObject, eventdata, handles)

if evalin('base','exist(''debugflag'')'), disp('rad_workers_callback'), end

% prevent the radio button from being de-selected
if get(hObject,'value') == 0
    set(hObject,'value',1)
end

% TODO: redraw the workers
hAxis = handles.ax_workers;

% update the plot and the current worker details
workers = findSelectedWorkers(handles);
plotworkermatrix(workers,hAxis)
update_currentworker(workers,handles)








% --- Executes when Figure is resized ------------------------------------
function Figure_ResizeFcn(hObject, eventdata, handles)

% Resize the figure, based on the initial coordinates stored in handles
% during the figure_opening function.

if evalin('base','exist(''debugflag'')'), disp('Resize'), end

% Work out how much the figure has been resized by
pos = get(hObject,'position');
try
    deltaY = pos(4) - handles.form_initheight;
    deltaX = pos(3) - handles.form_initwidth;
catch
    % e.g. if the figure opening function has not been called yet
    % (can occur on small screen resolutions) & handles.form_initheight
    % does not exist
    return
end

d = get(handles.Figure,'position');
handles.form_initwidth = d(3);
handles.form_initheight = d(4);

% If the figure is smaller than the initial (minimum) size, resize it back
% to the original size
resizefig = 0;
if deltaY<0
    pos(4) = handles.form_initheight;
    pos(2) = pos(2)+deltaY;
    resizefig = 1;
end
if deltaX<0
    pos(3) = handles.form_initwidth;
    resizefig = 1;
end
if resizefig == 1
    set(hObject,'position',pos);
end

% Ensure deltaX and Y ae positive
deltaY = (deltaY>0)*deltaY;
deltaX = (deltaX>0)*deltaX;

% Keep top 3 panels at the top of the figure
set(handles.panel_jms,'position',handles.pos_panel_jms + [0 deltaY 0 0]);
set(handles.panel_jobs,'position',handles.pos_panel_jobs + [0 deltaY 0 0]);
set(handles.panel_tasks,'position',handles.pos_panel_tasks + [0 deltaY 0 0]);

% Move LHS panels up to match
set(handles.panel_showworkers,'position',handles.pos_panel_showworkers + [0 deltaY 0 0]);
set(handles.panel_display,'position',handles.pos_panel_display + [0 deltaY 0 0]);
set(handles.panel_worker,'position',handles.pos_panel_worker + [0 deltaY 0 0]);

% Stretch large panel at the bottom and worker plot
set(handles.panel_workers,'position',handles.pos_panel_workers + [0 0 deltaX deltaY]);
set(handles.ax_workers,'position',handles.pos_ax_workers + [0 0 deltaX deltaY]);
set(handles.panel_axbox,'position',handles.pos_panel_axbox + [0 0 deltaX deltaY]);
set(handles.sld_Workers,'position',handles.pos_sld_Workers + [deltaX 0 0 deltaY]);

% Move and status bar stuff
set(handles.but_updateNow,'position',handles.pos_but_updateNow + [deltaX 0 0 0]);
set(handles.MasterUpdate,'position',handles.pos_MasterUpdate + [deltaX 0 0 0]);
set(handles.edit_status,'position',handles.pos_edit_status + [0 0 deltaX 0]);

% Force a redraw of all the workers, so we get sensible sizes
handles.numWorkersPlotted = -1;
guidata(hObject,handles);
UpdateWorkers(hObject)




% --- Callback for when update button pressed ----------------------------
function but_updateNow_Callback(hObject, eventdata, handles)

% Call to update everything
Update([],[],hObject)





% --- Executes on slider movement ----------------------------------------
function sld_Workers_Callback(hObject, eventdata, handles)

val = get(hObject,'value');
%origval = get(handles.ax_workers,'ylim')
origval = [0 1];
ylim_orig = handles.workerplot_ylim;

set(handles.ax_workers,'ylim',ylim_orig-val)





% --- Executes on button press in rad_zoom -------------------------------
function rad_zoom_Callback(hObject, eventdata, handles)

%TODO: fix this so that only 1 radio button can ever be lit

handles.numWorkersPlotted = -1;
guidata(hObject,handles);
UpdateWorkers(hObject)




