function varargout = tom_browse(varargin)
%TOM_BROWSE opens the new TOM browser.
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%     TOM_EMBROWSE
%
%   created by SN 07/23/09
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
gui_State = struct('gui_Name',          mfilename, ...
    'gui_Singleton',     gui_Singleton, ...
    'gui_OpeningFcn',    @tom_browse_OpeningFcn, ...
    'gui_OutputFcn',     @tom_browse_OutputFcn, ...
    'gui_LayoutFcn',     [], ...
    'gui_Callback',      []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before tom_browse is made visible.
function tom_browse_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for tom_browse
handles.output = hObject;
guidata(hObject, handles);

if nargin == 3,
    initial_dir = pwd;
    handles.file_filter='*';
    
elseif nargin > 4
    if strcmpi(varargin{1},'dir')
        if exist(varargin{2},'dir')
            initial_dir = varargin{2};
        else
            errordlg('Input argument must be a valid directory','Input Argument Error!')
            return
        end
    else
        errordlg('Unrecognized input argument','Input Argument Error!');
        return;
    end
end
%handles.ha = axes('Units','Pixels','Position',[50,60,200,185]);

handles.current_dir = [pwd filesep];
set(handles.tom_browse_figure,'Name',handles.current_dir);


% Populate the listbox

handles=load_listbox(hObject,initial_dir,handles);
handles.Position_orig=get(handles.display,'Position');
axis(handles.display,'off');
guidata(hObject, handles);

project_directory_list=which('tom_browse_project_dirs.mat','-all');
if isempty(project_directory_list)
    disp('cannot find project directory list.');
else
    disp('found project directory list(s):');
    disp(char(project_directory_list))
    disp('first file in list is used.');
end;
show_info(hObject, handles);
% Return figure handle as first output argument

% UIWAIT makes tom_browse wait for user response (see UIRESUME)
% uiwait(handles.tom_browse_figure);


% --- Outputs from this function are returned to the command line.
function varargout = tom_browse_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ------------------------------------------------------------
% Callback for list box - open .fig with guide, otherwise use open
% ------------------------------------------------------------
function tom_browse_listbox_files_Callback(hObject, eventdata, handles)
% Read the current directory and sort the names
index_selected = get(handles.tom_browse_listbox_files,'Value');
file_list = get(handles.tom_browse_listbox_files,'String');
filename = file_list{index_selected};
[path,name,ext] = fileparts(filename);
ver='';
set(handles.dims,'String','-');
if isequal(filename,['.' filesep]) || isequal(filename,'.'); handles=load_listbox(hObject,handles.current_dir,handles);return; end;
if (isequal(filename,['..' filesep]) || isequal(filename,'..')) && strcmp(get(handles.tom_browse_figure,'SelectionType'),'open')
    w=findstr(handles.current_dir,filesep);
    handles.current_dir=handles.current_dir(1:w(size(w,2)-1));
    %[path,name,ext,ver] = fileparts(handles.current_dir)
    handles=load_listbox(hObject,handles.current_dir,handles);
    set(handles.tom_browse_figure,'Name',handles.current_dir);
    set(handles.tom_browse_path,'String',handles.current_dir);
    guidata(hObject, handles);
    return;
end;
if handles.is_dir(index_selected) && strcmp(get(handles.tom_browse_figure,'SelectionType'),'open')    
%    if handles.is_dir(index_selected);
        change_path=filename;
        handles.current_dir=[handles.current_dir change_path];
        handles=load_listbox(hObject,handles.current_dir,handles);
        set(handles.tom_browse_figure,'Name',handles.current_dir);
        set(handles.tom_browse_path,'String',handles.current_dir);
        guidata(hObject, handles);
        return;
%    end;
end;
%        if handles.is_dir(index_selected); cd(filename); handles=load_listbox(hObject,'.',handles); return; end;
filename=[handles.current_dir filename];
[type]=tom_determine_file_type(filename);
handles.file_byte_M=handles.file_bytes_M{index_selected};
set(handles.dims,'String',[num2str(handles.file_byte_M) ' kB']);
if handles.is_dir(index_selected)
    type='directory';
    set(handles.dims,'String','-');
end;    
handles.filename=filename;
[success,message,messageid]=fileattrib(filename);
att='';
if success==0
    disp(messageid);
    att='-';
else
    handles.fileattribute=interpret_attribute(message);
end;
handles.filedate=handles.file_dates(index_selected);
handles.filetype=type;
set(handles.tom_browse_figure,'Name',handles.current_dir);
set(handles.tom_browse_path,'String',handles.current_dir);
set(handles.filetype_text,'String',handles.filetype);
set(handles.file_modified,'String',handles.filedate);
set(handles.file_attribute,'String',handles.fileattribute);
set(handles.display,'Position',handles.Position_orig);
guidata(hObject, handles);
if strcmp(get(handles.tom_browse_figure,'SelectionType'),'open')
    do_it(handles,filename, type);
else if strcmp(get(handles.tom_browse_figure,'SelectionType'),'normal')
        if ~isempty(imformats(type)) read_as_image(handles,filename);
        end;
        switch type
            case 'unknown'; %read_as_mat(handles,filename);
                show_info(hObject, handles);
                return;
            case 'em'; read_as_em(handles,filename); return;
            case 'mrc'; read_as_mrc(handles,filename); return;
            case 'imagic'; read_as_imagic(handles,filename); return;
            case 'spider'; read_as_spider(handles,filename); return;
            case 'MAT-File';  read_as_mat(handles,filename); return;
            case 'xmipp_fsc';  read_as_xmipp_fsc(handles,filename); return;
            otherwise;
                
        end
    end;
end;

guidata(hObject, handles);

function hand=load_listbox(hObject,dir_path,handles)
% Read the current directory and sort the names
dir_struct = dir([dir_path filesep handles.file_filter]);
if findstr(handles.file_filter,'.'); dir_struct_root=dir([dir_path filesep '*.' ]); dir_struct(size(dir_struct,1)+1)=dir_struct_root(1); dir_struct(size(dir_struct,1)+1)=dir_struct_root(2); end;
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
handles.file_names = sorted_names;
handles.is_dir = [dir_struct.isdir];
handles.sorted_index = sorted_index;
handles.is_dir = handles.is_dir(sorted_index);
for i=1:size(handles.is_dir,2);
    handles.file_dates{i} = dir_struct(sorted_index(i)).date;
    handles.file_bytes_M{i} = round(10.*dir_struct(sorted_index(i)).bytes./(1024))./10; %in kB
end;
for i=1:size(handles.is_dir,2); if  handles.is_dir(i); handles.file_names(i,:)=strcat(handles.file_names(i), filesep);  end;end;
guidata(handles.tom_browse_figure,handles);
set(handles.tom_browse_listbox_files,'String',handles.file_names,...
    'Value',1);
set(handles.tom_browse_path,'String',handles.current_dir);

hand=handles;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tom_browse_listbox_files_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','black');


% --- Executes during object creation, after setting all properties.
function tom_browse_figure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tom_browse_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Add the current directory to the path, as the pwd might change thru' the
% gui. Remove the directory from the path when gui is closed
% (See figure1_DeleteFcn)
setappdata(hObject, 'StartPath', pwd);
addpath(pwd);


% --- Executes during object deletion, before destroying properties.
function tom_browse_figure_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to tom_browse_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Remove the directory added to the path in the figure1_CreateFcn.
if isappdata(hObject, 'StartPath')
    rmpath(getappdata(hObject, 'StartPath'));
end


function tom_browse_filefilter_Callback(hObject, eventdata, handles)
% filefilter in the center bottom of the gui
contents = get(hObject,'String');
file_filter=contents{get(hObject,'Value')};
if isequal(file_filter,'*.*');file_filter='*'; end;
handles.file_filter=file_filter;
guidata(hObject, handles);
handles=load_listbox(hObject,pwd,handles);
guidata(hObject, handles);

function tom_browse_filefilter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function browse_dir_Callback(hObject, eventdata, handles)
% browse button 'Dir' in the center of the gui to change the current dir.
p=uigetdir;
if p~=0
    handles.current_dir=[p filesep];
    handles=load_listbox(hObject,handles.current_dir,handles);
    set(handles.tom_browse_figure,'Name',handles.current_dir);
    set(handles.tom_browse_path,'String',handles.current_dir);
    guidata(hObject, handles);
end;

function matlab_variable_Callback(hObject, eventdata, handles)
% load file to matlab workspace
filename=char(handles.file_names(get(handles.tom_browse_listbox_files,'Value')));
filename=[handles.current_dir filename];
[type]=tom_determine_file_type(filename);
switch type
    case 'unknown'; disp('type unknown');
    case 'em'; in=tom_emreadc(filename);
    case 'mrc'; in=tom_mrcread(filename);
    case 'imagic'; in=tom_imagicread(filename);
    case 'spider'; in=tom_spiderread(filename);
    case 'MAT-File'; s=load(filename);
        var_base = cellstr(evalin('base','who'));
        if isempty(var_base); var_base=1; end;
        for i=1:size(s,1);
            for ii=1:size(var_base,1)
                if isequal(char(fieldnames(s(i))),char(var_base(ii)))
                    ButtonName = questdlg('Variable exists! Overwrite variable?', 'Yes');
                    if isequal(ButtonName,'Yes')
                        assignin('base',char(fieldnames(s(i))),s(i));
                    else
                        break;
                    end;
                    
                else;
                    assignin('base',char(fieldnames(s(i))),s(i));
                end;
            end;
        end;
        return;
    otherwise in=imread(filename); type='img';
end
var_base = cellstr(evalin('base','who'));
if isequal(type,'img')
    var_exp='img_';%image
else
    if in.Header.Size(3)>1;var_exp='vol_';%volume
    else;var_exp='img_';%image
    end;
end;
if isempty(var_base)
    var1=[var_exp '1'];
else
    wos=0;ii=1;
    for i=1:size(var_base,1)
        wo=findstr(char(sort(var_base(i))),var_exp);
        if ~isempty(wo) wos=ii; ii=ii+1; else; end;
    end;
    var1=[var_exp num2str(wos+1)];
end;
assignin('base',var1,in);

function tom_browse_path_Callback(hObject, eventdata, handles)
% edit text field for directory, upper left corner
s=get(hObject,'String');
if exist(s)
    if ~isequal(s(end), filesep); s=[s filesep]; end;
    handles.current_dir=s;
    handles=load_listbox(hObject,handles.current_dir,handles);
    guidata(hObject, handles);
else
    set(hObject,'String',handles.current_dir);
end;

% --- Executes during object creation, after setting all properties.
function tom_browse_path_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function populate_header_info(handles);
% updates the header informations on gui
set(handles.angle,'String','-');
set(handles.defocus,'String','-');
set(handles.objpixsize,'String','-');
set(handles.comment,'String','-');
set(handles.filetype_text,'String',handles.filetype);
set(handles.dims,'String',[num2str(handles.Header.Size(1)) 'x' num2str(handles.Header.Size(2)) 'x' num2str(handles.Header.Size(3)) ' px']);
if isequal(handles.filetype,'directory')
    set(handles.dims,'String','-');
end;
%if isequal(handles.filetype,'unknown')
%    set(handles.dims,'String',num2str(file_byte_M));
%end;
if isequal(handles.filetype,'em')
    set(handles.angle,'String',[num2str(handles.Header.Tiltangle)]);
    set(handles.defocus,'String',[num2str(handles.Header.Defocus./10000)]);
    set(handles.objpixsize,'String',[num2str(handles.Header.Objectpixelsize)]);
    set(handles.comment,'String',[num2str(char(handles.Header.Comment))]);
end;

function do_it(handles,filename,type);

if isequal(type,'Matlab file') || isequal(type,'C file') ...
        || isequal(type,'Text file') || isequal(type,'C++ file') ...
        || isequal(type,'Header file')
    edit(filename);
end;

header=read_header(filename,type);
if isempty(header); return; end;
% check if volume or image data
if header.Size(3)==1
    actions{1}='Show in Imtool';
    actions{2}='Show in Imtool 1x binning';
    actions{3}='Show in Imtool 2x binning';
else
    actions{1}='Show in tom_volxyz';
    actions{2}='Show in tom_volxyz 1x binning';
    actions{3}='Show in tom_volxyz 2x binning';
end;

[s,v] = listdlg('PromptString','Select an action:',...
    'SelectionMode','single',...
    'ListString',actions);
% cancel - do nothing
if v==0 return; end;

do=actions{s};
if header.Size(3)==1; dim=2; else dim=3; end;
do_action(do,filename,type,dim);

function header=read_header(filename,type)
header='';
switch type
    case 'em'
        tmp=tom_reademheader(filename);
        header=tmp.Header;
    case 'spider'
        tmp=tom_readspiderheader(filename);
        header=tmp.Header;
end;
% other known file formats ,.jpg, .tif, etc.
if isempty(header)
    try
        tmp=imfinfo(filename);
        header.Size(1)=tmp.Width;
        header.Size(2)=tmp.Height;
        header.Size(3)=1;
    catch
        disp('cannot read this as an image file.');
    end;
end;

function do_action(action,filename,type,dim)

binning_factor=0;
in='';

if findstr(action,'binning');
    if findstr(action,'x');
        binning_factor=str2num(action(findstr(action,'x')-1));
    end;
end;
switch type
    case 'em'
        if dim==2
            in=tom_emreadc2(filename,'binning',[binning_factor binning_factor 1]);
        else
            in=tom_emreadc2(filename,'binning',[binning_factor binning_factor binning_factor]);
        end;
        in=in.Value;
        
    case 'spider'
        in=tom_spiderread(filename);
        in=tom_bin(in.Value,binning_factor);
    case 'image'
        in=imread(filename);
        in=imresize(in,1./(2.^binning_factor));
end;

% all the other known formats
if isempty(in)
    try
        in=imread(filename);
        in=imresize(in,1./(2.^binning_factor));
    catch
        disp('cannot read this file. Try tom_rawread.');
    end;
end;

m=mean2(in(1:4:end,1:4:end,:));
s=std2(in(1:4:end,1:4:end,:));

if findstr(action,'Show')
    if dim==2
        h=imtool(in','DisplayRange',[m-3.*s m+3.*s],'InitialMagnification','fit');
        set(h,'Name',filename);
    else
        tom_volxyz(in);
    end;
end;

function  show_info(hObject, handles)
% show the help on the right side in red
handles.image=imshow(zeros(256,256),'Parent',handles.display);
t0=text(20,20,'tom\_browse','Parent',handles.display);
set(t0,'Color','red');
set(t0,'Fontsize',16);
t1=text(30,40,'- Browse with mouse pointer or up&down key.','Parent',handles.display);
set(t1,'Color','red');
t2=text(30,60,'- Single click for display.','Parent',handles.display);
set(t2,'Color','red');
t3=text(30,80,'- Double click or ''enter'' for more options.','Parent',handles.display);
set(t3,'Color','red');
t4=text(30,100,'- Known formats: EM, Spider, Xmipp.','Parent',handles.display);
set(t4,'Color','red');
t5=text(30,120,'- Known formats: Matlab, Images(.tif,.jpg,...),C,C++,.txt','Parent',handles.display);
set(t5,'Color','red');
t6=text(30,140,'- Move mouse pointer over buttons to get help.','Parent',handles.display);
set(t6,'Color','red');
t7=text(30,160,'- The projects directory file: tom\_browse\_project\_dirs.mat','Parent',handles.display);
set(t7,'Color','red');
t8=text(30,170,'  must be saved somewhere in your Matlab path.','Parent',handles.display);
set(t8,'Color','red');
guidata(hObject, handles);




function add_button_Callback(hObject, eventdata, handles)
% add a project directory to the list
if exist('tom_browse_project_dirs.mat')
    load('tom_browse_project_dirs.mat');
    tom_browse_project_dirs_path=which('tom_browse_project_dirs.mat');
else
    % create new one
    ButtonName = questdlg('Create new one?', ...
        'Cannot find project directory list.', ...
        'yes', 'no', 'yes');
    switch ButtonName,
        case 'yes',
            project_dirs{1}=handles.current_dir;
            directoryname = uigetdir;
            save([directoryname filesep 'tom_browse_project_dirs.mat'],'project_dirs');
            disp('Project directory list created.');
        case 'no',
            disp('No project directory list created.'); return;
    end % switch
    return;
end;
project_dirs{size(project_dirs,2)+1}=handles.current_dir;
save(tom_browse_project_dirs_path,'project_dirs');
disp('Directory added to project directory list.');

function go_button_Callback(hObject, eventdata, handles)
% change to a project directory from the list
if exist('tom_browse_project_dirs.mat')
    load('tom_browse_project_dirs.mat');
else
    % create new one
    ButtonName = questdlg('Create new one?', ...
        'Cannot find project directory list.', ...
        'yes', 'no', 'yes');
    switch ButtonName,
        case 'yes',
            project_dirs{1}=handles.current_dir;
            directoryname = uigetdir;
            save([directoryname filesep 'tom_browse_project_dirs.mat'],'project_dirs');
            disp('Project directory list created.');
        case 'no',
            disp('No project directory list created.'); return;
    end 
    return;
end;

if isempty(project_dirs)
    ButtonName = questdlg('Project directory list is empty. Add a dir first.', ...
        'Project directory list.', ...
        'ok', 'ok');
else
    [s,v] = listdlg('PromptString','Change to directory:',...
        'Name','Go to directory',...
        'SelectionMode','single',...
        'ListString',project_dirs);
    if v==1
        cd(project_dirs{s});
        handles.current_dir = [pwd filesep];
        set(handles.tom_browse_figure,'Name',handles.current_dir);
        % Populate the listbox
        handles=load_listbox(hObject,handles.current_dir,handles);
    else return;
    end;
end;
guidata(hObject, handles);

function delete_button_Callback(hObject, eventdata, handles)
% delete a project directory from the list
if exist('tom_browse_project_dirs.mat')
    load('tom_browse_project_dirs.mat');
    tom_browse_project_dirs_path=which('tom_browse_project_dirs.mat');
else
    disp('create new one');return;
end;

[s,v] = listdlg('PromptString','Delete directory from list:',...
    'Name','Delete directory from list',...
    'SelectionMode','single',...
    'ListString',project_dirs);
if v==1
    project_dirs{s}=project_dirs{size(project_dirs,2)};
    tmp=project_dirs(1:size(project_dirs,2)-1);
    project_dirs=tmp;
else return;
end;
save(tom_browse_project_dirs_path,'project_dirs');
disp('Directory deleted from project directory list.');

function change_matlab_dir_button_Callback(hObject, eventdata, handles)
% change matlab path to current parh
cd(handles.current_dir);

function read_as_em(handles,filename)
in=tom_reademheader(filename);
r=ceil(in.Header.Size(1)/256);
if in.Header.Size(3)>1
    t=tom_emreadc2(filename,'subregion',[1 1 in.Header.Size(3)./2+1],[in.Header.Size(1)-1 in.Header.Size(2)-1 0]);
    t=t.Value;
    tr=t(1:4:end,1:4:end);
    m=mean2(tr);
    s=std2(double(tr));
    handles.image=imshow(t','Parent',handles.display,'DisplayRange',[m-3.*s m+3.*s]);
else
    t=tom_emreadc2(filename,'resample',[r r 1]);
    t=t.Value;
    tr=double(t(1:4:end,1:4:end));
    m=mean2(tr);
    s=std2(tr);
    handles.image=imshow(t','Parent',handles.display,'DisplayRange',[m-3.*s m+3.*s]);
end;
handles.Header=in.Header;
populate_header_info(handles);

function read_as_em_display(handles,filename)
in=tom_emreadc2(filename);
in=in.Value;
in2=in(1:4:end,1:4:end);
m=mean2(in2);
s=std2(double(in2));
imtool(in,[m-s.*3 m+s.*3]);

function read_as_image(handles,filename)
[X,MAP] = imread(filename);
handles.image=imshow(X,'Parent',handles.display);
handles.Header.Size(1)=size(X,2);handles.Header.Size(2)=size(X,1);handles.Header.Size(3)=1;
populate_header_info(handles);

function read_as_mat(handles,filename)
handles.Header.Size(1)=1;handles.Header.Size(2)=1;handles.Header.Size(3)=1;
populate_header_info(handles);

function read_as_xmipp_fsc(handles,filename)
fsc=importdata(filename);
plot(fsc.data(:,1),fsc.data(:,2),'Linewidth',2,'Parent',handles.display); hold(handles.display,'on');
plot(fsc.data(:,1),fsc.data(:,2),'ro','MarkerSize',2,'Linewidth',2,'Parent',handles.display); hold(handles.display,'off');
set(handles.display,'Position',[80 12 60 20]);
set(handles.display,'Ylim',[0 1.1]);
t=title(handles.display,['Fourier Shell Correlation, Nyquist @ ' num2str(fsc.data(end,4)) ' Ang']);
set(t,'Color','white');
xlabel(handles.display,'Resolution [1/Ang]');
set(handles.display,'Ytick',[0:0.1:1.1]);
set(handles.display,'Xtick',[0:fsc.data(end,1)./5:fsc.data(end,1)]);
ylabel(handles.display,'FSC');
set(handles.display,'XColor','white');
set(handles.display,'YColor','white');

function read_as_image_display(handles,filename)
[X,MAP] = imread(filename);
in2=X(1:4:end,1:4:end);
m=mean2(in2);
s=std2(double(in2));
imtool(X,[m-s.*3 m+s.*3]);
handles.Header.Size(1)=size(X,1);handles.Header.Size(2)=size(X,2);handles.Header.Size(3)=1;
populate_header_info(handles);

function read_as_spider(handles,filename)
in=tom_readspiderheader(filename);
t=tom_spiderread(filename);
t=t.Value;
tr=t(1:4:end,1:4:end);
m=mean2(tr);
s=std2(double(tr));
if in.Header.Size(3)>1
    handles.image=imshow(t(:,:,size(t,3)./2+1),'Parent',handles.display,'DisplayRange',[m-3.*s m+3.*s]);
else
    handles.image=imshow(t,'Parent',handles.display,'DisplayRange',[m-3.*s m+3.*s]);
end;
handles.Header=in.Header;
populate_header_info(handles);

function edit_header_button_Callback(hObject, eventdata, handles)
% open the EM header edit gui

filename = handles.filename;
if isequal(handles.filetype,'em')
    i=tom_emreadc2(filename);
    tom_editemheader(i);
end
handles=load_listbox(hObject,handles.current_dir,handles);
guidata(hObject, handles);
populate_header_info(handles);

function att=interpret_attribute(m);
att='';
switch m.UserRead
    case 0
        att=strcat(att,'-');
    case 1
        att=strcat(att,'r');
end;
if isnan(m.UserRead)
    att=strcat(att,'u');
end;
switch m.UserWrite
    case 0
        att=strcat(att,'-');
    case 1
        att=strcat(att,'w');
end;
if isnan(m.UserWrite)
    att=strcat(att,'u');
end;
switch m.UserExecute
    case 0
        att=strcat(att,'-');
    case 1
        att=strcat(att,'x');
end;
if isnan(m.UserExecute)
    att=strcat(att,'u');
end;
att=strcat(att,',');
switch m.GroupRead
    case 0
        att=strcat(att,'-');
    case 1
        att=strcat(att,'r');
end;
if isnan(m.GroupRead)
    att=strcat(att,'u');
end;
switch m.GroupWrite
    case 0
        att=strcat(att,'-');
    case 1
        att=strcat(att,'w');
end;
if isnan(m.GroupWrite)
    att=strcat(att,'u');
end;
switch m.GroupExecute
    case 0
        att=strcat(att,'-');
    case 1
        att=strcat(att,'x');
end;
if isnan(m.GroupExecute)
    att=strcat(att,'u');
end;
att=strcat(att,',');
switch m.OtherRead
    case 0
        att=strcat(att,'-');
    case 1
        att=strcat(att,'r');
end;
if isnan(m.OtherRead)
    att=strcat(att,'u');
end;
switch m.OtherWrite
    case 0
        att=strcat(att,'-');
    case 1
        att=strcat(att,'w');
end;
if isnan(m.OtherWrite)
    att=strcat(att,'u');
end;
switch m.OtherExecute
    case 0
        att=strcat(att,'-');
    case 1
        att=strcat(att,'x');
end;
if isnan(m.OtherExecute)
    att=strcat(att,'u');
end;
