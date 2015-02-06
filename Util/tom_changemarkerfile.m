function varargout = tom_changemarkerfile(varargin)
%TOM_CHANGEMARKERFILE deletes a column in a marker file
%
%   varargout = tom_changemarkerfile(varargin)
%
%   It's an interactive tool to delete a colomn in a marker file. 
%   The user has also the possibility to modify the tilt series.
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
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_SETMARK
%
%   created by WDN 07/06/04
%   updated by WDN 05/15/06
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
                   'gui_OpeningFcn', @tom_changemarkerfile_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_changemarkerfile_OutputFcn, ...
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


% --- Executes just before tom_changemarkerfile is made visible.
function tom_changemarkerfile_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = tom_changemarkerfile_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
init;
varargout{1} = handles.output;

% ---------------------------------
% --- BUTTON BROWSE MARKER FILE ---
% ---------------------------------
function browse_MF_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile({'*.em;';'*.*'}, 'Pick up a Marker file');
if isequal(filename,0) | isequal(pathname,0)
    return;
end;
handles.PathMF=pathname;
handles.NameMF=filename;
set(handles.mfdir,'String',handles.PathMF);
set(handles.mfname,'String',handles.NameMF);
a=tom_reademheader([handles.PathMF handles.NameMF]);
dimx=num2str(a.Header.Size(1));
dimy=num2str(a.Header.Size(2));
dimz=num2str(a.Header.Size(3));
set(handles.mfsize,'String',[ dimx ' x ' dimy ' x ' dimz ]);
set(handles.mf2,'String',dimy);
set(handles.mf3,'String',dimz);
set(handles.rmcol,'Enable','on');
handles.SizeMF_org=[str2num(dimx);str2num(dimy);str2num(dimz)];
handles.LastProj=str2num(dimy);
if ~isempty(findstr(filename,'.em'))
    handles.NewNameMF=[filename(1:findstr(filename,'.em')-1) '_mod.em'];
else
    handles.NewNameMF='MF_mod.em';
end
set(handles.new_mfname,'String',handles.NewNameMF);
set(handles.mfdoit,'Enable','on');
set(handles.check_renamemf,'Enable','on');
guidata(hObject,handles);
% ------------------------------------
% --- CHECKBOX SAVE MARKER FILE AS ---
% ------------------------------------
function check_renamemf_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    set(handles.new_mfname,'Enable','on');
    if ~isempty(findstr(handles.NameMF,'.em'))
        handles.NewNameMF=[handles.NameMF(1:findstr(handles.NameMF,'.em')-1) '_mod.em'];
    else
        handles.NewNameMF='MF_mod.em';
    end
else
    set(handles.new_mfname,'Enable','off');
    handles.NewNameMF='';
end
guidata(hObject,handles);
% ------------------------------------
% --- EDIT BOX SAVE MARKER FILE AS ---
% ------------------------------------
function new_mfname_Callback(hObject, eventdata, handles)
handles.NewNameMF=get(hObject,'String');
guidata(hObject,handles);
% -----------------------------
% --- EDIT BOX IMAGE NUMBER ---
% -----------------------------
function rmcol_Callback(hObject, eventdata, handles)
rmc=str2num(get(handles.rmcol,'String'));
if rmc<0 | rmc>handles.SizeMF_org(2)
    message=['Error!!! Out of range'];
    msgbox(message);
    set(handles.rmcol,'String','');
    return;
end
handles.Image2remove=rmc;
set(handles.mfdoit,'Enable','on');
guidata(hObject,handles);
% -----------------------
% --- BUTTON DO IT MF ---
% -----------------------
function mfdoit_Callback(hObject, eventdata, handles)
if isempty(get(handles.rmcol,'String'))
    message='Please, enter an image(column) number to delete';
    msgbox(message,'Error','error');
    return;
end
message=['Do you really want to delete image(column) number ',...
          get(handles.rmcol,'String') ' from the marker file?'];
Question=questdlg(message,'Delete image(column)','Yes','No','No');
if strcmp(Question,'No')
    set(handles.rmcol,'String','');
    return;
end
a=tom_emreadc([handles.PathMF handles.NameMF]);
a.Value(:,handles.Image2remove,:)=[];
a.Header.Size(2)=size(a.Value,2);
handles.SizeMF_mod=a.Header.Size;
org_path=pwd;
cd(handles.PathMF);
if get(handles.check_renamemf,'Value')
    tom_emwrite([handles.PathMF handles.NewNameMF],a);
else
    tom_emwrite([handles.PathMF handles.NameMF],a);
end
dimx=num2str(a.Header.Size(1));
dimy=num2str(a.Header.Size(2));
dimz=num2str(a.Header.Size(3));
set(handles.mfsize,'String',[ dimx ' x ' dimy ' x ' dimz ]);
set(handles.mf2,'String',dimy);
set(handles.mf3,'String',dimz);
%set(handles.rmcol,'String','');
message=['Marker file has been modified. Image (column) ' num2str(handles.Image2remove),...
        ' has been successfully removed. Do you want to modified the tilt series as well?'];
Question=questdlg(message,'Operation succeeded','Yes','No','No');
if strcmp(Question,'Yes')
    set(handles.browse_TS,'Enable','on');
    set(handles.tsdir,'Enable','on');
end    
guidata(hObject,handles);
% ---------------------------------
% --- BUTTON BROWSE TILT SERIES ---
% ---------------------------------
function browse_TS_Callback(hObject, eventdata, handles)
uiwait(tom_load_chmf);
m=findobj('Tag','chmf');Param=get(m,'Userdata');
switch Param.newproj_cancel
    case 'yes'
        return; %nothing because Cancel is ckicked 
end
handles.Path_TS=Param.PathName;
handles.Name_TS=Param.FileName;
handles.Ext_TS=Param.Ext;
handles.FirstNb_TS=Param.Firstnb;
handles.LastNb_TS=Param.Lastnb;
handles.RefImg_TS=Param.RefImage;
set(handles.tsdir,'Enable','on','String',Param.PathName);
set(handles.tsfirst,'String',[Param.FileName Param.Firstnb Param.Ext]);
set(handles.tslast,'String',[Param.FileName Param.Lastnb Param.Ext]);
set(handles.tsrefimg,'String',[Param.FileName Param.RefImage Param.Ext]);
file='';
for i=str2num(handles.FirstNb_TS):str2num(handles.LastNb_TS)    
    na=[Param.FileName num2str(i) Param.Ext];
    file=[file;cellstr(na)];
end
handles.List=file;
handles.NewNameTS=[handles.Name_TS 'mod_'];
handles.TSimg2remove=[handles.Name_TS '1' handles.Ext_TS];
set(handles.new_tsname,'Enable','on','String',handles.NewNameTS);
set(handles.tslist,'Enable','on','String',handles.List,'Value',handles.Image2remove);
set(handles.check_renamets,'Enable','on')
set(handles.tsdoit,'Enable','on');
guidata(hObject,handles);
% -----------------------------------
% --- CHECKBOX RENAME TILT SERIES ---
% -----------------------------------
function check_renamets_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    set(handles.new_tsname,'Enable','on');
    handles.NewNameTS=[handles.Name_TS 'mod_'];
else
    message=['It is strongly recommanded to keep the original tiltseries and rename ',...
             'the result with the deleted projection. Do you still want to proceed?'];
    Question=questdlg(message,'Rename tiltseries','Yes','No','No');
    if strcmp(Question,'Yes')     
        set(handles.new_tsname,'Enable','off');
        handles.NewNameTS='';
    else
       set(handles.check_renamets,'Value',1);
    end
end
guidata(hObject,handles);
% -----------------------------------
% --- EDIT BOX RENAME TILT SERIES ---
% -----------------------------------
function new_tsname_Callback(hObject, eventdata, handles)
a=get(hObject,'String');
aa=a(size(a,2));
if ~strcmp(aa,'_')
    a=[a '_'];
end
handles.NewNameTS=a;
guidata(hObject,handles);
% ------------------------------
% --- POPUP MENU TILT SERIES ---
% ------------------------------
function tslist_Callback(hObject, eventdata, handles)
a=get(handles.tslist,'String'); b=get(handles.tslist,'Value');
handles.TSimg2remove=a{b};
guidata(hObject,handles);
% -----------------------
% --- BUTTON DO IT TS ---
% -----------------------
function tsdoit_Callback(hObject, eventdata, handles)
message=['Are you sure you want to delete image ' handles.TSimg2remove ' from the tilt series? '];
Question=questdlg(message,'Delete image from tilt series','Yes','No','No');
if strcmp(Question,'No')
    return;
end
if get(handles.check_renamets,'Value')
end
path_org=pwd;
cd(handles.Path_TS);
j=str2num(handles.FirstNb_TS);
for i=str2num(handles.FirstNb_TS):str2num(handles.LastNb_TS)
    na=[handles.Name_TS num2str(i) handles.Ext_TS];
    if ~strcmp(na,handles.TSimg2remove)
        if get(handles.check_renamets,'Value')==1            
            copyfile(na,[handles.NewNameTS num2str(j) handles.Ext_TS]);
        else
            if exist([handles.Name_TS num2str(j) handles.Ext_TS])==0
                movefile(na,[handles.Name_TS num2str(j) handles.Ext_TS]);
            end
        end
        j=j+1;
    else
        if get(handles.check_renamets,'Value')==0
            delete(handles.TSimg2remove);
        end
    end
end
cd(path_org);
message=['The tilt series has been modified. Image ' num2str(handles.TSimg2remove),...
        ' has been successfully removed.'];
msgbox(message,'Success');
guidata(hObject,handles);
% ----------------------------
% --- BUTTON NEW SELECTION ---
% ----------------------------
function button_news_Callback(hObject, eventdata, handles)
init;
% -------------------
% --- BUTTON EXIT ---
% -------------------
function exit_menu_Callback(hObject, eventdata, handles)
delete(gcf);

% ------------------------
% ---  OTHER FUNCTION  ---
% ------------------------
function init
h_chmf=findobj(0,'Tag','chmf');handles=guidata(h_chmf);
set(handles.mfdir,'String','To select a Marker File --------- Click button Browse --------->');
set(handles.mfname,'String','');
set(handles.check_renamemf,'Value',0,'Enable','off');
set(handles.new_mfname,'Enable','off','String','');
set(handles.mfsize,'String','');
set(handles.mf2,'String','');
set(handles.mf3,'String','');
set(handles.rmcol,'Enable','off','String','');
set(handles.mfdoit,'Enable','off');
set(handles.tsdir,'Enable','off','String','To select a tilt series --------- Click button Browse --------->');
set(handles.browse_TS,'Enable','off');
set(handles.tsfirst,'String','');
set(handles.tslast,'String','');
set(handles.tsrefimg,'String','');
set(handles.check_renamets,'Value',1,'Enable','off');
set(handles.new_tsname,'String','');
set(handles.tslist,'Value',1,'String','-- none --','Enable','off')
set(handles.tsdoit,'Enable','off');

handles.PathMF='';
handles.NameMF='';
handles.NewNameMF='';
handles.SizeMF_org='';
handles.SizeMF_mod='';
handles.LastProj='';
handles.Image2remove='';
handles.Path_TS='';
handles.Name_TS='';
handles.Ext_TS='';
handles.FirstNb_TS='';
handles.LastNb_TS='';
handles.RefImg_TS='';
handles.TSimg2remove='';
handles.NewNameTS='';
handles.List='';
guidata(h_chmf,handles);

