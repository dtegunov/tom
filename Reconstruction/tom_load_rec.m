function varargout = tom_load_rec_tiltseries(varargin)
%TOM_LOAD_REC_TILTSERIES is called by tom_rec3d only
%
%   varargout = tom_load_rec_tiltseries(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   arargout    		...
%
%EXAMPLE
%   ... = tom_load_rec_tiltseries(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
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

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_load_rec_tiltseries_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_load_rec_tiltseries_OutputFcn, ...
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
% --- Executes just before untitled is made visible.
function tom_load_rec_tiltseries_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = tom_load_rec_tiltseries_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% -----------------------------------------------------------
% -------   BUTTON BROWSE   ---------------------------------
% -----------------------------------------------------------
function lb_browse_Callback(hObject, eventdata, handles)
set(handles.lb_path,'String','');
set(handles.lb_file,'String','');
set(handles.lb_firstnb,'String','');
set(handles.lb_lastnb,'String','');
set(handles.lb_ext,'String','');

[f, p] = uigetfile('*.*', '---- Click on the last number of your tilt series ----');
if isequal(f,0) | isequal(p,0) 
    %error('No data loaded.'); 
    return;%nothing because Cancel is ckicked
else
    myext='';myfile='';mynb='';        
    p=strrep(p,'\','/');
    set(handles.lb_path,'String',p);
    fext=findstr(f,'.');
    if ~isempty(fext)
        fext=fext(size(fext,2));
        myext=f(fext+1:size(f,2));
    else
        fext=size(f,2)+1;
    end
    set(handles.lb_ext,'String',myext);
    fnb=findstr(f,'_');
    if ~isempty(fnb)
        a=size(fnb,2);
        mynb=f(fnb(a)+1:fext-1);
    end
    set(handles.lb_lastnb,'String',mynb);
    myfile=f(1:fnb(a));
    set(handles.lb_file,'String',myfile);
    set(handles.lb_firstnb,'String','1');
end
% -----------------------------------------------------------
% -------   BUTTON OK   -------------------------------------
% -----------------------------------------------------------
function lb_ok_Callback(hObject, eventdata, handles)
Param.PathName=get(handles.lb_path,'String');
Param.FileName=get(handles.lb_file,'String');
Param.Firstnb=get(handles.lb_firstnb,'String');
Param.Lastnb=get(handles.lb_lastnb,'String');
Param.Ext=['.' get(handles.lb_ext,'String')];
aa=[Param.PathName Param.FileName '1' Param.Ext];
ah=tom_reademheader(aa);
tilt(1)=ah.Header.Parameter(19)/1000;
aa=[Param.PathName Param.FileName '2' Param.Ext];
ah=tom_reademheader(aa);
tilt(2)=ah.Header.Parameter(19)/1000;
if abs(tilt(1)) < abs(tilt(2))
    little=abs(tilt(1));
    No_ima=1;
else
    little=abs(tilt(2));
    No_ima=2;
end
for i=3:str2num(Param.Lastnb)
    aa=[Param.PathName Param.FileName num2str(i) Param.Ext];
    ah=tom_reademheader(aa);
    tilt(i)=ah.Header.Parameter(19)/1000;
    if abs(tilt(i)) < little
        little=abs(tilt(i));
        No_ima=i;
    end
end
Param.Angle=tilt;
Param.RefImage=num2str(No_ima);
Param.newproj_cancel='no';
set(findobj(0,'Tag','rec3d'),'Userdata',Param);
close(handles.tom_load_rec); % tom_load_rec will be closed

% -----------------------------------------------------------
% -------   BUTTON CANCEL   ---------------------------------
% -----------------------------------------------------------
function lb_cancel_Callback(hObject, eventdata, handles)
Param.newproj_cancel='yes';
set(findobj(0,'Tag','rec3d'),'Userdata',Param);
close(handles.tom_load_rec); % tom_load_rec will be closed

