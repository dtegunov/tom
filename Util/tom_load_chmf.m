function varargout = tom_load_chmf(varargin)
%TOM_LOAD_CHMF creates ...
%
%   varargout = tom_load_chmf(varargin)
%
%   tom_load_chmf_tiltseries is called by tom_rec3d only
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
%   ... = tom_load_chmf(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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
                   'gui_OpeningFcn', @tom_load_chmf_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_load_chmf_OutputFcn, ...
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
function tom_load_chmf_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = tom_load_chmf_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% -----------------------------------------------------------
% -------   BUTTON BROWSE   ---------------------------------
% -----------------------------------------------------------
function lb_browse_Callback(hObject, eventdata, handles)
[f, p] = uigetfile('*.*', '---- Click on the last number of your tilt series ----');
sp=size(p);
sp=sp(2);
mypath=[];
for i=1:sp
    a=p(i);    
    switch a
    case '\'        
        b='/';
        mypath=strcat(mypath,b);
    otherwise            
        mypath=[mypath a];
    end
end
hh=findobj('Tag','lb_path');set(hh,'String',mypath);
s=size(f);
s=s(2);
for i=1:s
    a=f(i);
    switch a
    case{'_'}
        c=[];
        for j=1:i
            c=strcat(c,f(j));
        end
        hh=findobj('Tag','lb_file');
        set(hh,'String',c);
        c=[];
        for j=1:3
            b=f(i+j);            
            switch b
            case{'0','1','2','3','4','5','6','7','8','9'}      
                c=strcat(c,b);
            case{'.'}
                hh=findobj('Tag','lb_firstnb');
                set(hh,'String','1');
                hh=findobj('Tag','lb_lastnb');
                set(hh,'String',c);                
                break;
            end
        end
    case{'.'}
        c=[];
        for j=i+1:s
            c=strcat(c,f(j));
        end
        hh=findobj('Tag','lb_ext');
        set(hh,'String',c);
    end
end
hh=findobj('Tag','lb_filemarker');
set(hh,'String','markfile');

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
set(findobj(0,'Tag','chmf'),'Userdata',Param);
close(handles.tom_load_chmf); % tom_load_chmf will be closed

% -----------------------------------------------------------
% -------   BUTTON CANCEL   ---------------------------------
% -----------------------------------------------------------
function lb_cancel_Callback(hObject, eventdata, handles)
Param.newproj_cancel='yes';
set(findobj(0,'Tag','chmf'),'Userdata',Param);
close(handles.tom_load_chmf); % tom_load_chmf will be closed

