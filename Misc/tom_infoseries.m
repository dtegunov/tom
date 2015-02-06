function varargout = tom_infoseries(varargin)
%TOM_INFOSERIES generates a tilt series documentation output
%
%   varargout = tom_infoseries(varargin)
%
%This documentation can be printed
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
%   tom_infoseries      opens a dialog box
%
%REFERENCES
%
%SEE ALSO
%   TOM_EMREAD, TOM_SETMARK
%
%   created by SN 10/14/02
%   updated by WDN 12/01/06
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
                   'gui_OpeningFcn', @tom_infoseries_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_infoseries_OutputFcn, ...
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


% --- Executes just before tom_infoseries is made visible.
function tom_infoseries_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for tom_infoseries
handles.output = hObject;

[filename, pathname] = uigetfile('*.alg', 'Load alignment file');
if isequal(filename,0) | isequal(pathname,0)
    disp('No data loaded.');
else
    mm=fopen([pathname filename],'r');
    textp=fscanf(mm,'%c',79);
    textp=fscanf(mm,'%s\n',1);param.mypathname=textp;
    textp=fscanf(mm,'%s\n',1);param.myfilename=textp;
    textp=fscanf(mm,'%s\n',1);param.myfirstnb=textp;
    textp=fscanf(mm,'%s\n',1);param.mylastnb=textp;
    textp=fscanf(mm,'%s\n',1);param.myext=textp;
    textp=fscanf(mm,'%s\n',1);param.myfilemarker_default=textp;
    textp=fscanf(mm,'%s\n',1);param.myfilemarker=textp;
    textp=fscanf(mm,'%s\n',1);param.image_ref=textp;
    textp=fscanf(mm,'%s\n',1);param.newproj_cancel=textp;
    fclose(mm);
    em_name=[param.mypathname param.myfilename param.image_ref param.myext];
    i=tom_emread(em_name);
    [meanv max min std]=tom_dev(i.Value);
        if (meanv-2*std)>=(meanv+3*std)
        imagesc(i.Value');
    else
        imagesc(i.Value',[meanv-(3*std) meanv+(3*std)]);
    end;
    a=regexprep(em_name,'_','\\_');
    title(a);
    colormap gray;
    set(handles.PathProjIma,'String',i.Header.Pathname);
    set(handles.RefProjIma,'String',i.Header.Filename);
    set(handles.SizeX,'String',i.Header.Size(1));
    set(handles.SizeY,'String',i.Header.Size(2));
    set(handles.SizeZ,'String',i.Header.Size(3));
    set(handles.Mag,'String',i.Header.Magnification);
    set(handles.ExpoTime,'String',i.Header.Exposuretime);
    set(handles.ObjPix,'String',i.Header.Objectpixelsize);
    set(handles.TiltAngle,'String',i.Header.Tiltangle);
    %reas marker file
    m=tom_emreadc(param.myfilemarker);
    handles.Matrixmark=m.Value;
    a=regexprep(num2str(handles.Matrixmark(1,:,1)),'      ',' / ');
    set(handles.NbAngles,'String', a);
    r=[0 0 0 ];
    [handles.Matrixmark, psi, sigma, x, y, z]  = tom_alignment3d(handles.Matrixmark, 1, param.image_ref, r, i.Header.Size(1));
    set(handles.NbTilts,'String',num2str(size(handles.Matrixmark,2)));
    set(handles.NbProj,'String',num2str(param.image_ref));
    set(handles.NbMk,'String',num2str(size(handles.Matrixmark,3)));
    set(handles.NbRef,'String',num2str(1));
    set(handles.TiltAzimuth,'String',num2str(psi.*180./pi));
    set(handles.RMS,'String',num2str(sigma));
    for lauf=1:size(handles.Matrixmark,3)%
        x=double(handles.Matrixmark(2,str2num(param.image_ref),lauf))-44;
        y=double(handles.Matrixmark(3,str2num(param.image_ref),lauf))+6;
        text(x,y,['O' num2str(lauf)],'FontName','FixedWidth','Fontsize',22,'Fontweight','bold')%text(x,y,'O','FontName','FixedWidth','Fontsize',22,'Fontweight','bold');
    end;
    handles.Image=i;
    handles.Param=param;
end;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = tom_infoseries_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;
%----------------------
%--- BUTTON DETAILS ---
%----------------------
function Details_Callback(hObject, eventdata, handles)
h=findobj(0,'Tag','infoseries');set(h,'Units','pixel');
pos_info=get(h,'Position');
set(h,'Units','normalized');
h=tom_editemheader(handles.Image);
set(h,'Units','pixel');
pos_edit=get(h,'Position');
set(h,'Position',[pos_info(1)+pos_info(3)+6 pos_info(2) pos_edit(3) pos_edit(4)]);
set(h,'Units','normalized');
%Set enable off all text of tom_editenheader
hh=guidata(h);
set(hh.load,'Enable','off');
set(hh.filename,'Enable','off');
set(hh.pathname,'Enable','off');
set(hh.sizex,'Enable','off');
set(hh.sizey,'Enable','off');
set(hh.sizez,'Enable','off');
set(hh.size_type,'Enable','off');
set(hh.voltage,'Enable','off');
set(hh.cs,'Enable','off');
set(hh.aperture,'Enable','off');
set(hh.mag,'Enable','off');
set(hh.postmag,'Enable','off');
set(hh.Expo_time,'Enable','off');
set(hh.Obj_pix,'Enable','off');
set(hh.microscope,'Enable','off');
set(hh.pixsize,'Enable','off');
set(hh.ccdarea,'Enable','off');
set(hh.defocus,'Enable','off');
set(hh.astigmatism,'Enable','off');
set(hh.astigmatism_angle,'Enable','off');
set(hh.focusincrement,'Enable','off');
set(hh.CountsPerElectron,'Enable','off');
set(hh.intensity,'Enable','off');
set(hh.EnergySlitwidth,'Enable','off');
set(hh.EnergyOffset,'Enable','off');
set(hh.tiltangle,'Enable','off');
set(hh.tiltaxis,'Enable','off');
set(hh.Marker_X,'Enable','off');
set(hh.Marker_Y,'Enable','off');
set(hh.comment,'Enable','off');
set(hh.apply,'Enable','off');
set(hh.Reset,'Enable','off');
%----------------------
%--- BUTTON PREVIEW ---
%----------------------
function Preview_Callback(hObject, eventdata, handles)
h=findobj(0,'Tag','infoseries');printpreview(h);
h=findobj(0,'Tag','EEH');
if ~isempty(h)
    printpreview(h);
end
%--------------------
%--- Button PRINT ---
%--------------------
function Print_Callback(hObject, eventdata, handles)
h=findobj(0,'Tag','infoseries');
print(h);
h=findobj(0,'Tag','EEH');
if ~isempty(h)
    print(h);
end



