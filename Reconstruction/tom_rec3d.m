function varargout = tom_rec3d(varargin)
%TOM_REC3D makes a reconstruction 3D
%
%   varargout = tom_rec3d(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout   		...
%
%EXAMPLE
%   ... = tom_rec3d(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_TILT_LINES, TOM_EMREAD, TOM_EMWRITE, TOM_REC3D
%
%   created by WDN 09/12/02
%   updated by WDN 06/14/05
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
                   'gui_OpeningFcn', @tom_rec3d_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_rec3d_OutputFcn, ...
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
function tom_rec3d_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
InitRec3d(handles);
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = tom_rec3d_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
%*********************************************************************
% --------------------------------------------------------------------
% -------   MENU NEW   -----------------------------------------------
% --------------------------------------------------------------------
function new_Callback(hObject, eventdata, handles)
InitRec3d(handles);
handles.Param=[];
set(handles.Do_reconstruction,'String','DO reconstruction');
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   MENU HELP   -----------------------------------------------
% --------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)
set(handles.label_pre,'ForegroundColor',[0 0 1]);
set(handles.label_post,'ForegroundColor',[0 0 1]);
set(handles.Binning,'ForegroundColor',[0 0 1]);
set(handles.Filter,'ForegroundColor',[0 0 1]);
set(handles.Weighting,'ForegroundColor',[0 0 1]);
set(handles.wexact,'ForegroundColor',[0 0 1]);
set(handles.wrweight,'ForegroundColor',[0 0 1]);
set(handles.label_size,'ForegroundColor',[0 0 1]);
set(handles.label_abs,'ForegroundColor',[0 0 1]);
uiwait(msgbox('Click with the right click of the mouse on the text in blue to have more information','Help','help'));
message_binning1=['Binning is used to reduced the number of pixel. For example binning 1'];
message_binning2=['means that an area of 2 adjacent pixels have been combined into one'];
message_binning3=['larger pixel, binning 2 with 4 adjacent pixels into one larger pixel and so on ...'];
message_prebinning1='Pre binning means that the image is binned first, ';
message_prebinning2='and then aligned, rotated and weighted.';
message_postbinning1='Post binning means that the image is aligned, ';
message_postbinning2='rotated, and weighted first and then binned.';
message_weighting='Adapt the contrast of the image in Fourrier Space before the backprojection.';
message_filter='The filter is used to remove the noise before the weighting.';
message_exactweighting1='Should be the same number as size Z of the volume.';
message_exactweighting2='Same weighting as in single particle analysis.';
message_rweighting1='Perform an analytical weighting. A r-weight is proportional. ';
message_rweighting2='to the frquency perpendiculary to the tilt axis.';
message_size='Size (x,y,z) of the reconstruction';
message_abscenter1='Coordinate (x,y) of Reference marker. The parameter Z is equal';
message_abscenter2='to the size X of the image + 1 pixel.';
%label binning
cmenu = uicontextmenu;
set(handles.Binning,'UIContextMenu',cmenu);
h=uimenu(cmenu,'Label',message_binning1,'Tag','help_bin');
h=uimenu(cmenu,'Label',message_binning2,'Tag','help_bin');
h=uimenu(cmenu,'Label',message_binning3,'Tag','help_bin');
%label pre
cmenu = uicontextmenu;
set(handles.label_pre,'UIContextMenu',cmenu);
h=uimenu(cmenu,'Label',message_prebinning1,'Tag','help_prebinning');
h=uimenu(cmenu,'Label',message_prebinning2,'Tag','help_prebinning');
%label post
cmenu = uicontextmenu;
set(handles.label_post,'UIContextMenu',cmenu);
h=uimenu(cmenu,'Label',message_postbinning1,'Tag','help_postbinning');
h=uimenu(cmenu,'Label',message_postbinning2,'Tag','help_postbinning');
%label filter
cmenu = uicontextmenu;
set(handles.Filter,'UIContextMenu',cmenu);
h=uimenu(cmenu,'Label',message_filter,'Tag','help_filter');
%label weighting
cmenu = uicontextmenu;
set(handles.Weighting,'UIContextMenu',cmenu);
h=uimenu(cmenu,'Label',message_weighting,'Tag','help_wei');
%label r-weight weighting
cmenu = uicontextmenu;
set(handles.wexact,'UIContextMenu',cmenu);
h=uimenu(cmenu,'Label',message_exactweighting1,'Tag','help_rwei');
h=uimenu(cmenu,'Label',message_exactweighting2,'Tag','help_rwei');
%label exact weighting
cmenu = uicontextmenu;
set(handles.wrweight,'UIContextMenu',cmenu);
h=uimenu(cmenu,'Label',message_rweighting1,'Tag','help_exactwei');
h=uimenu(cmenu,'Label',message_rweighting2,'Tag','help_exactwei');
%label Size (pix)
cmenu = uicontextmenu;
set(handles.label_size,'UIContextMenu',cmenu);
h=uimenu(cmenu,'Label',message_size,'Tag','help_size');
%label Absolute center
cmenu = uicontextmenu;
set(handles.label_abs,'UIContextMenu',cmenu);
h=uimenu(cmenu,'Label',message_abscenter1,'Tag','help_abs');
h=uimenu(cmenu,'Label',message_abscenter2,'Tag','help_abs');
% --------------------------------------------------------------------
% -------   BUTTON BROWSE ALIGNMENT FILE   ---------------------------
% --------------------------------------------------------------------
function browseAF_Callback(hObject, eventdata, handles)
[f, p] = uigetfile('*.alg','Load alignment file');  
if isequal(f,0)|isequal(p,0)
    return;%nothing because Cancel is ckicked    
else
    conf_file=[p f];
    %mm=fopen(conf_file,'r');
    %text=fscanf(mm,'%c',79);
    s=textread(conf_file,'%s','whitespace','\n');
    text=char(s(1));
    switch text
        case 'Alignment file use by tom_setmarker when the button ''Load'' a project is clicked'
            %text=fscanf(mm,'%s\n',1);
            %if (strcmp(text,'rel'))
                %Param.PathName=''; 
            %else
                %Param.PathName=text;
            %end;
            Param.PathName=char(s(2));
            %text=fscanf(mm,'%s\n',1);Param.FileName=text;
            Param.FileName=char(s(3));
            %text=fscanf(mm,'%s\n',1);Param.Firstnb=text;
            Param.Firstnb=char(s(4));
            %text=fscanf(mm,'%s\n',1);Param.Lastnb=text;
            Param.Lastnb=char(s(5));
            %text=fscanf(mm,'%s\n',1);Param.Ext=text;
            Param.Ext=char(s(6));
            %text=fscanf(mm,'%s\n',1);Param.MarkerFile_default=text;
            Param.MarkerFile_default=char(s(7));
            %text=fscanf(mm,'%s\n',1);Param.MarkerFiler=text;
            Param.MarkerFiler=char(s(8));
            %text=fscanf(mm,'%s\n',1);Param.RefImage=text;
            Param.RefImage=char(s(9));
            %fclose(mm);
            set(handles.AlignmentFileName,'String',[p f]); 
            set(handles.MarkerFileName,'String',Param.MarkerFiler);
            set(handles.TSPath,'String',Param.PathName);                
            set(handles.TSFirst,'String',[Param.FileName Param.Firstnb Param.Ext]);%      
            set(handles.TSLast,'String',[Param.FileName Param.Lastnb Param.Ext]);
            set(handles.TSRefim,'String',Param.RefImage);
            %Param.RefImage=str2num(Param.RefImage);                      
            reh=tom_reademheader([Param.PathName Param.FileName Param.Firstnb Param.Ext]);
            set(handles.sizex,'String',reh.Header.Size(1)./4);
            set(handles.sizey,'String',reh.Header.Size(2)./4);
            set(handles.sizez,'String',reh.Header.Size(1)./8);
            Param.Projection_Size_X=reh.Header.Size(1);
            Param.Projection_Size_Y=reh.Header.Size(2);
            set(handles.offsetx,'String','0');
            set(handles.offsety,'String','0');
            set(handles.offsetz,'String','0');
            temp=tom_emread(Param.MarkerFiler);
            Param.Matrixmark=temp.Value;
            Param.AbsoluteCenter=[(reh.Header.Size(1)./2 +1) (reh.Header.Size(2)./2 +1) (reh.Header.Size(1)./2 +1)];%changed FF
            set(handles.absx,'String',[Param.AbsoluteCenter(1)]);
            set(handles.absy,'String',[Param.AbsoluteCenter(2)]);        
            set(handles.absz,'String',[Param.AbsoluteCenter(3)]);
            set(handles.ExactWeighting,'String',reh.Header.Size(1)./8);       
            Param.Header=temp.Header;
            Nbm=size(Param.Matrixmark,3);
            j=1;
            for i=1:Nbm
                if Param.Matrixmark(2,str2num(Param.RefImage),i)~=-1
                    mat(j)=i;
                    j=j+1;
                end
            end
            if isempty(mat)
                mat=1
            end
            set(handles.Refmarker,'String',mat,'Value',1);
            %[tempmat,beta,error,x,y,z]=tom_alignment3d(Param.Matrixmark,1);
            [tempmat,beta,error,x,y,z]=tom_alignment3d(Param.Matrixmark,1,str2num('handles.Refmarker'),Param.AbsoluteCenter,Param.Projection_Size_X);
            set(handles.error,'String',num2str(error))
            %extract information for demomode
            psi_indeg=(beta.*(180./pi))+90;
            psi=psi_indeg.*pi./180;
            shift_tx = tempmat(7,str2num(Param.RefImage),1);
            shift_ty = tempmat(8,str2num(Param.RefImage),1);
            shiftx = cos(psi)*shift_tx + sin(psi)*shift_ty;
            shifty = -sin(psi)*shift_tx + cos(psi)*shift_ty;
            Param.Shift_RefImage(1)=shiftx;
            Param.Shift_RefImage(2)=shifty;
            Param.Rotation_RefImage=psi_indeg;
        otherwise
            uiwait(errordlg('The file selected is not a file for alignment. Please check your file. Your file should start with: "Alignment file use by tom_setmarker when the button ''Load'' a project is clicked"','File Error'));
            fclose(mm);
    end
end
Param.ReconstructionName=[Param.PathName 'reconstruction.vol'];
set(handles.VolumeName,'String',Param.ReconstructionName);
Param.WeightedFileName=[Param.PathName 'TEMP_BPP_mat_'];
set(handles.WeightedFileName,'String',Param.WeightedFileName);
handles.Param=Param;
ref_mark=1;
handles.Param.AbsoluteCenter=[handles.Param.Matrixmark(2,str2num(handles.Param.RefImage),ref_mark),...
    handles.Param.Matrixmark(3,str2num(handles.Param.RefImage),ref_mark),...
    (handles.Param.Projection_Size_X./2 +1)];
set(handles.absx,'String',[handles.Param.AbsoluteCenter(1)]);
set(handles.absy,'String',[handles.Param.AbsoluteCenter(2)]);            
set(handles.absz,'String',[handles.Param.AbsoluteCenter(3)]);
set(handles.LowpassFilter,'String','0.8');
EnablePanelRec3d('on');
EnablePanelMarkerFile('on');
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   BUTTON BROWSE TILTSERIES  ------------------------------------
% --------------------------------------------------------------------
function browseTS_Callback(hObject, eventdata, handles)
uiwait(tom_load_rec);
m=findobj('Tag','rec3d');
handles.Param=get(m,'Userdata');
switch handles.Param.newproj_cancel
    case 'yes'
        return; %nothing because Cancel is ckicked 
end
set(handles.TSPath,'String',handles.Param.PathName);
set(handles.TSFirst,'String',[handles.Param.FileName handles.Param.Firstnb handles.Param.Ext]); 
set(handles.TSLast,'String',[handles.Param.FileName handles.Param.Lastnb handles.Param.Ext]); 
set(handles.TSRefim,'String',handles.Param.RefImage); 
reh=tom_reademheader([handles.Param.PathName handles.Param.FileName handles.Param.Firstnb handles.Param.Ext]);
set(handles.sizex,'String',reh.Header.Size(1)./4);
set(handles.sizey,'String',reh.Header.Size(2)./4);
set(handles.sizez,'String',reh.Header.Size(1)./8);
handles.Param.Projection_Size_X=reh.Header.Size(1);
handles.Param.Projection_Size_Y=reh.Header.Size(2);
set(handles.offsetx,'String','0');
set(handles.offsety,'String','0');
set(handles.offsetz,'String','0');
set(handles.MarkerFileName,'String','');
set(handles.absx,'String','');
set(handles.absy,'String','');
set(handles.absz,'String','');
handles.Param.AbsoluteCenter=[(reh.Header.Size(1)./2 +1) (reh.Header.Size(2)./2 +1) (reh.Header.Size(1)./2 +1)];%changed FF
set(handles.ExactWeighting,'String',reh.Header.Size(1)./8);
set(handles.MarkerFileName,'String',' ');
EnablePanelMarkerFile('on');
guidata(hObject,handles);

% --------------------------------------------------------------------
% -------   BUTTON BROWSE MARKER FILE   ------------------------------
% --------------------------------------------------------------------
function browseMF_Callback(hObject, eventdata, handles)
actdir=pwd;
cd (handles.Param.PathName);
[f, p] = uigetfile('*.em','Load marker file');  
if isequal(f,0)|isequal(p,0)
    return;%nothing because Cancel is ckicked
else
    m=findobj('Tag','r3d_marker_file');set(m,'String',[p f]);                        
    temp=tom_emread([p f]);    
    set(handles.absx,'String',[handles.Param.AbsoluteCenter(1)]);
    set(handles.absy,'String',[handles.Param.AbsoluteCenter(2)]);      
    set(handles.absz,'String',[handles.Param.AbsoluteCenter(3)]);
    handles.Param.Matrixmark=temp.Value;
    handles.Param.Header=temp.Header;
    handles.Param.MarkerFiler=[p f];
    set(handles.MarkerFileName,'String',handles.Param.MarkerFiler)
    Nbm=size(handles.Param.Matrixmark,3);
    j=1;
    for i=1:Nbm
        if handles.Param.Matrixmark(2,str2num(handles.Param.RefImage),i)~=-1
            mat(j)=i;
            j=j+1;
        end
    end
    if isempty(mat)
        mat=1
    end
    set(handles.Refmarker,'String',mat,'Value',1);
    set(handles.WeightedFileName,'String',[handles.Param.PathName 'TEMP_BPP_mat_']);
    set(handles.VolumeName,'String',[handles.Param.PathName 'reconstruction.vol']);
    set(handles.LowpassFilter,'String','0.8');
    EnablePanelRec3d('on');
    [tempmat,beta,error,x,y,z]=tom_alignment3d(handles.Param.Matrixmark,1,str2num('handles.Refmarker'),handles.Param.AbsoluteCenter,handles.Param.Projection_Size_X);
    %extract information for demomode
    psi_indeg=(beta.*(180./pi))+90;
    psi=psi_indeg.*pi./180;
    shift_tx = tempmat(7,str2num(handles.Param.RefImage),1);
    shift_ty = tempmat(8,str2num(handles.Param.RefImage),1);
    shiftx = cos(psi)*shift_tx + sin(psi)*shift_ty;
    shifty = -sin(psi)*shift_tx + cos(psi)*shift_ty;
    handles.Param.Shift_RefImage(1)=shiftx;
    handles.Param.Shift_RefImage(2)=shifty;
    handles.Param.Rotation_RefImage=psi_indeg;
    set(handles.error,'String',num2str(error))
    ref_mark=1; 
    handles.Param.AbsoluteCenter=[handles.Param.Matrixmark(2,str2num(handles.Param.RefImage),ref_mark),...
    handles.Param.Matrixmark(3,str2num(handles.Param.RefImage),ref_mark),...
    (handles.Param.Projection_Size_X./2 +1)];
    set(handles.absx,'String',[handles.Param.AbsoluteCenter(1)]);
    set(handles.absy,'String',[handles.Param.AbsoluteCenter(2)]);            
    set(handles.absz,'String',[handles.Param.AbsoluteCenter(3)]);

end
guidata(hObject,handles);
cd (actdir);
% --------------------------------------------------------------------
% -------   POPUP MENU REFERENCE MARKER   ----------------------------
% --------------------------------------------------------------------
function Refmarker_Callback(hObject, eventdata, handles)
a=get(hObject,'String');
b=get(hObject,'Value');
ref_mark=str2num(a(b));
if ischar(handles.Param.RefImage)    
    handles.Param.AbsoluteCenter=[handles.Param.Matrixmark(2,str2num(handles.Param.RefImage),ref_mark),...
        handles.Param.Matrixmark(3,str2num(handles.Param.RefImage),ref_mark) ,...
        (handles.Param.Projection_Size_X./2 +1)];
else
    handles.Param.AbsoluteCenter=[handles.Param.Matrixmark(2,handles.Param.RefImage,ref_mark) handles.Param.Matrixmark(3, ...
            handles.Param.RefImage,ref_mark) (handles.Param.Projection_Size_X./2 +1)];
end    
set(handles.absx,'String',[handles.Param.AbsoluteCenter(1)]);
set(handles.absy,'String',[handles.Param.AbsoluteCenter(2)]);            
set(handles.absz,'String',[handles.Param.AbsoluteCenter(3)]);
[tempmat,psi,error,x,y,z]=tom_alignment3d(handles.Param.Matrixmark,ref_mark);
set(handles.error,'String',num2str(error))
guidata(hObject,handles);

% --------------------------------------------------------------------
% -------   EDIT SIZE X   --------------------------------------------
% --------------------------------------------------------------------
function sizex_Callback(hObject, eventdata, handles)
if get(handles.demo,'Value')
    h_axe=findall(findobj(0,'Tag','demo1'),'Type','axe');
    if size(h_axe,1)>1
        handles=InitDemoImage(handles);
    end
    DemoRedrec(handles);
end
% --------------------------------------------------------------------
% -------   EDIT SIZE Y   --------------------------------------------
% --------------------------------------------------------------------
function sizey_Callback(hObject, eventdata, handles)
if get(handles.demo,'Value')
    h_axe=findall(findobj(0,'Tag','demo1'),'Type','axe');
    if size(h_axe,1)>1
        handles=InitDemoImage(handles);
    end
    DemoRedrec(handles);
end
% --------------------------------------------------------------------
% -------   EDIT SIZE Z   --------------------------------------------
% --------------------------------------------------------------------
function sizez_Callback(hObject, eventdata, handles)
set(handles.ExactWeighting,'String',get(hObject,'String'));
% --------------------------------------------------------------------
% -------   EDIT OFFSET X   ------------------------------------------
% --------------------------------------------------------------------
function offsetx_Callback(hObject, eventdata, handles)
if get(handles.demo,'Value')
    h_axe=findall(findobj(0,'Tag','demo1'),'Type','axe');
    if size(h_axe,1)>1
        handles=InitDemoImage(handles);
    end
    DemoRedrec(handles);
end
% --------------------------------------------------------------------
% -------   EDIT OFFSET Y   ------------------------------------------
% --------------------------------------------------------------------
function offsety_Callback(hObject, eventdata, handles)
if get(handles.demo,'Value')
    h_axe=findall(findobj(0,'Tag','demo1'),'Type','axe');
    if size(h_axe,1)>1
        handles=InitDemoImage(handles);
    end
    DemoRedrec(handles);
end
% --------------------------------------------------------------------
% -------   BUTTON SAVE AS 3D RECONSTRUCTION FILE   ------------------
% --------------------------------------------------------------------
function browseVOL_Callback(hObject, eventdata, handles)
actdir=pwd;
cd (handles.Param.PathName);
[f, p] = uiputfile('*.vol','SAVE RECONSTRUCTION AS');
if isequal(f,0)|isequal(p,0)
    %nothing because Cancel is ckicked
else
    myfile=[p f];
    if isempty(findstr(myfile,'.vol'))
        myfile=strcat(myfile,'.vol');
    end
    set(handles.VolumeName,'String',myfile);                       
    handles.Param.ReconstructionName=myfile;
end
guidata(hObject,handles);
cd (actdir);
% --------------------------------------------------------------------
% -------   EDIT NAME OF WEIGHTED FILE   -----------------------------
% --------------------------------------------------------------------
function WeightedFileName_Callback(hObject, eventdata, handles)
wfn=get(hObject,'String');
if strcmp(wfn(size(wfn,2)),'_')==0;
    wfn=[wfn '_'];
    set(hObject,'String',wfn);
end
% --------------------------------------------------------------------
% -------   BUTTON SAVE AS WEIGHTED FILE   ---------------------------
% --------------------------------------------------------------------
function browseWEI_Callback(hObject, eventdata, handles)
actdir=pwd;
cd (handles.Param.PathName);
[f, p] = uiputfile('*.*','NAME OF THE WEIGHTED FILES');
if isequal(f,0)|isequal(p,0)
    %nothing because Cancel is ckicked
else
    if strcmp(f(size(f,2)),'_')==0;
        f=[f '_'];        
    end
    set(handles.WeightedFileName,'String',[p f]);
end
cd (actdir);
% --------------------------------------------------------------------
% -------   CHECKBOX WEIGHTED FILE ONLY   ----------------------------
% --------------------------------------------------------------------
function OnlyWeigthedFile_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
% -------   POPUP MENU PRE BINNING   ---------------------------------
% --------------------------------------------------------------------
function prebinning_Callback(hObject, eventdata, handles)
if get(handles.demo,'Value')
    h_axe=findall(findobj(0,'Tag','demo1'),'Type','axe');
    if size(h_axe,1)>1
        handles=InitDemoImage(handles);
    end
    DemoRedrec(handles);
end
% --------------------------------------------------------------------
% -------   POPUP MENU POST BINNING   --------------------------------
% --------------------------------------------------------------------
function postbinning_Callback(hObject, eventdata, handles)
if get(handles.demo,'Value')
    h_axe=findall(findobj(0,'Tag','demo1'),'Type','axe');
    if size(h_axe,1)>1
        handles=InitDemoImage(handles);
    end
    DemoRedrec(handles);
end
% --------------------------------------------------------------------
% -------   EDIT FILTER   --------------------------------------------
% --------------------------------------------------------------------
function LowpassFilter_Callback(hObject, eventdata, handles)
asd=get(handles.LowpassFilter,'String');
s=size(asd);
s=s(2);
filter=[];
for i=1:s
    a=asd(i);
    switch a
        case ','
            b='.';
            filter=strcat(filter,b);
        otherwise
            filter=strcat(filter,a);
    end
end
filter=str2num(filter);
set(handles.LowpassFilter,'String',filter);
% --------------------------------------------------------------------
% -------   RADIO WEIGHTING NONE   -----------------------------------
% --------------------------------------------------------------------
function wnone_Callback(hObject, eventdata, handles)
set(handles.wnone,'Value',1);
set(handles.wrweight,'Value',0);
set(handles.wexact,'Value',0);
set(handles.ExactWeighting,'Enable','off');
% --------------------------------------------------------------------
% -------   RADIO WEIGHTING R_WEIGHT   -------------------------------
% --------------------------------------------------------------------
function wrweight_Callback(hObject, eventdata, handles)
set(handles.wnone,'Value',0);
set(handles.wrweight,'Value',1);
set(handles.wexact,'Value',0);%
set(handles.ExactWeighting,'Enable','off');
% --------------------------------------------------------------------
% -------   RADIO WEIGHTING EXACT   ----------------------------------
% --------------------------------------------------------------------
function wexact_Callback(hObject, eventdata, handles)
set(handles.wnone,'Value',0);
set(handles.wrweight,'Value',0);
set(handles.wexact,'Value',1);
set(handles.ExactWeighting,'Enable','on');
% --------------------------------------------------------------------
% -------   CHECKBOX DEMO MODE   -------------------------------------
% --------------------------------------------------------------------
function demo_Callback(hObject, eventdata, handles)
if get(handles.demo,'Value')
    m=findobj(0,'Tag','rec3d');
    set(m,'Units','Pixel');pr3d=get(m,'Position');set(m,'Units','Normalized')
    h_demo=figure('Name','Reconstruction','Tag','demo1'); pdemo=get(h_demo,'Position');
    set(h_demo,...
        'Menubar','none',...
        'Toolbar','none',...
        'CloseRequestFcn','',...
        'Position',[(pr3d(1)+pr3d(3)+10) pr3d(2) pdemo(3) pr3d(4)],...
        'doublebuffer','on',...
        'Units','Normalized');    
    handles=InitDemoImage(handles);
    DemoRedrec(handles);
else
     m=findobj(0,'Tag','demo1');
     delete(m);
end    
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   CHECKBOX PARALLEL MODE   -------------------------------------
% --------------------------------------------------------------------
function parallel_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
% -------   BUTTON DO RECONSTRUCTION   -------------------------------
% --------------------------------------------------------------------
function Do_reconstruction_Callback(hObject, eventdata, handles)

if strcmp('DO reconstruction',get(hObject,'String'))
    set(hObject,'String','Stop');
else
    %set(hObject,'String','DO reconstruction','Enable','off');
    set(hObject,'String','DO reconstruction');
    drawnow;
    return
end% if

if isempty(get(handles.VolumeName,'String'))
    message=['ERROR!!!   Enter the file''s name of the reconstruction volume.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end
asd=size(get(handles.Refmarker,'String'));
for i=1:asd(1)
    if get(handles.Refmarker,'Value')==i
        ref_mark=i;
    end
end
size_rec(1)=str2num(get(handles.sizex,'String'));
size_rec(2)=str2num(get(handles.sizey,'String'));
size_rec(3)=str2num(get(handles.sizez,'String'));
offset_rec(1)=str2num(get(handles.offsetx,'String'));
offset_rec(2)=str2num(get(handles.offsety,'String'));
offset_rec(3)=str2num(get(handles.offsetz,'String'));
abs_center_rec(1)=str2num(get(handles.absx,'String'));
abs_center_rec(2)=str2num(get(handles.absy,'String'));
abs_center_rec(3)=str2num(get(handles.absz,'String'));
asd=size(get(handles.postbinning,'String'));
for i=1:asd(1)
    w=[0,1,2,3,4,5];
    if get(handles.postbinning,'Value')==i
        postbinning=w(i);
    end
end
asd=size(get(handles.prebinning,'String'));
for i=1:asd(1)
    w=[0,1,2,3,4,5];
    if get(handles.prebinning,'Value')==i
        prebinning=w(i);
    end
end
demo_mode=get(handles.demo,'Value');
weighting_size=1;
if get(handles.wnone,'Value')==1
    weighting='none';
elseif get(handles.wrweight,'Value')==1
    weighting='analytical';
elseif get(handles.wexact,'Value')==1
    weighting='exact';
    weighting_size=str2num(get(handles.ExactWeighting,'String'));
end
ref_image=str2num(get(handles.TSRefim,'String'));
lowpass_filter=str2num(get(handles.LowpassFilter,'String'));
handles.Param.AbsoluteCenter(1)=str2num(get(handles.absx,'String'));
handles.Param.AbsoluteCenter(2)=str2num(get(handles.absy,'String'));
handles.Param.AbsoluteCenter(3)=str2num(get(handles.absz,'String'));
handles.Param.ReconstructionName=get(handles.VolumeName,'String');
weightedfile_name=get(handles.WeightedFileName,'String');
wfo=get(handles.OnlyWeigthedFile,'Value');
parallel=get(handles.parallel,'Value');
volume = tom_reconstruction3d([handles.Param.PathName handles.Param.FileName],...
        [handles.Param.Ext],[handles.Param.MarkerFiler],ref_mark,...
        offset_rec,prebinning,postbinning,weighting,handles.Param.ReconstructionName,...
        lowpass_filter,size_rec,demo_mode,weighting_size, ...
        handles.Param.AbsoluteCenter,ref_image,weightedfile_name,wfo,parallel);
    uiwait(msgbox('Reconstruction done','Finished','help'));
    set(hObject,'String','DO reconstruction','Enable','on');

% --------------------------------------------------------------------
% -------   BUTTON EXIT   --------------------------------------------
% --------------------------------------------------------------------
function Exit_Callback(hObject, eventdata, handles)
message=['Do you really want to quit 3D Reconstruction?'];
Question=questdlg(message,'QUIT','Yes','No','Yes');
if strcmp(Question,'Yes')
    m=findobj('Tag','rec3d');
    if ~isempty(m)
        delete(m);
    end
    m=findobj('Tag','demo1');
    if ~isempty(m)
        delete(m);
    end
    
end
% ********************************************************************
% **********               OTHER FUNCTIONS                  **********
% ********************************************************************

% -------   Enable Panel Reconstruction 3D   -------
% Set the property 'Enable' to on or off
function EnablePanelRec3d(Status)
h_r3d=findobj(0,'Tag','rec3d');
handles=guidata(h_r3d);
set(handles.sizex,'Enable',Status);%editbox size X
set(handles.sizey,'Enable',Status);%editbox size Y
set(handles.sizez,'Enable',Status);%editbox size Z
set(handles.offsetx,'Enable',Status);%editbox offset X
set(handles.offsety,'Enable',Status);%editbox offset Y
set(handles.offsetz,'Enable',Status);%editbox offset Z
set(handles.absx,'Enable',Status);%editbox absolute center X
set(handles.absy,'Enable',Status);%editbox absolute center Y
set(handles.absz,'Enable',Status);%editbox absolute center Z
set(handles.VolumeName,'Enable',Status);%editbox volume name
set(handles.prebinning,'Enable',Status);%popupmenu pre binning
set(handles.postbinning,'Enable',Status);%popupmenu post binning
set(handles.LowpassFilter,'Enable',Status);%editbox filter number
set(handles.ExactWeighting,'Enable',Status);%editbox exactweighting
set(handles.demo,'Enable',Status);%checkbox demo
set(handles.OnlyWeigthedFile,'Enable',Status);%checkbox Weigthed file only
set(handles.label_x,'Enable',Status);%label X
set(handles.label_y,'Enable',Status);%label Y
set(handles.label_z,'Enable',Status);%label Z
set(handles.label_size,'Enable',Status);%label size
set(handles.label_offset,'Enable',Status);%label offset
set(handles.label_abs,'Enable',Status);%label absolute center
set(handles.label_name,'Enable',Status);%label volume name
set(handles.label_pre,'Enable',Status);%label pre binning
set(handles.label_post,'Enable',Status);%label post binning
set(handles.label_times1,'Enable',Status);%label times
set(handles.label_times2,'Enable',Status);%label times
set(handles.label_lowpass,'Enable',Status);%label lowpass filter
set(handles.Do_reconstruction,'Enable',Status);%Button do reconstruction
set(handles.browseVOL,'Enable',Status);%Button browse volume
set(handles.browseWEI,'Enable',Status);%Button browse weighted file
set(handles.wnone,'Enable',Status);%radiobutton none
set(handles.wrweight,'Enable',Status);%radiobutton r-weight
set(handles.wexact,'Enable',Status);%radiobutton exact

% -------   Enable Panel Marker file   -------
% Set the property 'Enable' to on or off
function EnablePanelMarkerFile(Status)
h_r3d=findobj(0,'Tag','rec3d');
handles=guidata(h_r3d);
set(handles.MarkerFileName,'Enable',Status); %editbox marker file 
set(handles.browseMF,'Enable',Status); %button browse

% -------   Initialize the menu   -------
function InitRec3d(handles)
set(handles.AlignmentFileName,'String','');%editbox alignment file
set(handles.TSPath,'String','');%editbox path
set(handles.TSFirst,'String','');%editbox from image
set(handles.TSLast,'String','');%editbox to image
set(handles.TSRefim,'String','');%editbox ref image
set(handles.Refmarker,'Value',1);%popupmenu ref marker
set(handles.MarkerFileName,'String','');%editbox marker file
set(handles.error,'String','');%editbox error
set(handles.WeightedFileName,'String','');%editbox Weighted File Name
set(handles.sizex,'String','');%editbox size X
set(handles.sizey,'String','');%editbox size Y
set(handles.sizez,'String','');%editbox size Z
set(handles.offsetx,'String','');%editbox offset X
set(handles.offsety,'String','');%editbox offset Y
set(handles.offsetz,'String','');%editbox offset Z
set(handles.absx,'String','');%editbox absolute center X
set(handles.absy,'String','');%editbox absolute center Y
set(handles.absz,'String','');%editbox absolute center Z
set(handles.VolumeName,'String','');%editbox volume name
set(handles.prebinning,'Value',3);%popupmenu pre binning
set(handles.postbinning,'Value',1);%popupmenu post binning
set(handles.LowpassFilter,'String','');%editbox filter number
set(handles.ExactWeighting,'String','');%editbox exactweighting
set(handles.demo,'Value',0);%checkbox demo
set(handles.OnlyWeigthedFile,'Value',0);%checkbox Weigthed file only
set(handles.label_pre,'ForegroundColor','k','UIContextMenu','');%label pre binning
set(handles.label_post,'ForegroundColor','k','UIContextMenu','');%label post binning
set(handles.Binning,'ForegroundColor','k','UIContextMenu','');%label binning
set(handles.Filter,'ForegroundColor','k','UIContextMenu','');%label filter
set(handles.Weighting,'ForegroundColor','k','UIContextMenu','');%label weighting
set(handles.wexact,'ForegroundColor','k','UIContextMenu','');%label exact weighting
set(handles.wrweight,'ForegroundColor','k','UIContextMenu','');%label r weighting
set(handles.label_size,'ForegroundColor','k','UIContextMenu','');%label size
set(handles.label_abs,'ForegroundColor','k','UIContextMenu','');%label absolute center
EnablePanelRec3d('off');
EnablePanelMarkerFile('off');

% -------   InitDemoImage   -------
function handles=InitDemoImage(handles)
myimage=[handles.Param.PathName handles.Param.FileName num2str(handles.Param.RefImage) handles.Param.Ext];
i=tom_emreadc(myimage);
%apply the alignment to the projection
%i=tom_rotate(i.Value,-handles.Param.Rotation_RefImage,'linear');
%i=tom_move(i,[-fix(handles.Param.Shift_RefImage(1)) -fix(handles.Param.Shift_RefImage(2))]);
%i=tom_emheader(i);
%end of hack
handles.DemoImage=i;
%delete plot
h_axe=findall(findobj(0,'Tag','demo1'),'Type','axe');
if ~isempty(h_axe)
    for rn=1:size(h_axe,1)
        delete(h_axe(rn));
    end
end
%
m=findobj('Tag','demo1');set(0,'CurrentFigure',m);
subplot(2,2,1);
set(gca,'position',[0.15 0.581098 0.327023 0.373902]);
m=findobj('Tag','demo1');set(0,'CurrentFigure',m);
tom_imagesc(i.Value,'noinfo');
title('IMAGE','FontWeight','bold','Units','normalized');

% -------   DemoRedrec   -------
function DemoRedrec(handles)
size_rec(1)=str2num(get(handles.sizex,'String'));
size_rec(2)=str2num(get(handles.sizey,'String'));
offset(1)=str2num(get(handles.offsetx,'String'));
offset(2)=str2num(get(handles.offsety,'String'));
a=get(handles.postbinning,'String');
b=get(handles.postbinning,'Value');
postbinning=str2num(a{b});
a=get(handles.prebinning,'String');
b=get(handles.prebinning,'Value');
prebinning=str2num(a{b});
a=round(size(handles.DemoImage.Value,1)./2)-((size_rec(1)./2)*(2^postbinning)*(2^prebinning));
a=a+offset(1);
b=round(size(handles.DemoImage.Value,1)./2)-((size_rec(2)./2)*(2^postbinning)*(2^prebinning));
b=b+offset(2);
w=round(size_rec(1))*(2^postbinning)*(2^prebinning);
h=round(size_rec(2))*(2^postbinning)*(2^prebinning);
if isempty(findobj(0,'Tag','rectan1'))
    m=findobj('Tag','demo1');set(0,'Currentfigure',m);
    rectan1=rectangle('Position',[a,b,w h],...
        'Tag','rectan1',...
        'LineWidth',2.5,...
        'EdgeColor',[1 0 0]);
else
    set(findobj(0,'Tag','rectan1'),'Position',[a,b,w h],...
        'Tag','rectan1',...
        'LineWidth',2.5,...
        'EdgeColor',[1 0 0]);
end




