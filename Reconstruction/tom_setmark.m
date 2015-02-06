function varargout = tom_setmark(varargin)
%TOM_SETMARK is a tool used to make an alignment file
%
%   varargout = tom_setmark(varargin)
%
%  The user load his tilt series picture and by clicking with the mouse, 
%  he puts makers on the picture. He can save his work by clicking 
%  'Save project' (Creation of a text file with general information 
%  (*.alg)) or he can create a Marker file (Creation of a EM-format with 
%  coordinate of the  marker). The coordinates are stored in a matrix of 
%  dimensions 12 x no.of image x no. of markers. The most important 
%  thing about the matrix is:
%   -1st row: the tilt angle is stored
%   -2nd row: x coordinate of the marker
%   -3rd row: y coordinate of the marker
%   -4th row: deviation of the marker with the line indicated the tilt axis azimut
%   -5th row: x coordinate of the alignment shift vector of the individual marker
%   -6th row: y coordinate of the alignment shift vector of the individual marker
%   -7th row: average of all the x shift
%   -8th row: average of all the y shift   
%   -9th row: x deviation of the mean shift
%   -10th row: y deviation of the mean shift
%   -11th row: running projection number
%   -12th row: size of the projection image
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
%   ... = tom_setmark(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by WDN 09/12/02
%   updated by WDN 03/01/05
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

%   last change: display actual file number in tom_tilt_lines 
%   FF 09/02/2010
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_setmark_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_setmark_OutputFcn, ...
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
function tom_setmark_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
EnableAllButton('off');
set(handles.new_conf,'Enable','on');
set(handles.load_conf,'Enable','on');
set(handles.save_conf,'Enable','on');
set(findobj(0,'Tag','menu'),'Units','normalized',...
    'Position',[0.0039 0.0557 0.2734 0.9141]);
VisibleFrameCC('off');
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = tom_setmark_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --------------------------------------------------------------------
% -------   BUTTON NEW PROJECT   -------------------------------------
% --------------------------------------------------------------------
function varargout = new_conf_Callback(hObject, eventdata, handles)
uiwait(tom_load_tiltseries);
m=findobj('Tag','menu');
param=get(m,'Userdata');
switch param.newproj_cancel
case 'no' %when the button Cancel on tom_load_tiltseries is not clicked
    set(handles.menu_number_ima,'String',param.mylastnb);
    set(handles.menu_ref_ima,'String',param.image_ref); 
    message=['...' param.myfilename param.image_ref param.myext];
    set(handles.menu_name_ima,'String',message);%set the name of the working image
    set(handles.menu_wi,'String',param.image_ref);
    InitMenu(handles);% refresh tom_setmark Menu
    n=str2num(param.image_ref);
    set(handles.menu_wi_angle,'String',param.angle(n));
    Figw=figure('Units','normalized',... 
        'menubar','figure',...
        'Toolbar','figure',...
        'Tag','Figwo',...  
        'Name','working image',...
        'NumberTitle','off',...
        'Position',[0.6031 0.0566 0.3125 0.390625],...
        'DoubleBuffer','on',...
        'Units','normalized',...
        'BackingStore','off',...
        'CloseRequestFcn','',...
        'Visible','on'); 
    menu0 = uimenu (Figw,'label', 'Process');    	    
    uimenu (menu0,'label', 'Contrast','Callback','n=gca;tom_smk_kontrast(n)');
    colormap(gray);%set the colormap
    Figw_z=figure('Units','normalized',... 
        'menubar','figure',...
        'Toolbar','figure',...
        'Tag','Figwo_z',...  
        'Name','High mag working image',...
        'NumberTitle','off',...
        'Position',[0.6195 0.5264 0.2953 0.3125],...
        'DoubleBuffer','on',...
        'Units','normalized',...
        'BackingStore','off',...
        'CloseRequestFcn','',...
        'Visible','on');
    menu1 = uimenu (Figw_z,'label', 'Process');    	    
    uimenu (menu1,'label', 'Contrast','Callback','n=gca;tom_smk_kontrast(n)');
    colormap(gray);%set the colormap  
    Figw_z_boxinfo2 = uicontrol('Units','normalized',...
        'Style', 'text',... 
        'Tag','figw_z_boxinfo2',...
        'String','',...
        'FontName','default',...
        'FontWeight','bold',...
        'BackgroundColor',[0.8 0.8 0.8],...
        'ForegroundColor',[1 0.45 0],...
        'HorizontalAlignment','left',...        
        'Position',[0.00142466 0.928356 0.256849 0.0763699],...
        'Units','normalized',...
        'Visible','on');    
    Fig_prev =figure('Units','normalized',...
        'MenuBar','figure',...
        'toolbar','figure',...
        'Tag','Fig_previ',...  
        'Name','previous image',...
        'NumberTitle','off',...        
        'DoubleBuffer','on',...
        'Units','normalized',...
        'Position',[0.2836 0.0566 0.3125 0.390625],...
        'BackingStore','off',...
        'CloseRequestFcn','',...
        'Visible','off');
    colormap(gray);
    Fig_prev_z=figure('Units','normalized',...
        'menubar','figure',...
        'Toolbar','figure',...
        'Tag','Fig_previ_z',...  
        'Name','High mag previous image',...
        'NumberTitle','off',...
        'Position',[0.3008 0.5264 0.2953 0.3125],...
        'DoubleBuffer','on',...
        'Units','normalized',...
        'BackingStore','off',...
        'CloseRequestFcn','',...
        'Visible','off');
    colormap(gray);%set the colormap  
    %%%%%%%%%%  working image  %%%%%%%%%%%%%
    aa=[param.mypathname param.myfilename param.image_ref param.myext];
    handles.WorkingImage=tom_emreadc(aa);%handles.WorkingImage.Value=double(handles.WorkingImage.Value);   
    set(0,'CurrentFigure',Figw);
    tom_setmark_imagesc(handles.WorkingImage.Value);
    title('WORKING IMAGE','FontWeight','bold','Units','normalized');
    message=num2str(param.angle(n));
    set(Figw,'Name',['...' param.myfilename param.image_ref param.myext '    Angle: ' num2str(param.angle(n))]);
    Xz=(handles.WorkingImage.Header.Size(1)/2)-((handles.WorkingImage.Header.Size(1)*0.1953125)/2);
    Yz=(handles.WorkingImage.Header.Size(1)/2)-((handles.WorkingImage.Header.Size(1)*0.1953125)/2);
    if handles.WorkingImage.Header.Size(1)>256
        r1=rectangle('Position',[Xz,Yz,256,256],'Tag','r1');
    else
        r1=rectangle('Position',[Xz,Yz,handles.WorkingImage.Header.Size(1),handles.WorkingImage.Header.Size(1)],'Tag','r1');
    end
    rect1=get(r1,'Position');
    am=findobj('Tag','menu');param=get(am,'Userdata');   
    set(am,'UserData',param);
    %%%%%%%%%%  working image zoomed %%%%%%%%%%%%%
    set(0,'CurrentFigure',Figw_z);
    set(Figw_z,'Name',['...' param.myfilename param.image_ref param.myext '    Angle: ' num2str(param.angle(n))]);
    temp=handles.WorkingImage.Value(rect1(1):(rect1(1)+rect1(3)),rect1(2):(rect1(2)+rect1(4)));
    tom_setmark_imagesc(temp);%display the file
    title('WORKING IMAGE ZOOMED','FontWeight','bold','Units','normalized');  
    handles.rectangle=rect1;   
    handles.FilterName='bandpass filter';
    guidata(hObject, handles);
    menu_contrast_Callback(handles.menu_contrast,'',handles);
    menu_filter_Callback(handles.menu_filter, '', handles);
end

% --------------------------------------------------------------------
% -------   BUTTON LOAD PROJECT   ------------------------------------
% --------------------------------------------------------------------
function varargout = load_conf_Callback(hObject, eventdata, handles)
rectangle=[0 0 0 0];
[f, p] = uigetfile('*.alg', 'LOAD PROJECT');  
if isequal(f,0)|isequal(p,0)
    %nothing because Cancel is ckicked
else
    conf_file=[p f];
    s=textread(conf_file,'%s','whitespace','\n');
    text=char(s(1));
    switch text
        case 'Alignment file use by tom_setmarker when the button ''Load'' a project is clicked'
            InitMenu(handles);
            Figw=figure('Units','normalized',... 
                'menubar','figure',...
                'Toolbar','figure',...
                'Tag','Figwo',...  
                'Name','working image',...
                'NumberTitle','off',...
                'Position',[0.6031 0.0566 0.3125 0.390625],...
                'DoubleBuffer','on',...
                'Units','normalized',...
                'BackingStore','off',...
                'CloseRequestFcn','',...
                'Visible','on');
            colormap(gray);%set the colormap
            Figw_z=figure('Units','normalized',... 
                'menubar','figure',...
                'Toolbar','figure',...
                'Tag','Figwo_z',...  
                'Name','High mag working image',...
                'NumberTitle','off',...
                'Position',[0.6195 0.5264 0.2953 0.3125],...
                'DoubleBuffer','on',...
                'Units','normalized',...
                'BackingStore','off',...
                'CloseRequestFcn','',...
                'Visible','on');
            colormap(gray);%set the colormap  
            Figw_z_boxinfo2 = uicontrol('Units','normalized',...
                'Style', 'text',... 
                'Tag','figw_z_boxinfo2',...
                'String','',...
                'FontName','default',...
                'FontWeight','bold',...
                'BackgroundColor',[0.8 0.8 0.8],...
                'ForegroundColor',[1 0.45 0],...
                'HorizontalAlignment','left',...               
                'Position',[0.00142466 0.928356 0.256849 0.0763699],...
                'Units','normalized',...
                'Visible','on');   
            Fig_prev =figure('Units','normalized',...
                'MenuBar','figure',...
                'toolbar','figure',...
                'Tag','Fig_previ',...  
                'Name','previous image',...
                'NumberTitle','off',...
                'Position',[0.2836 0.0566 0.3125 0.390625],...
                'DoubleBuffer','on',...
                'Units','normalized',...
                'BackingStore','off',...
                'CloseRequestFcn','',...
                'Visible','on');
            colormap(gray);
            Fig_prev_z=figure('Units','normalized',...
                'menubar','figure',...
                'Toolbar','figure',...
                'Tag','Fig_previ_z',...  
                'Name','High mag previous image',...
                'NumberTitle','off',...
                'Position',[0.3008 0.5264 0.2953 0.3125],...
                'DoubleBuffer','on',...
                'Units','normalized',...
                'BackingStore','off',...
                'CloseRequestFcn','',...
                'Visible','on');
            colormap(gray);%set the colormap
            am=findobj('Tag','menu');param=get(am,'Userdata');                       
            param.mypathname=char(s(2));%text=fscanf(mm,'%s\n',1);param.mypathname=text;            
            param.myfilename=char(s(3));
            param.myfirstnb=char(s(4));
            param.mylastnb=char(s(5));
            param.myext=char(s(6));
            param.myfilemarker_default=char(s(7));
            param.myfilemarker=char(s(8));
            param.image_ref=char(s(9));
            if size(s,1)==15 %old alg format
                param.newproj_cancel=char(s(10));
                re=char(s(11));rectangle(1)=str2num(re);
                re=char(s(12));rectangle(2)=str2num(re);
                re=char(s(13));rectangle(3)=str2num(re);
                re=char(s(14));rectangle(4)=str2num(re);
            else %new format
                param.newproj_cancel='no';
                re=char(s(10));rectangle(1)=str2num(re);
                re=char(s(11));rectangle(2)=str2num(re);
                re=char(s(12));rectangle(3)=str2num(re);
                re=char(s(13));rectangle(4)=str2num(re);
            end
            param.rectangle=rectangle;
            handles.rectangle=rectangle;
            showima=(param.image_ref);
            ImageNr=str2num(showima);
            for i=1:str2num(param.mylastnb)
                aa=[param.mypathname param.myfilename num2str(i) param.myext];
                ah=tom_reademheader(aa);
                tilt(i)=ah.Header.Parameter(19)./1000.0;
            end
            param.Projection_size=ah.Header.Size(1);            
            param.angle=tilt;            
            %%%%%%%%%%  MatrixMark  %%%%%%%%%%%%%
            cc=param.myfilemarker;
            temp=tom_emreadc(cc);            
            if size(temp.Value,1)==10 % transform old marker file (10 rows) to new one (12 rows) 
                temp2=-1*ones(12,str2num(param.mylastnb),size(temp.Value,3));
                temp2(1:10,:,:)=temp.Value;
                aa=[param.mypathname param.myfilename '1' param.myext];
                ah=tom_reademheader(aa);
                temp2(12,:,:)=ah.Header.Size(1);
                for i=1:str2num(param.mylastnb)
                    temp2(11,i,:)=i;
                end
                temp.Value=temp2;
            end   
            param.Matrixmark=double(temp.Value);
            Nbm=size(param.Matrixmark,3);        
            set(handles.menu_activemk,'String',num2str(Nbm));
            param.stop_auto_proc=0;
            param.currentmark=Nbm;
            set(handles.menu_name_mf,'String',cc);
            set(handles.menu_number_ima,'String',param.mylastnb); 
            set(handles.menu_ref_ima,'String',param.image_ref);
            am=findobj('Tag','menu');set(am,'Userdata',param);
            %%%%%%%%%%  working image  %%%%%%%%%%%%%
            aa=[param.mypathname param.myfilename param.image_ref param.myext]; 
            handles.WorkingImage=tom_emreadc(aa);%handles.WorkingImage.Value=double(handles.WorkingImage.Value);
            handles=tom_refresh_ima(showima,'work',handles);
            %%%%%%%%%%  information image   %%%%%%%%%%%%%
            bb=[param.mypathname param.myfilename num2str(ImageNr-1) param.myext];
            handles.PresentationImage=tom_emreadc(bb);%handles.PresentationImage.Value=double(handles.PresentationImage.Value);
            %Hack fb
               % handles.PresentationImage.Value=tom_filter(double(handles.PresentationImage.Value),7);
            %end Hack fb
            
            handles=tom_refresh_ima(num2str(ImageNr-1),'prev',handles);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            handles.FilterName='bandpass filter';
            guidata(hObject, handles); 
            loop=1;            
            menu_contrast_Callback(handles.menu_contrast,'',handles);
            menu_filter_Callback(handles.menu_filter, '', handles);
        otherwise
            uiwait(errordlg('The file selected is not a file for alignment. Please check your file. Your file should start with: "Alignment file use by tom_setmarker when the button ''Load'' a project is clicked"','File Error'));
    end
    guidata(hObject, handles);
end
% --------------------------------------------------------------------
% -------   BUTTON SAVE PROJECT   ------------------------------------
% --------------------------------------------------------------------
function varargout = save_conf_Callback(hObject, eventdata, handles)
am=findobj('Tag','menu');param=get(am,'Userdata');
actdir=pwd;
cd (param.mypathname);
[myname, mypathname] = uiputfile('*.alg', 'SAVE PROJECT AS'); 
myproject=[mypathname myname];
if myproject(1)~=0 & myproject(2)~=0   
    if isempty(findstr(myproject,'.alg'))
        myproject=strcat(myproject,'.alg');
    end
    mf=get(handles.menu_name_mf,'String');
    if isempty(mf)
        tom_savemark(handles);
    else
        message=['This marker file ' mf ' already exists. Do you want to replace it?'];
        Question=questdlg(message,'Save marker file','Yes','No','No');
        if strcmp(Question,'Yes') 
            am=findobj('Tag','menu');param=get(am,'Userdata');
            tom_emwrite(param.myfilemarker,param.Matrixmark);
        elseif strcmp(Question,'No')
            asd=tom_savemark(handles);
            switch asd
                case 'no'
                    return;%break;
            end
        end
    end        
    mm=fopen(myproject,'w');
    am=findobj('Tag','menu');
    param=get(am,'Userdata');
    fprintf(mm,'%s\n','Alignment file use by tom_setmarker when the button ''Load'' a project is clicked'); 
    fprintf(mm,'%s\n',param.mypathname); 
    fprintf(mm,'%s\n',param.myfilename);  
    fprintf(mm,'%s\n',param.myfirstnb);
    fprintf(mm,'%s\n',param.mylastnb); 
    fprintf(mm,'%s\n',param.myext);
    fprintf(mm,'%s\n',param.myfilemarker_default);
    fprintf(mm,'%s\n',param.myfilemarker);
    fprintf(mm,'%s\n',param.image_ref);  
    %fprintf(mm,'%s\n',param.newproj_cancel);  
    fprintf(mm,'%d\n',[handles.rectangle(1) handles.rectangle(2) handles.rectangle(3) handles.rectangle(4)]);
    %fprintf(mm,'%d\n',param.Projection_size);    
    fclose(mm); 
end 
cd (actdir);
% --------------------------------------------------------------------
% -------   BUTTON PREVIOUS IMAGE   ----------------------------------
% --------------------------------------------------------------------
function varargout = menu_prev_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','');
am=findobj('Tag','menu');param=get(am,'Userdata');
m=findobj('Tag','Fig_previ');asd=get(m,'Visible');
switch asd
case 'on' %INFORMATION IMAGE  is visible
    if get(handles.menu_radioup,'Value')==1 %up
        Ima_wi=get(handles.menu_wi,'String');Ima_wi_num=str2num(Ima_wi)-1;Ima_wi=num2str(Ima_wi_num);
        Ima_pi_num=Ima_wi_num-1;Ima_pi=num2str(Ima_pi_num);
        if Ima_wi_num<str2num(param.myfirstnb) 
            %No action when Image>last picture or when Image<first picture
        else
            handles.WorkingImage=handles.PresentationImage;
            if Ima_pi_num<=0                                              
                handles.PresentationImage.Value=0;
            else                
                aa=[param.mypathname param.myfilename Ima_pi param.myext]; 
                handles.PresentationImage=tom_emreadc(aa);

            end
            %%%%%%%%%%  information image   %%%%%%%%%%%%%
            handles=tom_refresh_ima(Ima_pi,'prev',handles);
        end
    else %down
        Ima_wi=get(handles.menu_infoim,'String');
        if isempty(Ima_wi)
            return;
        end
        Ima_wi_num=str2num(Ima_wi);
        Ima_pi_num=Ima_wi_num+1;Ima_pi=num2str(Ima_pi_num);
        if Ima_wi_num>str2num(param.mylastnb)
            %No action when Image>last picture or when Image<first picture
        else
            handles.WorkingImage=handles.PresentationImage;
            if Ima_pi_num>=str2num(param.mylastnb)+1                                             
                handles.PresentationImage.Value=0;
            else                
                aa=[param.mypathname param.myfilename Ima_pi param.myext]; 
                handles.PresentationImage=tom_emreadc(aa);         
            end
            %%%%%%%%%%  information image   %%%%%%%%%%%%%
            handles=tom_refresh_ima(Ima_pi,'prev',handles);
        end        
    end 
case 'off'%INFORMATION IMAGE  is NOT visible
    if get(handles.menu_radioup,'Value')==1 %up
        Ima_wi=get(handles.menu_wi,'String');Ima_wi_num=str2num(Ima_wi)-1;Ima_wi=num2str(Ima_wi_num);
        if Ima_wi_num<str2num(param.myfirstnb)
            return;
        end
    else %down
        Ima_wi=get(handles.menu_wi,'String');Ima_wi_num=str2num(Ima_wi)+1;Ima_wi=num2str(Ima_wi_num);
        if Ima_wi_num>str2num(param.mylastnb)
            return;
        end        
    end
    aa=[param.mypathname param.myfilename Ima_wi param.myext]; 
    handles.WorkingImage=tom_emreadc(aa);  
end 
%%%%%%%%%%  working image  %%%%%%%%%%%%%
handles=tom_refresh_ima(Ima_wi,'work',handles);               
guidata(hObject, handles);
% --------------------------------------------------------------------
% -------   BUTTON NEXT IMAGE   --------------------------------------
% --------------------------------------------------------------------
function varargout = menu_next_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','');
am=findobj('Tag','menu');param=get(am,'Userdata');
showima=get(handles.menu_wi,'String');
m=findobj('Tag','Fig_previ');asd=get(m,'Visible');
switch asd
case 'on' %INFORMATION IMAGE  is visible
    if get(handles.menu_radioup,'Value')==1 %up
        Ima_wi=get(handles.menu_wi,'String');Ima_wi_num=str2num(Ima_wi)+1;Ima_wi=num2str(Ima_wi_num);
        Ima_pi_num=Ima_wi_num-1;Ima_pi=num2str(Ima_pi_num);
        if Ima_wi_num>str2num(param.mylastnb) 
            return;
        else
            handles.PresentationImage=handles.WorkingImage;
            aa=[param.mypathname param.myfilename Ima_wi param.myext]; 
            handles.WorkingImage=tom_emreadc(aa);
            %%%%%%%%%%  information image   %%%%%%%%%%%%%
            handles=tom_refresh_ima(Ima_pi,'prev',handles);            
        end
    else %down
        Ima_wi=get(handles.menu_wi,'String');Ima_wi_num=str2num(Ima_wi)-1;Ima_wi=num2str(Ima_wi_num);
        Ima_pi_num=Ima_wi_num+1;Ima_pi=num2str(Ima_pi_num);
        if Ima_wi_num<str2num(param.myfirstnb)
            return;
        else
            handles.PresentationImage=handles.WorkingImage;
            aa=[param.mypathname param.myfilename Ima_wi param.myext]; 
            handles.WorkingImage=tom_emreadc(aa);         
            %%%%%%%%%%  information image   %%%%%%%%%%%%%
            handles=tom_refresh_ima(Ima_pi,'prev',handles);            
        end        
    end 
case 'off'%INFORMATION IMAGE  is NOT visible
    if get(handles.menu_radioup,'Value')==1 %up
        Ima_wi=get(handles.menu_wi,'String');Ima_wi_num=str2num(Ima_wi)+1;Ima_wi=num2str(Ima_wi_num);
        if Ima_wi_num>str2num(param.mylastnb)
            return;
        end
    else %down
        Ima_wi=get(handles.menu_wi,'String');Ima_wi_num=str2num(Ima_wi)-1;Ima_wi=num2str(Ima_wi_num);
        if Ima_wi_num<str2num(param.myfirstnb)
            return;
        end
    end
    aa=[param.mypathname param.myfilename Ima_wi param.myext]; 
    handles.WorkingImage=tom_emreadc(aa);   
end 
%%%%%%%%%%  working image  %%%%%%%%%%%%%
handles=tom_refresh_ima(Ima_wi,'work',handles);               
guidata(hObject, handles);
% --------------------------------------------------------------------
% -------   BUTTON ZOOM - Change the positon of the square on W.I. ---
% --------------------------------------------------------------------
function varargout = menu_zoom_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','ZOOM');
am=findobj('Tag','menu');param=get(am,'Userdata');
m=findobj('Tag','Figwo');
set(0,'CurrentFigure',m);
SetZoomPointer;%change shape of pointer on WI
title('CLICK ON THE IMAGE','FontWeight','bold','Units','normalized','Color','r');
handle_image=findall(m,'Type','image');
SetButtonDownFcn_wi('tom_setmark(''MouseDown_zoom'',gcbo,'''')'); 
EnableAllButton('off');
set(handles.menu_autostop,'Enable','on');

% --------------------------------------------------------------------
% -------   BUTTON SHOW IMAGE ON WORKING I.   ------------------------
% --------------------------------------------------------------------
function varargout = menu_showworkima_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','');
m=findobj('Tag','menu');
param=get(m,'Userdata');
showima=get(handles.menu_showwi,'String');
ImageNr=str2num(showima);
if isempty(ImageNr)
    message=['ERROR!!!   Please enter the image number you want to see.'];
    uiwait(msgbox(message,'Error','error'));
else    
    if ImageNr>str2num(param.mylastnb) | ImageNr<str2num(param.myfirstnb)
        message=['ERROR!!!   You want to see the image n??? ' showima ' ???                           This image doesn''t exist. Please select an other one.'];
        uiwait(msgbox(message,'Error','error'));
    else
        %%%%%%%%%%  working image  %%%%%%%%%%%%%
        aa=[param.mypathname param.myfilename showima param.myext]; 
        handles.WorkingImage=tom_emreadc(aa);
        handles=tom_refresh_ima(showima,'work',handles);
        %%%%%%%%%%  information image   %%%%%%%%%%%%%
        m=findobj('Tag','Fig_previ');asd=get(m,'Visible');
        switch asd
        case 'on' %INFORMATION IMAGE  is visible 
            if ImageNr-1>0
                bb=[param.mypathname param.myfilename num2str(ImageNr-1) param.myext];
                handles.PresentationImage=tom_emreadc(bb);
                handles=tom_refresh_ima(num2str(ImageNr-1),'prev',handles);
            else
                handles=tom_refresh_ima(num2str(ImageNr-1),'prev',handles);
            end
            guidata(hObject, handles);
        end
    end
end
guidata(hObject, handles);
% --------------------------------------------------------------------
% -------   BUTTON SHOW IMAGE ON PRESENTATION I.   -------------------
% --------------------------------------------------------------------
function varargout = menu_showqwim_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','');
m=findobj('Tag','menu');
param=get(m,'Userdata');
showima=get(handles.menu_showii,'String');
ImageNr=str2num(showima);
if isempty(ImageNr)
    message=['ERROR!!!   Please enter the image number you want to see.'];
    uiwait(msgbox(message,'Error','error'));
else    
    if ImageNr>str2num(param.mylastnb) | ImageNr<str2num(param.myfirstnb)
        message=['ERROR!!!   You want to see the image n??? ' showima ' ??? This image doesn''t exist. Please select an other one.'];
        h=msgbox(message,'Error','error');
        m=get(h,'Position');
        set(h,'Position',[100 m(2) m(3) m(4)]);
        %uiwait;
    else
        %%%%%%%%%%  information image   %%%%%%%%%%%%%    
        bb=[param.mypathname param.myfilename showima param.myext];
        handles.PresentationImage=tom_emreadc(bb);
        handles=tom_refresh_ima(showima,'prev',handles);                   
    end
    set(handles.menu_showii,'String','');
end
% --------------------------------------------------------------------
% -------   RADIOBUTTON APPLY CONTRAST -------------------------------
% --------------------------------------------------------------------
function menu_contrast_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.menu_filter,'Value',0);
    EnableApplyFilter('off');
    EnableApplyContrast('on');    
    if get(handles.cont_auto,'Value')
        av=mean2(handles.WorkingImage.Value);
        standart=std2(handles.WorkingImage.Value);
        maxi=max(max(handles.WorkingImage.Value));
        mini=min(min(handles.WorkingImage.Value));
        set(handles.cont_std,'String',num2str(standart));
        set(handles.cont_mean,'String',num2str(av));
        set(handles.scale1,'String',num2str(mini));
        set(handles.scale2,'String',num2str(maxi));
        handles.Filter.mean=av;
        handles.Filter.std=standart;
        handles.Filter.min=mini;
        handles.Filter.max=maxi;
        guidata(hObject, handles);
    end
else
    EnableApplyContrast('off');
    set(handles.cont_auto,'Value',1);
    set(handles.cont_manual,'Value',0);
    set(handles.cont_mean,'String','' );
    set(handles.cont_std,'String','' );
    set(handles.scale1,...
        'Style','text',...
        'BackgroundColor',[0.878431 0.87451 0.890196],...
        'String','');
    set(handles.scale2,...
        'Style','text',...
        'BackgroundColor',[0.878431 0.87451 0.890196],...
        'String','');
    a=findall(findobj(0,'Tag','Figwo'),'Type','text');
    flag=1;
    for ii=1:size(a,1)
        switch get(a(ii),'String')
            case 'CLICK ON THE IMAGE'
                flag=0;
                break;
        end
    end
    if flag==1
        handles=tom_refresh_ima(get(handles.menu_wi,'String'),'work',handles);
        guidata(hObject, handles);
    end    
end

% --------------------------------------------------------------------
% -------   RADIOBUTTON AUTOMATIC CONTRAST ---------------------------
% --------------------------------------------------------------------
function cont_auto_Callback(hObject, eventdata, handles)
set(handles.cont_manual,'Value',0);
if get(hObject,'Value')==0
    set(hObject,'Value',1);
end
set(handles.scale1,...
    'String',handles.Filter.min,...
    'Style','text',...
    'BackgroundColor',[0.878431 0.87451 0.890196],...
    'Enable','on');
set(handles.scale2,...
    'String',handles.Filter.max,...
    'Style','text',...
    'BackgroundColor',[0.878431 0.87451 0.890196],...
    'Enable','on');

% --------------------------------------------------------------------
% -------   RADIOBUTTON MANUAL CONTRAST ------------------------------
% --------------------------------------------------------------------
function cont_manual_Callback(hObject, eventdata, handles)
set(handles.cont_auto,'Value',0)
if get(hObject,'Value')==0
    set(hObject,'Value',1);
end
set(handles.scale1,...
    'Style','edit',...
    'BackgroundColor',[1 1 1],...
    'Enable','on');
set(handles.scale2,...
    'Style','edit',...
    'BackgroundColor',[1 1 1],...
    'Enable','on');
% --------------------------------------------------------------------
% -------   BUTTON DO IT CONTRAST ------------------------------------
% --------------------------------------------------------------------
function menu_cont_doit_Callback(hObject, eventdata, handles)
handles=tom_refresh_ima(get(handles.menu_wi,'String'),'work',handles);
guidata(hObject, handles);

% --------------------------------------------------------------------
% -------   RADIOBUTTON APPLY FILTER ---------------------------------
% --------------------------------------------------------------------
function menu_filter_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.menu_contrast,'Value',0);
    EnableApplyContrast('off');
    EnableApplyFilter('on');
else
    EnableApplyFilter('off');
    set(handles.filtertype,'Value',3);%filter type to 'tom_bandpass'
    set(handles.filter_val_low,'String','2');%editbox filter low
    set(handles.filter_val_hi,'String','70');%editbox filter hi
    handles.FilterName='bandpass filter';
    guidata(hObject, handles);
    a=findall(findobj(0,'Tag','Figwo'),'Type','text');
    flag=1;
    for ii=1:size(a,1)
        switch get(a(ii),'String')
            case 'CLICK ON THE IMAGE'
                flag=0;
                break;
        end
    end
    if flag==1
        handles=tom_refresh_ima(get(handles.menu_wi,'String'),'work',handles);
        guidata(hObject, handles);
    end
end

% --------------------------------------------------------------------
% -------   POPUPMENU  FILTER ----------------------------------------
% --------------------------------------------------------------------
function filtertype_Callback(hObject, eventdata, handles)
a=get(hObject,'Value');
b=get(hObject,'String');
filtername=b{a};
switch filtername
    case 'lowpass filter'
        set(handles.filter_low,'Enable','on');
        set(handles.filter_val_low,'Enable','on','String','50');
        set(handles.filter_hi,'Enable','off');
        set(handles.filter_val_hi,'Enable','off','String','');        
    case 'highpass filter'
        set(handles.filter_low,'Enable','off');
        set(handles.filter_val_low,'Enable','off','String','');
        set(handles.filter_hi,'Enable','on');
        set(handles.filter_val_hi,'Enable','on','String','20');
    case 'bandpass filter'
        set(handles.filter_low,'Enable','on');
        set(handles.filter_val_low,'Enable','on','String','2');
        set(handles.filter_hi,'Enable','on');
        set(handles.filter_val_hi,'Enable','on','String','70');
end
handles.FilterName=filtername;
guidata(hObject, handles);

% --------------------------------------------------------------------
% -------   BUTTON DO IT FILTER --------------------------------------
% --------------------------------------------------------------------
function menu_filt_doit_Callback(hObject, eventdata, handles)
handles=tom_refresh_ima(get(handles.menu_wi,'String'),'work',handles);
guidata(hObject, handles);
% --------------------------------------------------------------------
% -------   BUTTON ACTIVE -set an existing active marker -------------
% --------------------------------------------------------------------
function menu_active_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','');
am=findobj('Tag','menu');param=get(am,'Userdata');
actmk=get(handles.menu_setmk,'String');
Nbm=size(param.Matrixmark,3);
if str2num(actmk)<=Nbm    
    set(handles.menu_activemk,'String',actmk);
    set(handles.menu_setmk,'String','');
    param.currentmark=str2num(actmk);
    set(am,'UserData',param);
    %%%%%%%%%%  working image  %%%%%%%%%%%%%
    handles=tom_refresh_ima(get(handles.menu_wi,'String'),'work',handles);
    %%%%%%%%%%  informafion image  %%%%%%%%%%%%% 
    handles=tom_refresh_ima(get(handles.menu_infoim,'String'),'prev',handles);
else
    message=['ERROR!!!   This marker doesn''t exist.'];
    uiwait(msgbox(message,'Error','error'));
    set(handles.menu_setmk,'String','');
end   
% --------------------------------------------------------------------
% -------   BUTTON NEW MARKER -increase param.Matrixmark   -----------
% --------------------------------------------------------------------
function menu_newmk_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','');
am=findobj('Tag','menu');param=get(am,'Userdata');
Nbm=size(param.Matrixmark,3)+1;
param.Matrixmark(:,:,Nbm)=0;
param.Matrixmark(2:3,:,Nbm)=-1;
param.Matrixmark(11:12,:,Nbm)=param.Matrixmark(11:12,:,1);
param.Matrixmark(1,:,Nbm)=param.Matrixmark(1,:,1);
actmk=num2str(Nbm);
set(handles.menu_activemk,'String',actmk);
param.currentmark=str2num(actmk);
set(am,'UserData',param);
set(handles.menu_setmk,'String','');
% --------------------------------------------------------------------
% -------   BUTTON ADD A MARKER   ------------------------------------
% --------------------------------------------------------------------
function varargout = marker_add_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','ADD MARKER');
SetClickPointer;
SetZoomPointer;
EnableAllButton('off');
set(handles.menu_autostop,'Enable','on');
SetButtonDownFcn_wi('tom_setmark(''MouseDown_WI'',gcbo,''S'')'); 
SetButtonDownFcn_wiz('tom_setmark(''MouseDown_AddWIZ'',gcbo)');

% --------------------------------------------------------------------
% -------   BUTTON DELETE A MARKER   ---------------------------------
% --------------------------------------------------------------------
function varargout = marker_del_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','DELETE MARKER');
set(findobj(0,'Tag','Figwo_z'),'Pointer','fullcrosshair');
SetButtonDownFcn_wiz('tom_setmark(''MouseDown_DelMK'',gcbo)');
SetZoomPointer;%change shape of pointer on WI
SetButtonDownFcn_wi('tom_setmark(''MouseDown_zoom'',gcbo,''D'')');  
EnableAllButton('off');
set(handles.menu_autostop,'Enable','on');
msgbox('Please, select the marker you want to delete','Delete marker','help');

% --------------------------------------------------------------------
% -------   BUTTON START -Automatic procedure to add markers   -------
% --------------------------------------------------------------------
function varargout = menu_autostart_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','AUTOMATIC');
SetClickPointer;
am=findobj('Tag','menu');param=get(am,'Userdata');
showima=get(handles.menu_wi,'String');ImageNr=str2num(showima);j=ImageNr;
EnableAllButton('off');
set(handles.menu_autostop,'Enable','on');
if get(handles.menu_radioup,'Value')==1 %up
    j=ImageNr+1;
    if j>size(param.Matrixmark,2)
        message=('Impossible to start the automatic procedure. Last projection reached');
        msgbox(message,'Error','warn');
        menu_autostop_Callback(handles.menu_autostop,'', handles);
        return;
    end
    if param.currentmark>1%More than one mark
        pt=prediction(param,j,'up');
    else%less than one mark
        pt(1)=param.Matrixmark(2,j-1,param.currentmark);% extract x from Matrimark.So no prediction
        pt(2)=param.Matrixmark(3,j-1,param.currentmark);% extract y from Matrimark.So no prediction
    end
    Xz= pt(1)-(handles.rectangle(3)/2);Yz=pt(2)-(handles.rectangle(4)/2);%pt is the coordinate of the prediction, Xz Yz coordinates of rectangle with pt in the center
    aa=[param.mypathname param.myfilename num2str(j) param.myext];%read next picture
    handles.PresentationImage=handles.WorkingImage; %direction up
    if get(handles.menu_turbo,'Value')==0 %turbo is off
        SetZoomPointer;
        %%%%%%%%%%  working image  %%%%%%%%%%%%%
        handles.rectangle=[Xz Yz handles.rectangle(3) handles.rectangle(4)];
        handles.WorkingImage=tom_emreadc(aa);
        handles=tom_refresh_ima(num2str(j),'work',handles);
        %%%%%%%%%%  information image  %%%%%%%%%%%%%
        handles=tom_refresh_ima(num2str(j-1),'prev',handles);
    else%turbo is on        
        %%%%%%%%%%  information image  %%%%%%%%%%%%%
        m=findobj('Tag','Fig_previ');set(0,'CurrentFigure',m);%Let information and working image in this order
        tit=get(m,'Name');                                    %otherwise, there is a shift on I.image due to the rectangle
        if isempty(findstr(tit,'MODE TURBO'));
            bb=[param.mypathname param.myfilename num2str(j-1) param.myext];
            r1=round(param.Matrixmark(2,j-1,param.currentmark)-(handles.rectangle(3)/2));%x coordinate
            r2=round(param.Matrixmark(3,j-1,param.currentmark)-(handles.rectangle(4)/2));%y coordinate
            w=round(handles.rectangle(3));h=round(handles.rectangle(4));
            [handles.PresentationImage,handles.rectangle]=fillbymean(bb,r1,r2,w,h);
        end
        handles=tom_turbo(num2str(j-1),'prev',handles);
        %%%%%%%%%%  working image  %%%%%%%%%%%%%
        handles.rectangle=[Xz Yz handles.rectangle(3) handles.rectangle(4)];
        r1=round(handles.rectangle(1));r2=round(handles.rectangle(2));
        w=round(handles.rectangle(3));h=round(handles.rectangle(4));
        [handles.WorkingImage,handles.rectangle]=fillbymean(aa,r1,r2,w,h);
        handles=tom_turbo(num2str(j),'work',handles);
    end
else %down
    j=ImageNr-1;
    if j<1
        message=('Impossible to start the automatic procedure. First projection reached');
        msgbox(message,'Error','warn');
        menu_autostop_Callback(handles.menu_autostop,'', handles);
        return;
    end    
    if param.currentmark>1%More than one mark
        pt=prediction(param,j,'down');
    else%less than one mark
        pt(1)=param.Matrixmark(2,j+1,param.currentmark);% extract x from Matrimark.So no prediction
        pt(2)=param.Matrixmark(3,j+1,param.currentmark);% extract y from Matrimark.So no prediction
    end
    Xz= pt(1)-(handles.rectangle(3)/2);Yz=pt(2)-(handles.rectangle(4)/2);%pt is the coordinate of the prediction, Xz Yz coordinates of rectangle with pt in the center
    aa=[param.mypathname param.myfilename num2str(j) param.myext];%read previous picture
    handles.PresentationImage=handles.WorkingImage; %direction down
    if get(handles.menu_turbo,'Value')==0 %turbo is off
        SetZoomPointer;
        %%%%%%%%%%  working image  %%%%%%%%%%%%%
        handles.rectangle=[Xz Yz handles.rectangle(3) handles.rectangle(4)];
        handles.WorkingImage=tom_emreadc(aa);
        handles=tom_refresh_ima(num2str(j),'work',handles);
        %%%%%%%%%%  information image  %%%%%%%%%%%%%
        handles=tom_refresh_ima(num2str(j+1),'prev',handles);
    else%turbo is on        
        %%%%%%%%%%  information image  %%%%%%%%%%%%%
        m=findobj('Tag','Fig_previ');set(0,'CurrentFigure',m);%Let information and working image in this order
        tit=get(m,'Name');                                    %otherwise, there is a shift on I.image due to the rectangle
        if isempty(findstr(tit,'MODE TURBO'));
            bb=[param.mypathname param.myfilename num2str(j+1) param.myext];
            r1=round(param.Matrixmark(2,j+1,param.currentmark)-(handles.rectangle(3)/2));%x coordinate
            r2=round(param.Matrixmark(3,j+1,param.currentmark)-(handles.rectangle(4)/2));%y coordinate
            w=round(handles.rectangle(3));h=round(handles.rectangle(4));
            [handles.PresentationImage,handles.rectangle]=fillbymean(bb,r1,r2,w,h);
        end
        handles=tom_turbo(num2str(j+1),'prev',handles);
        %%%%%%%%%%  working image  %%%%%%%%%%%%%
        handles.rectangle=[Xz Yz handles.rectangle(3) handles.rectangle(4)];
        r1=round(handles.rectangle(1));r2=round(handles.rectangle(2));
        w=round(handles.rectangle(3));h=round(handles.rectangle(4));
        [handles.WorkingImage,handles.rectangle]=fillbymean(aa,r1,r2,w,h);
        handles=tom_turbo(num2str(j),'work',handles);
    end
end
guidata(hObject, handles);
%handle_image=findall(findobj(0,'Tag','Figwo_z'),'Type','image');
%set(handle_image,'ButtonDownFcn','tom_setmark(''MouseDown_WIZ'',gcbo)');
%handle_image=findall(findobj(0,'Tag','Figwo'),'Type','image');
%set(handle_image,'ButtonDownFcn','tom_setmark(''MouseDown_WI'',gcbo,''A'')');
 
SetButtonDownFcn_wiz('tom_setmark(''MouseDown_WIZ'',gcbo)');
SetButtonDownFcn_wi('tom_setmark(''MouseDown_WI'',gcbo,''A'')');

% --------------------------------------------------------------------
% -------   BUTTON STOP -no action   ---------------------------------
% --------------------------------------------------------------------
function varargout = menu_autostop_Callback(hObject, eventdata, handles)
set(findobj('Tag','figw_z_boxinfo2'),'String','STOP');
set(findobj('Tag','Figwo'),'Pointer','arrow');
set(findobj('Tag','Figwo_z'),'Pointer','arrow');
if get(handles.menu_turbo,'Value')
    set(handles.menu_turbo,'Value',0);
    am=findobj('Tag','menu');param=get(am,'Userdata');
    handles=update(param,handles);
    guidata(hObject, handles);
end
showima=get(handles.menu_wi,'String');
SetButtonDownFcn_wi(''); %clear ButtonDownFcn
SetButtonDownFcn_wiz('');
EnableAllButton('on');
menu_contrast_Callback(handles.menu_contrast,'',handles);
menu_filter_Callback(handles.menu_filter, '', handles);
set(0,'CurrentFigure',findobj(0,'Tag','Figwo'));
title('WORKING IMAGE','FontWeight','bold','Units','normalized','Color','k');

% --------------------------------------------------------------------
% -------   RADIO UP   -----------------------------------------------
% --------------------------------------------------------------------
function varargout = menu_radioup_Callback(hObject, eventdata, handles)
set(handles.menu_radioup,'Value',1);
set(handles.menu_radiodown,'Value',0);
% --------------------------------------------------------------------
% -------   RADIO DOWN   ---------------------------------------------
% --------------------------------------------------------------------
function varargout = menu_radiodown_Callback(hObject, eventdata, handles)
set(handles.menu_radiodown,'Value',1);
set(handles.menu_radioup,'Value',0);
% --------------------------------------------------------------------
% -------   CHECKBOX CROSS CORRELATION   -----------------------------
% --------------------------------------------------------------------
function menu_cross_Callback(hObject, eventdata, handles)
if get(handles.menu_cross,'Value')
    EnableCrossCorrelation('off')
    set(findobj(0,'Tag','figw_z_boxinfo2'),'String','CORRELATION');
else
    EnableCrossCorrelation('on')
    menu_contrast_Callback(handles.menu_contrast,'',handles);
    menu_filter_Callback(handles.menu_filter, '', handles);    
    set(findobj(0,'Tag','figw_z_boxinfo2'),'String','');
end  
% --------------------------------------------------------------------
% -------   RADIO BUTTON UP   ----------------------------------------
% --------------------------------------------------------------------
function dir_up_Callback(hObject, eventdata, handles)
set(hObject,'Value',1);
set(handles.dir_down,'Value',0);
set(handles.dir_both,'Value',0);
% --------------------------------------------------------------------
% -------   RADIO BUTTON DOWN   --------------------------------------
% --------------------------------------------------------------------
function dir_down_Callback(hObject, eventdata, handles)
set(hObject,'Value',1);
set(handles.dir_up,'Value',0);
set(handles.dir_both,'Value',0);
% --------------------------------------------------------------------
% -------   RADIO BUTTON BOTH   --------------------------------------
% --------------------------------------------------------------------
function dir_both_Callback(hObject, eventdata, handles)
set(hObject,'Value',1);
set(handles.dir_down,'Value',0);
set(handles.dir_up,'Value',0);
% --------------------------------------------------------------------
% -------   BUTTON DO IT CROSS CORRELATION   -------------------------
% --------------------------------------------------------------------
function menu_cc_doit_Callback(hObject, eventdata, handles)
am=findobj('Tag','menu');param=get(am,'Userdata');
showima=get(handles.menu_wi,'String');ImageNr=str2num(showima);j=ImageNr;
for i=1:size(param.Matrixmark,2)
    if param.Matrixmark(2:3,i,1)==-1
        message='Marker 1 has to be in all projection. Can not process a cross correlation';
        msgbox(message,'Error','warn')
        return;
    end
end
if get(handles.dir_up,'Value')
    param.Matrixmark(2,ImageNr+1:str2num(param.mylastnb),param.currentmark)=-1; %delete previous CC values
    param.Matrixmark(3,ImageNr+1:str2num(param.mylastnb),param.currentmark)=-1;
    for i=ImageNr:str2num(param.mylastnb)-1
        [status,param]=cross_corr(param,handles,i,'up');
        if ~isempty(status)
            msgbox(status);
            resetcc(hObject, eventdata, handles); %reinitialize the button after cross-correlation
            return;
        end
    end
    message=('Alignment by cross correlation in upper direction done');
    msgbox(message);
    resetcc(hObject, eventdata, handles);
elseif get(handles.dir_down,'Value')
    param.Matrixmark(2,str2num(param.myfirstnb):ImageNr-1,param.currentmark)=-1; %delete previous CC values
    param.Matrixmark(3,str2num(param.myfirstnb):ImageNr-1,param.currentmark)=-1;
    for i=ImageNr:-1:str2num(param.myfirstnb)+1
        [status,param]=cross_corr(param,handles,i,'down');
        if ~isempty(status)
            msgbox(status);
            resetcc(hObject, eventdata, handles);
            %%%%%%%%%%  information image  %%%%%%%%%%%%%
            handles=tom_refresh_ima(num2str(ImageNr-1),'prev',handles);
            return;
        end
    end
    message=('Alignment by cross correlation in down direction done');
    msgbox(message);
    resetcc(hObject, eventdata, handles);
    %%%%%%%%%%  information image  %%%%%%%%%%%%%
    handles=tom_refresh_ima(num2str(ImageNr-1),'prev',handles);
elseif get(handles.dir_both,'Value')
    param.Matrixmark(2,ImageNr+1:str2num(param.mylastnb),param.currentmark)=-1; %delete previous CC values
    param.Matrixmark(3,ImageNr+1:str2num(param.mylastnb),param.currentmark)=-1;
    param.Matrixmark(2,str2num(param.myfirstnb):ImageNr-1,param.currentmark)=-1; %delete previous CC values
    param.Matrixmark(3,str2num(param.myfirstnb):ImageNr-1,param.currentmark)=-1;
    for i=ImageNr:str2num(param.mylastnb)-1
        [status1,param]=cross_corr(param,handles,i,'up');
        if ~isempty(status1)
            break;
        end
    end
    for i=ImageNr:-1:str2num(param.myfirstnb)+1
        [status2,param]=cross_corr(param,handles,i,'down');
        if ~isempty(status2)
            if ~isempty(status1)
                message=[status2 '. ' status1];
            else
                message=['Alignment by cross correlation in upper direction done. ' status2 ];
            end
            msgbox(message);
            resetcc(hObject, eventdata, handles);
            %%%%%%%%%%  information image  %%%%%%%%%%%%%
            handles=tom_refresh_ima(num2str(ImageNr-1),'prev',handles);
            return;
        end
    end
    if ~isempty(status1)
        message=(['Alignment by cross correlation in down direction done. ' status1])
    else
        message=('Alignment by cross correlation in upper and down direction done');
    end
    msgbox(message);
    resetcc(hObject, eventdata, handles);
    %%%%%%%%%%  information image  %%%%%%%%%%%%%
    handles=tom_refresh_ima(num2str(ImageNr-1),'prev',handles);
end


% --------------------------------------------------------------------
% -------   BUTTON SHOW MARKER ON PRESENTATION I.   ------------------
% --------------------------------------------------------------------
function varargout = menu_showmark_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','');
am=findobj('Tag','menu');param=get(am,'Userdata');
numima=get(handles.menu_no_ima,'String');
showima=str2num(numima);
numm=str2num(get(handles.menu_no_mark,'String'));
Nbm=size(param.Matrixmark,3);
if isempty(numm)
    message=['ERROR!!!   Please enter the marker number you want to show.'];
    uiwait(msgbox(message,'Error','error'));
else
    if numm>Nbm | numm==0
        message=['ERROR!!!   This mark doesn''t exist. Please select an other one.'];
        uiwait(msgbox(message,'Error','error'));
    else
        m=findobj('Tag','Fig_previ_z');set(0,'CurrentFigure',m);
        pt(1)=param.Matrixmark(2,showima,numm);
        pt(2)=param.Matrixmark(3,showima,numm);        
        rect2=handles.rectangle;
        Xz= pt(1)-(rect2(3)/2);Yz=pt(2)-(rect2(4)/2);
        m=findobj('Tag','r2');set(m,'Position',[Xz Yz rect2(3) rect2(4)]);
        handles.rectangle=get(m,'Position');
        handles=tom_limz('prev',handles);rect2=handles.rectangle;
        temp=handles.PresentationImage.Value(round(rect2(1)):(round(rect2(1)+rect2(3))),round(rect2(2)):round((rect2(2)+rect2(4))));
        tom_setmark_imagesc(temp);%display the file                
        handles=tom_refresh_mark(numima,param.Matrixmark,'prevz',handles);
    end
end
% --------------------------------------------------------------------
% -------   BUTTON MOVIE ON PRESENTATION I.   ------------------------
% --------------------------------------------------------------------
function varargout = menu_movie_Callback(hObject, eventdata, handles)
am=findobj('Tag','menu');param=get(am,'Userdata');
view_movie=str2num(get(handles.menu_no_mark,'String'));
if isempty(view_movie)
    message=['ERROR!!!   Please enter the marker number you want to see the movie.'];
    uiwait(msgbox(message,'Error','error'));
else    
    if view_movie>size(param.Matrixmark,3)| view_movie<=0
        message=['ERROR!!!   This mark doesn''t exist. Please select an other one.'];
        uiwait(msgbox(message,'Error','error'));       
    else
        m=findobj('Tag','Fig_previ_z');
        set(0,'CurrentFigure',m);
        set(m,'Name',['Movie centered on Marker Nr.: ' num2str(view_movie)]); 
        for index=1:size(param.Matrixmark,2)
            Radius=10;
            w=256;h=256;
            u=w/2;v=h/2;
            uu = [u u u u-Radius u+Radius];vv = [v-Radius v+Radius v v v] ;      
            set(0,'CurrentFigure',m);
            r1=param.Matrixmark(2,index,view_movie);
            r2=param.Matrixmark(3,index,view_movie);
            bb=[param.mypathname param.myfilename num2str(index) param.myext];
            if r1~=-1 & r2~=-1  
                r1=double(r1)-256./2;r1=round(r1);
                r2=double(r2)-256./2;r2=round(r2);
                if index==18
                    rrr='';
                end
                title(['MOVIE: image nr. ' num2str(index)] ,'FontWeight','bold');
                [img_tmp,movie_rectangle]=fillbymean(bb,r1,r2,w,h);
                tom_setmark_imagesc(img_tmp.Value);
                title(['MOVIE: image nr. ' num2str(index)] ,'FontWeight','bold');
                line(uu,vv,'LineWidth',1,'Color',[1 0 0])%light red 
                drawnow;
           end
        end
    end
end

% --------------------------------------------------------------------
% -------   BUTTON LOAD MARKER FILE   --------------------------------
% --------------------------------------------------------------------
function varargout = menu_loadmarkfile_Callback(hObject, eventdata, handles)
am=findobj('Tag','menu');param=get(am,'Userdata');
actdir=pwd;
cd (param.mypathname);
[myfile, mypathname]=uigetfile({'*.em'},'Load a marker file');
if isequal(myfile,0)|isequal(mypathname,0)
    %nothing because Cancel is ckicked
else
	temp=tom_emread([mypathname myfile ]);
    if size(temp.Value,1)==10 % transform old marker file (10 rows) to new one (12 rows) 
        temp2=-1*ones(12,str2num(param.mylastnb),size(temp.Value,3));
        temp2(1:10,:,:)=temp.Value;
        aa=[param.mypathname param.myfilename '1' param.myext];
        ah=tom_reademheader(aa);
        temp2(12,:,:)=ah.Header.Size(1);
        for i=1:str2num(param.mylastnb)
            temp2(11,i,:)=i;
        end
        temp.Value=temp2;
    end
	param.Matrixmark=temp.Value;
    param.myfilemarker=[mypathname myfile];
	Nbm=size(param.Matrixmark,3);
    param.currentmark=Nbm;
	set(handles.menu_activemk,'String',num2str(Nbm))
	set(handles.menu_name_mf,'String',[mypathname myfile ]);
	set(am,'Userdata',param);
	%%%%%%%%%%  working image  %%%%%%%%%%%%%
	showima=get(handles.menu_wi,'String');
	handles=tom_refresh_ima(showima,'work',handles);
	%%%%%%%%%%  information image   %%%%%%%%%%%%%
	m=findobj('Tag','Fig_previ');asd=get(m,'Visible');
	switch asd
	case 'off'
        set(m,'Visible','on');
        ImageNr=str2num(get(handles.menu_wi,'String'));
        bb=[param.mypathname param.myfilename num2str(ImageNr-1) param.myext];
        handles.PresentationImage=tom_emreadc(bb);
        m=findobj('Tag','Fig_previ_z');
        set(m,'Visible','on'); 
        handles=tom_refresh_ima(num2str(ImageNr-1),'prev',handles);
        guidata(hObject, handles);
	otherwise
        showima=get(handles.menu_infoim,'String');
        handles=tom_refresh_ima(showima,'prev',handles);    
	end
	cd (actdir);
end
% --------------------------------------------------------------------
% -------   BUTTON SAVE MARKER FILE   --------------------------------
% --------------------------------------------------------------------
function varargout = menu_savemarkfile_Callback(hObject, eventdata, handles)
tom_savemark(handles);

% --------------------------------------------------------------------
% -------   BUTTON VIEW MARKER FILE   --------------------------------
% --------------------------------------------------------------------
function menu_view_Callback(hObject, eventdata, handles)
am=findobj('Tag','menu');param=get(am,'Userdata');
if get(handles.AllOrOneGM,'Value')
    param.Matrixmark
else
    disp(['*****  Marker Number ' get(handles.whichGM,'String') ' *****' ]);
    gom=str2num(get(handles.whichGM,'String'));
    disp(param.Matrixmark(:,:,gom))
end
assignin('base','Matrixmark',param.Matrixmark);
% --------------------------------------------------------------------
% -------   CHECKBOX ALL   -------------------------------------------
% --------------------------------------------------------------------
function AllOrOneGM_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.whichGM,'Enable','off');
else
    set(handles.whichGM,'Enable','on');
end
% --------------------------------------------------------------------
% -------   BUTTON TILT LINE ON PRESENTATION I.   --------------------
% --------------------------------------------------------------------
function varargout = menu_tilt_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','');
am=findobj('Tag','menu');param=get(am,'Userdata');
ref=get(handles.menu_tilt1,'String');
pkt=get(handles.menu_tilt2,'String');
Nbm=size(param.Matrixmark,3);
switch ref
    case''
        message=['ERROR!!!  No mark selected. Please enter one in the 1st case.'];
        uiwait(msgbox(message,'Error','error'));         
    otherwise
        f=str2num(ref);
        if f>=1&f<=Nbm            
            switch pkt
                case''
                    message=['ERROR!!!  No mark selected. Please enter one in the 2nd case.'];
                    uiwait(msgbox(message,'Error','error'));
                otherwise
                    g=str2num(pkt);
                    if g>=1&g<=Nbm                
                        m=findobj('Tag','Fig_tilt');
                        if isempty(m)
                            m =figure('Units','normalized',...
                            'MenuBar','none',...
                            'toolbar','figure',...
                            'Tag','Fig_tilt',...  
                            'Name','TILT LINES',...
                            'NumberTitle','off',...
                            'Position',[0.4 0.303711 0.3125 0.390625],...
                            'DoubleBuffer','on',...
                            'Units','normalized',...
                            'BackingStore','off',...
                            'Visible','on');          
                        end
                        set(0,'CurrentFigure',m);
                        clf;                        
                        pkt=str2num(pkt);ref=str2num(ref);
                        tom_tilt_lines(param.Matrixmark,ref,pkt);
                    else
                        message=['ERROR!!!  The mark number ' pkt ' doesn''t exist. Please enter an other one.'];
                        uiwait(msgbox(message,'Error','error'));
                    end                    
            end
        else
            message=['ERROR!!!  The mark number ' ref ' doesn''t exist. Please enter an other one.'];
            uiwait(msgbox(message,'Error','error'));
        end
end                      
% --------------------------------------------------------------------
% -------   BUTTON ALIGNMENT   ---------------------------------------
% --------------------------------------------------------------------                
function varargout = menu_align_Callback(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','');
am=findobj('Tag','menu');
param=get(am,'Userdata');
m=findobj('Tag','menu_definemarker_tilt');
f=get(m,'String');
mk=str2num(f);
Nbm=size(param.Matrixmark,3);
switch f
    case ''
        message=['ERROR!!!  No mark selected. Please enter one.'];
        uiwait(msgbox(message,'Error','error'));
    otherwise
        if Nbm<=1
            message=['ERROR!!!  You need at least 2 markers to make a alignment. So please put enought marker before trying to make a alignment'];
            uiwait(msgbox(message,'Error','error'));
        elseif param.Matrixmark(2:3,:,Nbm)==-1
            act_mk=get(handles.menu_activemk,'String');
            message=['ERROR!!!  the active marker is number ' act_mk ' and there is no marker ' act_mk ' clicked in all the images. So please click one. '];
            uiwait(msgbox(message,'Error','error'));                      
        else             
            if mk>=1&mk<=Nbm
                [param.Matrixmark,psi,error,x,y,z]=tom_alignment3d(param.Matrixmark,mk);               
                m=findobj('Tag','menu_error');set(m,'String',num2str(error));
                param.Matrixmark(7,3,:)=param.Projection_size;%needed by EM program for rec                                                                
                set(am,'Userdata',param);
            else
                message=['ERROR!!!  The mark number ' f ' doesn''t exist. Please enter an other one.'];
                uiwait(msgbox(message,'Error','error'));
            end
        end
end
m=findobj('Tag','menu');
param=get(m,'Userdata');
temp=param.Matrixmark;
temp=tom_emheader(temp);
tom_emwrite(param.myfilemarker_default,temp);
% --------------------------------------------------------------------
% -------   BUTTON 3D RECONSTRUTION   --------------------------------
% --------------------------------------------------------------------
function varargout = menu_3d_Callback(hObject, eventdata, handles)
tom_rec3d;
% --------------------------------------------------------------------
% -------   BUTTON QUIT PROGRAMME   ----------------------------------
% --------------------------------------------------------------------
function varargout = menu_quit_Callback(hObject, eventdata, handles)
message=['Do you really want to quit tom_setmark?'];
Question=questdlg(message,'QUIT','Yes','No','Yes');
if strcmp(Question,'Yes')
    m=findobj('Tag','Figwo');
    if ~isempty(m)
        delete(m);
    end
    m=findobj('Tag','Figwo_z');
    if ~isempty(m)    
        delete(m);
    end
    m=findobj('Tag','Fig_previ');
    if ~isempty(m)    
        delete(m);
    end
    m=findobj('Tag','Fig_previ_z');
    if ~isempty(m)    
        delete(m);
    end 
    m=findobj('Tag','menu');
    if ~isempty(m)    
        delete(m);
    end 
end
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% ----------- Function needs by tom_setmark -------------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------

% ----------- Function update -----------
function handles=update(param,handles)
% UPDATE is used to refresh the working image and 
% information image when the mode turbo is finished 
%%%%%%%%%%  working image  %%%%%%%%%%%%%
wi=get(handles.menu_wi,'String');
aa=[param.mypathname param.myfilename wi param.myext]; 
handles.WorkingImage=tom_emreadc(aa);
handles=tom_refresh_ima(wi,'work',handles);
%%%%%%%%%%  information image  %%%%%%%%%%%%%
ii=get(handles.menu_infoim,'String');
bb=[param.mypathname param.myfilename ii param.myext];
handles.PresentationImage=tom_emreadc(bb);
handles=tom_refresh_ima(ii,'prev',handles);                        

% ----------- Function tom_setmark_imagesc -----------
function tom_setmark_imagesc(in,parameter)
% TOM_SETMARK_IMAGESC is the same as tom_imagesc but easier
%  Input:  -in: Array of data of the image
%          -parameter: 'fixed' to have one pixel on screen equal to one
%                       pixel on file
%  Output:
%           - 
in = tom_filter(in,3,'circ');
in_red=imresize(in,.1);
[meanv max min std]=tom_devinternal(in_red);
if (meanv-3*std)>=(meanv+3*std)
    imagesc(in');
else
    imagesc(in',[meanv-3*std meanv+3*std]);colormap gray;axis image;
end;
colormap gray;   
if nargin==2
    switch parameter
        case 'fixed'
            set(gca,'Units','pixels');
            pp=get(gca,'Position');sf=size(in);            
            set(gca,'Position',[pp(1) pp(2) sf(1) sf(2)]);
        otherwise
            axis image; axis ij; colormap gray; %nothing changed, as nargin=1
    end
elseif nargin==1
    axis image; axis ij; %colormap gray;
end  

% ----------- Function tom_devinternal -----------
function [a,b,c,d,e]=tom_devinternal(A);
%TOM_DEVINTERNAL is only used  by TOM_SETMARK_IMAGESC 
A=double(A);
[s1,s2,s3]=size(A);
a=sum(sum(sum(A)))/(s1*s2*s3);
b=max(max(max(A)));
c=min(min(min(A)));
d=std2(A);
e=d^2;
% ----------- Function setmark_circle -----------
function [X, Y] = setmark_circle(w,r,n)
%SETMARK_CIRCLE is used to calculate the coordinate of a circle.
%Input
%-w: Center of the circle. W must be a complex number as w=x + yi 
%    (x and y are the coordinate)
%-r: Radius of the circle
%-n: The width of the line         
%Output:
%-X: it is a matrix of coordinate to draw the circle
%-Y: it is a matrix of coordinate to draw the circle
 w1 = real(w);
 w2 = imag(w);
        for k = 1:n
           t = k*pi/n;
           X(k) = w1 + r*cos(t);
           Y(k) = w2 + r*sin(t);
           X(n+k) = w1 - r*cos(t);
           Y(n+k) = w2 - r*sin(t);
           X(2*n+1) = X(1);
           Y(2*n+1) = Y(1);
        end

% ----------- Function tom_drawmark -----------
function tom_drawmark(x,y,i)
% TOM_DRAWMARK draw a circle and a cross to represent the marker 
%Input:
% -x: coordinate x
% -y: coordinate y
% -i: text to print (must be a string)
hold on;
Center= x + y*sqrt(-1);
Radius = 10;Gridpt = 100;
[u,v]=setmark_circle(Center,Radius,Gridpt);
uu = [x x x x-Radius x+Radius];
vv = [y-Radius y+Radius y y y];
info=get(findobj(0,'Tag','figw_z_boxinfo2'),'String');
switch info %set ButtonDownFcn on the line and text. Otherwise no zoom when click on them
    case 'AUTOMATIC'
    line(u,v,'LineWidth',1,'Color',[1 0 0]);%red dark
    L=line(uu,vv,'LineWidth',1,'color',[1 0 0],'ButtonDownFcn','tom_setmark(''MouseDown_WI'',gcbo,''A'')');%red dark
    T=text(x+25,y+3,i,'FontWeight','bold','Color',[1 0.75 0],'Fontsize',20,'ButtonDownFcn','tom_setmark(''MouseDown_WI'',gcbo,''A'')');%orange
    case 'ADD MARKER'
    line(u,v,'LineWidth',1,'Color',[1 0 0]);%red dark    
    L=line(uu,vv,'LineWidth',1,'color',[1 0 0],'ButtonDownFcn','tom_setmark(''MouseDown_WI'',gcbo,''S'')');%red dark
    T=text(x+25,y+3,i,'FontWeight','bold','Color',[1 0.75 0],'Fontsize',20,'ButtonDownFcn','tom_setmark(''MouseDown_WI'',gcbo,''S'')');%orange
    otherwise
    line(u,v,'LineWidth',1,'Color',[1 0 0]);%red dark 
    L=line(uu,vv,'LineWidth',1,'color',[1 0 0]);%red dark
    T=text(x+25,y+3,i,'FontWeight','bold','Color',[1 0.75 0],'Fontsize',20);%orange        
end 
hold off;

% ----------- Function tom_delmark -----------
function tom_delmark(Im,W,Ma,Cm,h)
% TOM_DELMARK delete a mark on the current picture or in all pictures
%Input:
%-Im: The image number to refresh. Must be a string
%-W: If 'work' then refresh the working image
%    if 'prev' then refresh the information image
%-Ma: the value of param.Matrixmark
%-h: all the handles
%-Cm: the mark to delete
am=findobj('Tag','menu');param=get(am,'Userdata');
switch W
case 'Current'    
    Ma(2,Im,Cm)=-1;
    Ma(3,Im,Cm)=-1;
    if Ma(2:3,:,Cm)==-1 
        Nbm=size(param.Matrixmark,3);
        if Nbm~=1                
            Ma(:,:,Cm)=[];
            set(h.menu_activemk,'String',num2str(Nbm-1));
        end
    end   
    showima=num2str(Im);                
    param.Matrixmark=Ma;
    set(am,'Userdata',param);
    handles=tom_refresh_ima(showima,'work',h);
    temp=Ma;
    temp=tom_emheader(temp);
    tom_emwrite(param.myfilemarker_default,temp);
case 'All'
    Nbm=size(Ma,3);
    if Nbm==1
        Ma(2:3,:,1)=-1;
        set(h.menu_activemk,'String',num2str(Nbm));
    else
        for ii=Cm:Nbm
            if ii~=Nbm
                Ma(2:3,:,ii)=Ma(2:3,:,ii+1);
            end
        end
        Ma(:,:,Nbm)=[];
        set(h.menu_activemk,'String',num2str(Nbm-1));
    end            
    am=findobj('Tag','menu');param=get(am,'Userdata');
    param.Matrixmark=Ma;
    param.currentmark=str2num(get(h.menu_activemk,'String'));
    set(am,'Userdata',param);
    showima=num2str(Im);
    handles=tom_refresh_ima(showima,'work',h);
    showima=num2str(Im-1); 
    handles=tom_refresh_ima(showima,'prev',h);
    temp=param.Matrixmark;
    temp=tom_emheader(temp);
    tom_emwrite(param.myfilemarker_default,temp);
end

% ----------- Function tom_savemark -----------
function fiex=tom_savemark(h)
% TOM_SAVEMARK save param.Matrixmark into a file   
%   Syntax: fiex=tom_savemark(h)
%       Input:
%           h: all the handles
%       Output:
%           fiex: If 'yes' then write into a file
%                 If 'no' then cancel pressed 
am=findobj('Tag','menu');param=get(am,'Userdata');
temp=param.Matrixmark;
temp=tom_emheader(temp);
actdir1=pwd;
cd (param.mypathname);
[myname, mypathname] = uiputfile('*.em', 'SAVE YOUR MARKER FILE AS');
myfile=[mypathname myname];
if myfile(1)~=0 & myfile(2)~=0 %myfile= 0 0 when 'cancel' is clicked
    if isempty(findstr(myfile,'.em'))
        myfile=strcat(myfile,'.em');
    end
    tom_emwrite(myfile,temp);
    set(h.menu_name_mf,'String',myfile);
    param.myfilemarker=myfile;
    fiex='yes';
    set(am,'Userdata',param);
else
    fiex='no';
end
cd (actdir1);
% ----------- Function tom_limz -----------
function h=tom_limz(W,h)
% TOM_LIMZ check if the square is in the picture. If the square 
%   is off limit, then the function place correctly the square.   
%   Syntax: tom_limz(W,handles)
%       Input:
%           W: Where to replace
%               - If 'work' replace in the working image
%               - If 'prev' replace in the information image
%           Rec: rectangle
%           h: all the handles
%       Output:
%           h: all the handles
ss=size(h.WorkingImage.Value);
if h.rectangle(1)<=0
    h.rectangle(1)=1;               
end
if h.rectangle(2)<=0;
    h.rectangle(2)=1;
end
if (h.rectangle(1)+h.rectangle(3))>ss(1)
    h.rectangle(1)=ss(1)-h.rectangle(3)-1;
end
if (h.rectangle(2)+h.rectangle(4))>ss(2)
    h.rectangle(2)=ss(2)-h.rectangle(4)-1;
end
switch W
case 'work'
    rr=findobj('Tag','r1');set(rr,'Position',h.rectangle)
case 'prev'
    rr=findobj('Tag','r2');set(rr,'Position',h.rectangle)
end

% ----------- Function tom_tilt_lines -----------
function tom_tilt_lines(alig,ref,pkt);
%  TOM_TILT_LINES plots marker points 
%   TOM_TILT_LINES(ALIG,REF,PKT) 
%
%   alig=   Alignment array (from marker file)
%   ref=    Reference point
%   pkt=    Which Point?
txt=1;
all=0;
if ref==pkt
    all=1;
    txt=0;
end
[s1,s2,s3]=size(alig);
if all==0
  plot(0,0,'r+');hold on;zoom on;
  inds=find((alig(2,:,ref) > -1) & (alig(3,:,ref) > -1) & (alig(2,:,pkt) > -1) & (alig(3,:,pkt) > -1));
  alig_x=alig(2,inds,ref)-alig(2,inds,pkt);
  alig_y=alig(3,inds,ref)-alig(3,inds,pkt);
  for kk=1:size(inds,2)
    plot(alig_x(kk),alig_y(kk),'r+');
    if txt==1
      text('Position',[alig_x(kk)+1 alig_y(kk)+1],'String',alig(11,inds(kk),1)); 
    end
  end
grid off;
hold on;
end
if all==1
  clf; plot(0,0,'r+');hold on;
  for ref=1:s3
    for pkt=ref+1:s3
        inds=find((alig(2,:,ref) > -1) & (alig(3,:,ref) > -1) & (alig(2,:,pkt) > -1) & (alig(3,:,pkt) > -1));
        alig_x=alig(2,inds,ref)-alig(2,inds,pkt);
        alig_y=alig(3,inds,ref)-alig(3,inds,pkt);
        plot(alig_x,alig_y,'+');
        text('Position',[alig_x(1) alig_y(1)],'String',ref,'Fontsize',20);   
        text('Position',[alig_x(1)+40 alig_y(1)-40],'String',pkt,'Fontsize',20);   
    end
  end
end

% ----------- Function tom_refresh_ima -----------
function handles=tom_refresh_ima(ImageNr,Where,handles)
% TOM_REFRESH_IMA refresh the working image and the information image. It
%   displays the markers of the current pictures. This function is called 
%   when the user click on the button:
%       -'Load' a project
%       -'Show' information image
%       -'Add' a marker
%       -'Delete'a marker
%       -'Start' an automatic procedure  
am=findobj('Tag','menu');param=get(am,'Userdata');
Nam=get(handles.menu_activemk,'String');
for n=1:2
    switch Where
    case 'work'  %refresh working image    
        if n==1
            if str2num(ImageNr)<=0
                break;
            end
            rect1=handles.rectangle;
            m=findobj('Tag','Figwo');
            set(0,'CurrentFigure',m);
            tom_setmark_imagesc(handles.WorkingImage.Value);
            title('WORKING IMAGE','FontWeight','bold','Units','normalized');
            Xz=double(param.Matrixmark(2,str2num(ImageNr),str2num(Nam)));
            Yz=double(param.Matrixmark(3,str2num(ImageNr),str2num(Nam)));
            showima=str2num(ImageNr);
            if Xz~=-1 & Yz~=-1
                r1=rectangle('Position',[Xz-rect1(3)/2,Yz-rect1(4)/2,rect1(3) rect1(4)],'Tag','r1');
            else
                r1=rectangle('Position',[rect1(1),rect1(2),rect1(3),rect1(4)],'Tag','r1');               
            end
            rect1=get(r1,'Position');                 
            handles.rectangle=rect1;
            set(m,'Name',['...' param.myfilename ImageNr param.myext '  Angle: ' num2str(param.angle(showima))]);
            set(handles.menu_wi,'String',ImageNr);
            set(handles.menu_wi_angle,'String',num2str(param.angle(str2num(ImageNr))));
            message=['...' param.myfilename ImageNr param.myext];
            set(handles.menu_name_ima,'String',message);
            set(handles.menu_showwi,'String','');            
        else
            m=findobj('Tag','Figwo_z');
            set(0,'CurrentFigure',m);
            showima=str2num(ImageNr);
            set(m,'Name',['...' param.myfilename ImageNr param.myext '  Angle: ' num2str(param.angle(showima))]);
            handles.rectangle=round(handles.rectangle);
            handles=tom_limz('work',handles);rect1=handles.rectangle;
            temp=handles.WorkingImage.Value(rect1(1):(rect1(1)+rect1(3)),rect1(2):(rect1(2)+rect1(4)));
            if get(handles.menu_contrast,'Value')%contrast
                if  get(handles.cont_auto,'Value') %automatic
                    av=mean2(temp);standart=std2(temp);                    
                    maxi=max(max(temp));mini=min(min(temp));                    
                    low=av-(3*standart);hi=av+(3*standart);
                    set(handles.cont_mean,'String',av);
                    set(handles.cont_std,'String',standart);
                    set(handles.scale1,'String',low);
                    set(handles.scale2,'String',hi);                    
                else%manual
                    low=str2num(get(handles.scale1,'String'));
                    hi=str2num(get(handles.scale2,'String'));
                end
                tom_imagesc(temp,'noinfo','range',[low hi]);
            elseif get(handles.menu_filter,'Value')%filter
                low=str2num(get(handles.filter_val_low,'String'));
                hi=str2num(get(handles.filter_val_hi,'String'));
                switch handles.FilterName
                    case 'lowpass filter'
                        rect1=round(get(findobj(0,'Tag','r1'),'Position'));
                        temp=handles.WorkingImage.Value(rect1(1):(rect1(1)+rect1(3)),rect1(2):(rect1(2)+rect1(4)));
                        im=tom_bandpass(double(temp),0,low);
                        tom_imagesc(im,'noinfo');
                    case 'highpass filter'
                        rect1=round(get(findobj(0,'Tag','r1'),'Position'));
                        temp=handles.WorkingImage.Value(rect1(1):(rect1(1)+rect1(3)),rect1(2):(rect1(2)+rect1(4)));
                        im=tom_bandpass(double(temp),hi,500);
                        tom_imagesc(im,'noinfo');
                    case 'bandpass filter'
                        rect1=round(get(findobj(0,'Tag','r1'),'Position'));
                        temp=handles.WorkingImage.Value(rect1(1):(rect1(1)+rect1(3)),rect1(2):(rect1(2)+rect1(4)));
                        im=tom_bandpass(double(temp),low,hi);
                        tom_imagesc(im,'noinfo');
                end
            else
                tom_setmark_imagesc(temp)%display the file
            end                        
            title('WORKING IMAGE ZOOMED','FontWeight','bold','Units','normalized');
        end        
        if param.Matrixmark(2:3,:,:)==-1
            %no action
        else
            Currmark=size(param.Matrixmark,3);
            for i=1:Currmark
                if n==1                    
                    Xz=double(param.Matrixmark(2,str2num(ImageNr),i));x=Xz;
                    Yz=double(param.Matrixmark(3,str2num(ImageNr),i));y=Yz;
                else
                    Xz=double(param.Matrixmark(2,str2num(ImageNr),i));
                    Yz=double(param.Matrixmark(3,str2num(ImageNr),i));
                    x=Xz-handles.rectangle(1);
                    y=Yz-handles.rectangle(2);  
                end
                if Xz~=-1 & Yz~=-1
                    tom_drawmark(x,y,num2str(i));
                end    
            end 
        end        
    case 'prev' %refresh information image
        if n==1
            rect2=handles.rectangle;
            m=findobj('Tag','Fig_previ');
            set(0,'CurrentFigure',m);
            if str2num(ImageNr)==0 | str2num(ImageNr)>str2num(get(handles.menu_number_ima,'String'))
                tom_setmark_imagesc(0);
                set(m,'Name','No picture' );
                title('INFORMATION IMAGE','FontWeight','bold','Units','normalized');
                set(handles.menu_infoim,'String','');     
                set(handles.menu_infoim_angle,'String','');
                set(handles.menu_no_ima,'String','');
                m=findobj('Tag','Fig_previ_z');
                set(0,'CurrentFigure',m);
                set(m,'Name','No picture' );
                tom_setmark_imagesc(0);
                title('INFORMATION IMAGE ZOOMED','FontWeight','bold','Units','normalized');
                break;
            else  
                tom_setmark_imagesc(handles.PresentationImage.Value); 
                tit=title('INFORMATION IMAGE','FontWeight','bold','Units','normalized');
                Xz=double(param.Matrixmark(2,str2num(ImageNr),str2num(Nam)));
                Yz=double(param.Matrixmark(3,str2num(ImageNr),str2num(Nam)));
                showima=str2num(ImageNr);
                if Xz~=-1 & Yz~=-1
                    r2=rectangle('Position',[Xz-rect2(3)/2,Yz-rect2(4)/2,rect2(3) rect2(4)],'Tag','r2');
                else
                    r2=rectangle('Position',[rect2(1),rect2(2),rect2(3),rect2(4)],'Tag','r2');
                end
                rect2=get(r2,'Position');
                handles.rectangle=rect2;                
                set(m,'Name',['...' param.myfilename ImageNr param.myext '  Angle: ' num2str(param.angle(showima))]);
                set(handles.menu_infoim,'String',ImageNr);     
                set(handles.menu_infoim_angle,'String',num2str(param.angle(str2num(ImageNr))));
                set(handles.menu_no_ima,'String',ImageNr);
            end
        else      
            showima=str2num(ImageNr);
            m=findobj('Tag','Fig_previ_z');
            set(0,'CurrentFigure',m);
            set(m,'Name',['...' param.myfilename ImageNr param.myext '  Angle: ' num2str(param.angle(showima))]);
            handles.rectangle=round(handles.rectangle);
            handles=tom_limz('prev',handles);rect2=handles.rectangle;
            temp=handles.PresentationImage.Value(rect2(1):(rect2(1)+rect2(3)),rect2(2):(rect2(2)+rect2(4)));
            tom_setmark_imagesc(temp);%display the file
            title('INFORMATION IMAGE ZOOMED','FontWeight','bold','Units','normalized');
        end
        Currmark=size(param.Matrixmark,3);
        for i=1:Currmark
            if n==1
                Xz=double(param.Matrixmark(2,str2num(ImageNr),i));x=Xz;
                Yz=double(param.Matrixmark(3,str2num(ImageNr),i));y=Yz;
            elseif n==2
                Xz=double(param.Matrixmark(2,str2num(ImageNr),i));
                Yz=double(param.Matrixmark(3,str2num(ImageNr),i));
                x=Xz-handles.rectangle(1);
                y=Yz-handles.rectangle(2);
            end
            if Xz~=-1 & Yz~=-1
                tom_drawmark(x,y,num2str(i));
            end            
        end       
    end 
end

% ----------- Function tom_turbo -----------
function handles=tom_turbo(ImageNr,Where,handles)
% TOM_TURBO is the same routine as tom_refresch_ima but 
%   with option turbo  
%   Syntax: tom_turbo(ImageNr,Where, handles)
%       Input:
%           ImageNr: The image number to refresh. Must be a string
%           Where: If 'work' then refresh the working image
%                  If 'prev' then refresh the information image
%           handles: all the handles
%       Output:
%           handles: all the handles 
am=findobj('Tag','menu');param=get(am,'Userdata');
Nam=get(handles.menu_activemk,'String');
switch Where
case 'work'  %refresh working image   
    m=findobj('Tag','Figwo');
    set(0,'CurrentFigure',m);
    clf;
    set(m,'Name','           MODE TURBO' );
    m=findobj('Tag','Figwo_z');
    set(0,'CurrentFigure',m);
    showima=str2num(ImageNr);
    set(m,'Name',['...' param.myfilename ImageNr param.myext '  Angle: ' num2str(param.angle(showima))]);    
    if get(handles.menu_contrast,'Value')%contrast
        if  get(handles.cont_auto,'Value') %automatic
            av=mean2(handles.WorkingImage.Value);standart=std2(handles.WorkingImage.Value);
            maxi=max(max(handles.WorkingImage.Value));mini=min(min(handles.WorkingImage.Value));
            low=av-(3*standart);hi=av+(3*standart);
            set(handles.cont_mean,'String',av);
            set(handles.cont_std,'String',standart);
            set(handles.scale1,'String',low);
            set(handles.scale2,'String',hi);
        else%manual
            low=str2num(get(handles.scale1,'String'));
            hi=str2num(get(handles.scale2,'String'));
        end
        tom_imagesc(handles.WorkingImage.Value,'noinfo','range',[low hi]);
    elseif get(handles.menu_filter,'Value')%filter
        low=str2num(get(handles.filter_val_low,'String'));
        hi=str2num(get(handles.filter_val_hi,'String'));
        switch handles.FilterName
            case 'lowpass filter'
                im=tom_bandpass(double(handles.WorkingImage.Value),0,low);
                tom_imagesc(im,'noinfo');
            case 'highpass filter'
                im=tom_bandpass(double(handles.WorkingImage.Value),hi,500);
                tom_imagesc(im,'noinfo');
            case 'bandpass filter'
                im=tom_bandpass(double(handles.WorkingImage.Value),low,hi);
                tom_imagesc(im,'noinfo');
        end
    else
        tom_setmark_imagesc(handles.WorkingImage.Value)%display the file
    end
    title('WORKING IMAGE ZOOMED','FontWeight','bold','Units','normalized');
    set(handles.menu_wi,'String',ImageNr);
    set(handles.menu_wi_angle,'String',num2str(param.angle(str2num(ImageNr))));
    message=['...' param.myfilename ImageNr param.myext];
    set(handles.menu_name_ima,'String',message);
    set(handles.menu_showwi,'String','');               
    if param.Matrixmark(2:3,:,:)==-1
        %no action
    else
        Currmark=size(param.Matrixmark,3);
        for i=1:Currmark
            Xz=double(param.Matrixmark(2,str2num(ImageNr),i));
            Yz=double(param.Matrixmark(3,str2num(ImageNr),i));
            x=Xz-handles.rectangle(1);
            y=Yz-handles.rectangle(2);  
            if Xz~=-1 & Yz~=-1
                tom_drawmark(x,y,num2str(i));
            end    
        end 
    end        
case 'prev' %refresh information image  
    m=findobj('Tag','Fig_previ');
    set(0,'CurrentFigure',m);
    clf;
    set(m,'Name','           MODE TURBO' );
    if str2num(ImageNr)==0 | str2num(ImageNr)>str2num(get(handles.menu_number_ima,'String'))
        set(handles.menu_infoim,'String','');     
        set(handles.menu_infoim_angle,'String','');
        set(handles.menu_no_ima,'String','');
        m=findobj('Tag','Fig_previ_z');
        set(0,'CurrentFigure',m);
        set(m,'Name','No picture' );
        tom_setmark_imagesc(0);
        return;
    else         
        m=findobj('Tag','Fig_previ_z');
        set(0,'CurrentFigure',m);       
        showima=str2num(ImageNr);
        set(m,'Name',['...' param.myfilename ImageNr param.myext '  Angle: ' num2str(param.angle(showima))]);
        tom_setmark_imagesc(handles.PresentationImage.Value);%tom_setmark_imagesc(handles.PresentationImage.Value,'fixed');%display the file
        title('INFORMATION IMAGE ZOOMED','FontWeight','bold','Units','normalized');
        set(handles.menu_infoim,'String',ImageNr);     
        set(handles.menu_infoim_angle,'String',num2str(param.angle(str2num(ImageNr))));
        set(handles.menu_no_ima,'String',ImageNr);        
    end
    Currmark=size(param.Matrixmark,3);
    for i=1:Currmark
        Xz=double(param.Matrixmark(2,str2num(ImageNr),i));
        Yz=double(param.Matrixmark(3,str2num(ImageNr),i));
        x=Xz-handles.rectangle(1);
        y=Yz-handles.rectangle(2);
        if Xz~=-1 & Yz~=-1
            tom_drawmark(x,y,num2str(i));
        end            
    end
end       

% ----------- Function tom_refresh_mark -----------
function handles=tom_refresh_mark(ImageNr,Data,Where,handles)
% TOM_REFRESH_MARK displays the markers of the current pictures.
%   Syntax: handles=tom_refresh_mark(ImageNr,Data,Where,handles)
%       Input:
%           ImageNr: The image number to refresh. Must be a string
%           Data: Matrixmark
%           Where: If 'work' then refresh the working image
%                  If 'workz' then refresh the working image
%                  If 'prev' then refresh the information image
%                  If 'prevz' then refresh the information image
%           handles: all of the handles
%       Output:
%           handles: all of the handles 
am=findobj('Tag','menu');param=get(am,'Userdata');
Nam=get(handles.menu_activemk,'String');
Currmark=size(Data,3);
switch Where
case 'work'  %refresh working image    
    m=findobj('Tag','Figwo');
    set(0,'CurrentFigure',m); 
    Xz=double(Data(2,str2num(ImageNr),str2num(Nam)));
    Yz=double(Data(3,str2num(ImageNr),str2num(Nam)));
    rect1=handles.rectangle;
    if Xz~=-1 & Yz~=-1
        m=findobj('Tag','r1'); 
        set(m,'Position',[Xz-rect1(3)/2,Yz-rect1(4)/2,rect1(3) rect1(4)]);
        rect1=get(m,'Position');                 
        handles.rectangle=rect1;
    else
        m=findobj('Tag','r1'); 
        r1=rectangle('Position',[rect1(1),rect1(2),rect1(3),rect1(4)],'Tag','r1');               
    end                    
    for i=1:Currmark
        if Xz~=-1 & Yz~=-1
            x=Xz;y=Yz;
            tom_drawmark(x,y,num2str(i));
        end    
    end
case 'workz'
    m=findobj('Tag','Figwo_z');
    set(0,'CurrentFigure',m);
    for i=1:Currmark
        Xz=double(Data(2,str2num(ImageNr),i));
        Yz=double(Data(3,str2num(ImageNr),i));
        x=Xz-handles.rectangle(1);
        y=Yz-handles.rectangle(2);
        if Xz~=-1 & Yz~=-1
            tom_drawmark(x,y,num2str(i));
        end    
    end
case 'prev' %refresh information image
    m=findobj('Tag','Fig_previ');
    set(0,'CurrentFigure',m);
    Xz=double(param.Matrixmark(2,str2num(ImageNr),str2num(Nam)));
    Yz=double(param.Matrixmark(3,str2num(ImageNr),str2num(Nam)));
    rect2=handles.rectangle;
    if Xz~=-1 & Yz~=-1
        m=findobj('Tag','r2'); 
        set(m,'Position',[Xz-rect2(3)/2,Yz-rect2(4)/2,rect2(3) rect2(4)]);
        rect2=get(m,'Position');                 
        handles.rectangle=rect2;
    else
        m=findobj('Tag','r2'); 
        set(m,'Position',[rect2(1),rect2(2),rect2(3),rect2(4)]);               
    end                        
    for i=1:Currmark
        if Xz~=-1 & Yz~=-1
            x=Xz;y=Yz;
            tom_drawmark(x,y,num2str(i));
        end    
    end
case 'prevz'
    m=findobj('Tag','Fig_previ_z');
    set(0,'CurrentFigure',m);
    for i=1:Currmark
        Xz=double(Data(2,str2num(ImageNr),i));
        Yz=double(Data(3,str2num(ImageNr),i));
        x=Xz-handles.rectangle(1);
        y=Yz-handles.rectangle(2);
        if Xz~=-1 & Yz~=-1
            tom_drawmark(x,y,num2str(i));
        end    
    end
end
% ----------- Funcion Prediction -----------
function cp=prediction(p,im,dir)
% Predict where x & y of the square will be
% syntax: cp=prediction(p,im)
% Input
%   p: parameter
%   im: current image
%   dir: direction. 'up' or 'down'
% Output
%   cp: cp(1) is x and cp(2) is y
cm=p.currentmark;ci=im;flag=0;
switch dir
    case 'up'
        %prev mk doesn't exist in prev. & curr. image, or curr. mk doesn't exist in prev.image
        if (p.Matrixmark(2,ci,cm-1)==-1) | (p.Matrixmark(2,ci-1,cm-1)==-1) | (p.Matrixmark(2,ci-1,cm)==-1)
            while p.Matrixmark(2,ci-1,cm)==-1
                ci=ci-1;
            end
            while (p.Matrixmark(2,ci,cm-1)==-1) | (p.Matrixmark(2,ci-1,cm-1)==-1)
                cm=cm-1;
                if cm-1<=0%no vector found between prev & curr image
                    cm=p.currentmark;ni=ci;flag=2;
                    while p.Matrixmark(2,ci,cm-1)==-1%search marker in curr. image
                        cm=cm-1;
                        if cm-1<=0
                            flag=1;break;
                        end
                    end
                    if (cm-1<=0) & (flag==1)
                        break;
                    end;
                    while p.Matrixmark(2,ni-1,cm-1)==-1% search marker in prev. image
                        ni=ni-1;
                        if ni<=str2num(p.image_ref)%no vector found up to ref. image
                            ni=im;flag=1; break;                           
                        end                        
                    end
                    break;
                end
            end
            if flag==1%no vector found so keep the postion of the square
                if get(findobj(0,'Tag','menu_turbo'),'Value')
                    h_menu=findobj(0,'Tag','menu');h=guidata(h_menu);                
                    rect=h.rectangle; 
                else
                    rect=get(findobj(0,'Tag','r1'),'Position');
                end
                cp(1)=rect(1)+rect(3)/2;
                cp(2)=rect(2)+rect(4)/2;
            elseif flag==2%calculate vector on image, not utterly on prev. image
                cp(1)=double(p.Matrixmark(2,ni-1,p.currentmark)) - double(p.Matrixmark(2,ni-1,cm-1)) + double(p.Matrixmark(2,ci,cm-1));
                cp(2)=double(p.Matrixmark(3,ni-1,p.currentmark)) - double(p.Matrixmark(3,ni-1,cm-1)) + double(p.Matrixmark(3,ci,cm-1));
            else %calculate vector on prev. image
                cp(1)=double(p.Matrixmark(2,ci-1,p.currentmark)) - double(p.Matrixmark(2,ci-1,cm-1)) + double(p.Matrixmark(2,ci,cm-1));
                cp(2)=double(p.Matrixmark(3,ci-1,p.currentmark)) - double(p.Matrixmark(3,ci-1,cm-1)) + double(p.Matrixmark(3,ci,cm-1));
            end
        else %prev mk exist in prev. & curr. image
            cp(1)=double(p.Matrixmark(2,im-1,cm)) - double(p.Matrixmark(2,ci-1,cm-1)) + double(p.Matrixmark(2,ci,cm-1));
            cp(2)=double(p.Matrixmark(3,im-1,cm)) - double(p.Matrixmark(3,ci-1,cm-1)) + double(p.Matrixmark(3,ci,cm-1));
        end
    case 'down'
        if (p.Matrixmark(2,ci,cm-1)==-1) | (p.Matrixmark(2,ci+1,cm-1)==-1)%prev mk doesn't exist in prev. & curr. image
            while (p.Matrixmark(2,ci,cm-1)==-1) | (p.Matrixmark(2,ci+1,cm-1)==-1)
                cm=cm-1;
                if cm-1<=0%no vector found between prev & curr image
                    cm=p.currentmark;ni=ci;flag=2;
                    while (p.Matrixmark(2,ci,cm-1)==-1)% 
                        cm=cm-1;
                        if cm-1<=0
                            flag=1;break;
                        end
                    end            
                    if (cm-1<=0) & (flag==1)
                        break;
                    end
                    while p.Matrixmark(2,ni+1,cm-1)==-1
                        ni=ni+1;
                        if ni+1>=str2num(p.image_ref)
                            ni=im;flag=1;break;
                        end
                    end
                    break;
                end
            end
            if flag==1
                if get(findobj(0,'Tag','menu_turbo'),'Value')
                    h_menu=findobj(0,'Tag','menu');h=guidata(h_menu);                
                    rect=h.rectangle; 
                else
                    rect=get(findobj(0,'Tag','r1'),'Position');
                end
                cp(1)=rect(1)+rect(3)/2;
                cp(2)=rect(2)+rect(4)/2;
            elseif flag==2
                cp(1)=double(p.Matrixmark(2,ni+1,p.currentmark)) - double(p.Matrixmark(2,ni+1,cm-1)) + double(p.Matrixmark(2,ci,cm-1));
                cp(2)=double(p.Matrixmark(3,ni+1,p.currentmark)) - double(p.Matrixmark(3,ni+1,cm-1)) + double(p.Matrixmark(3,ci,cm-1));
            else
                cp(1)=double(p.Matrixmark(2,im+1,p.currentmark)) - double(p.Matrixmark(2,ci+1,cm-1)) + double(p.Matrixmark(2,ci,cm-1));
                cp(2)=double(p.Matrixmark(3,im+1,p.currentmark)) - double(p.Matrixmark(3,ci+1,cm-1)) + double(p.Matrixmark(3,ci,cm-1));
            end
        else
            cp(1)=double(p.Matrixmark(2,im+1,p.currentmark)) - double(p.Matrixmark(2,ci+1,cm-1)) + double(p.Matrixmark(2,ci,cm-1));
            cp(2)=double(p.Matrixmark(3,im+1,p.currentmark)) - double(p.Matrixmark(3,ci+1,cm-1)) + double(p.Matrixmark(3,ci,cm-1));
        end
end
cp=round(cp);
% ----------- Functionfillbymean -----------
function [fill,rect]=fillbymean(name_ima,cx,cy,w,h)
% FILLBYMEAN fill the matrix by the mean of the image when 
% tom_emread(...,'subregion'...) is used and the area subregion
% is outside of the image. 
% Syntax: [fill,rec]=fillbymean(name_ima,cx,cy,w,h)
%     Input:
%         name_ima: name of the image (with path)
%         cx: coordinate x of the region to read
%         cy: coordinate y of the region to read
%         w: width of the region to read
%         h: height of the region to read
%     Output:
%         fill: Struct of the image corrected by the cont_mean
%         rect: Coordinate of the rectangle
w_org=w;h_org=h;XX=1;yy=1;flag=0;
rh=tom_reademheader(name_ima);
if (cy<1) & ((cx+w)>rh.Header.Size(1))%case of coord y < size(im) AND rectangle width > size(im). Must be before tom_emread
    XX=(cx+w)-rh.Header.Size(1);w=w-XX;
    YY=1-cy;h=h-YY;
elseif cy<1 %case of coord y  < size(im). Must be before tom_emread
    YY=1-cy;h=h-YY;
end
fill=tom_emread(name_ima,'subregion',[cx cy 1],[w-1 h-1 0]);%fill.Value=double(fill.Value);  
if ((cy+h)>rh.Header.Size(2)) & (cx<1)%case of rectangle height > size(im) AND x coord < size(im)  
    tt=h_org -((cy+h)-rh.Header.Size(2));
    av=mean(mean(fill.Value(-cx:w_org,1:tt)));% mean(mean(..)) to have 1 value
    fill.Value(1:(1-cx),:)=av;%fill the empty area by the mean
    fill.Value(:,tt:h_org)=av;
    rect=[cx cy w_org h_org];
    return;
elseif ((cy+h)>rh.Header.Size(2)) & ((cx+w_org)>rh.Header.Size(1))%case of rectangle height > size(im) AND rectangle width > size(im)
    ttx=w_org -((cx+w_org)-rh.Header.Size(1));
    tty=h_org -((cy+h)-rh.Header.Size(2));
    av=mean(mean(fill.Value(1:ttx,1:tty)));
    fill.Value(ttx:w_org,:)=av;
    fill.Value(:,tty:h_org)=av;
    rect=[cx cy w_org h_org];  
    return;
elseif (cx<1) & (cy<1)%case of coord x < size(im) AND coord y < size(im)
    av=mean(mean(fill.Value((1-cx):w_org,:)));
    ma=zeros(round(w_org),round(h_org));ma(:)=av;    
    fill.Value=tom_paste(ma,fill.Value,[1,(1-cy)]);%fill the empty area by the mean
    fill.Value(1:(1-cx),:)=av;
    rect=[cx cy w_org h_org];  
    return;       
elseif cx<0 %case of coord x < size(im)
    av=mean(mean(fill.Value(-cx:w_org,:)));
    fill.Value(1:(1-cx),:)=av; 
elseif cy<1 %case of coor y < size(im). Also for coor y < size(im) AND rectangle width > size(im)
    av=mean(mean(fill.Value));
    ma=zeros(round(w_org),round(h_org));ma(:)=av;
    fill.Value=tom_paste(ma,fill.Value,[1,(1-cy)]);
elseif (cx+w_org)>rh.Header.Size(1)%case of rectangle width > size(im) 
    tt=w_org -((cx+w_org)-rh.Header.Size(1));
    av=mean(mean(fill.Value(1:tt,:)));
    ma=zeros(round(w_org),round(h_org));ma(:)=av;
    fill.Value=tom_paste(fill.Value,ma,[tt,1]);
elseif (cy+h)>rh.Header.Size(2); %case of rectangle height > size(im)
    tt=h_org -((cy+h)-rh.Header.Size(2));
    av=mean(mean(fill.Value(:,1:tt)));
    ma=zeros(round(w_org),round(h_org));ma(:)=av;
    fill.Value=tom_paste(fill.Value,ma,[1,tt]);
end
rect=[cx cy w_org h_org];

% ----------- Function SetZoomPointer -----------
% change the cursor in WI.
function SetZoomPointer
P = ones(16);P(:,:)=NaN;P(:,:)=2;
P(1,:) = 2; P(16,:) = 2;
P(:,1) = 2; P(:,16) = 2;
P(4:13,4:13) = NaN; % Create a transparent region in the center
P(1:4,8:9) = 1; P(13:16,8:9) = 1;
P(8:9,1:4) = 1; P(8:9,13:16) = 1;
set(findobj(0,'Tag','Figwo'),'Pointer','custom',...
'PointerShapeCData',P,'PointerShapeHotSpot',[9 9]);%pointer shape WI

% ----------- Function SetClickPointer -----------
% change the zoom in WIZ. 
function SetClickPointer
P = ones(16);P(:,:)=NaN;
P(8,1:6)=2; P(8,10:16)=2;P(1:6,8)=2;P(10:16,8)=2;
set(findobj(0,'Tag','Figwo_z'),'Pointer','custom',...
'PointerShapeCData',P,'PointerShapeHotSpot',[8 8]);

% ----------- Function InitMenu -----------
% Iniatialze the button, text and all the object in the menu
function InitMenu(handles);
m=findobj('Tag','Figwo');delete(m);
m=findobj('Tag','Figwo_z');delete(m);
m=findobj('Tag','Fig_previ');delete(m);
m=findobj('Tag','Fig_previ_z');delete(m);
set(handles.menu_infoim,'String','');% text info. image
set(handles.menu_infoim_angle,'String','');%text angle info. image
set(handles.menu_showii,'String','');%editbox show info. image
set(handles.menu_showwi,'String','');%editbox show working image
set(handles.menu_contrast,'Value',0,'Enable','on');%checkbutton apply contrast
set(handles.menu_filter,'Value',0,'Enable','on');%checkbutton apply filter
set(handles.menu_activemk,'String','1');%text active marker
set(handles.menu_setmk,'String','');%editbox active existing marker
set(handles.menu_radioup,'Value',1,'Enable','on');%radiobutton up
set(handles.menu_radiodown,'Value',0,'Enable','on');%radiobutton down
set(handles.menu_turbo,'Value',0,'Enable','on');%checkbutton turbo
set(handles.menu_no_mark,'String','');%edit box marker number in movie panel
set(handles.menu_no_ima,'String','');%text info. image n??? in movie panel
set(handles.menu_name_mf,'String','');%text marker file
set(handles.menu_definemarker_tilt,'String','');%editbox markerpoint alignment
set(handles.menu_tilt1,'String','');%editbox tilt line 1
set(handles.menu_tilt2,'String','');%editbox tilt line 2
set(handles.menu_error,'String','');%text error
set(handles.menu_turbo,'Enable','on');%checkbox turbo
set(handles.menu_cross,'Enable','on');%checkbox cross correlation
EnableAllButton('on');%enable all button


% ----------- Function EnableApplyContrast -----------
% set the property 'Enable' (on or off) of the panel Contrast 
function EnableApplyContrast(Status);
h_menu=findobj(0,'Tag','menu');
handles=guidata(h_menu);
set(handles.cont_auto,'Enable',Status);%radiobutton automatic
set(handles.cont_manual,'Enable',Status);%radiobutton manual
set(handles.text_mean,'Enable',Status);%text mean
set(handles.text_std,'Enable',Status);%text std
set(handles.text_min,'Enable',Status);%text min
set(handles.text_max,'Enable',Status);%text max
set(handles.cont_mean,'Enable',Status);%value mean
set(handles.cont_std,'Enable',Status);%value std
set(handles.scale1,'Enable',Status);%value min
set(handles.scale2,'Enable',Status);%value max
set(handles.menu_cont_doit,'Enable',Status);%button do it

% ----------- Function EnableApplyFilter -----------
% set the property 'Enable' (on or off) of the panel Filter (pix)
function EnableApplyFilter(Status);
h_menu=findobj(0,'Tag','menu');
handles=guidata(h_menu);
set(handles.filtertype,'Enable',Status);%popup menu filter
set(handles.filter_val_low,'Enable',Status);%editbox filter low
set(handles.filter_val_hi,'Enable',Status);%editbox filter hi
set(handles.filter_low,'Enable',Status);%text low
set(handles.filter_hi,'Enable',Status);%text hi
set(handles.menu_filt_doit,'Enable',Status);%button do it

% ----------- Function EnableCrossCorrelation -----------
%set the property 'Enable' and 'Visible' (on or off) of the panel 
%Cross Correlation parameter 
function EnableCrossCorrelation(Status)
h_menu=findobj(0,'Tag','menu');
handles=guidata(h_menu);
set(handles.menu_radioup,'Enable',Status);%radiobutton up
set(handles.menu_radiodown,'Enable',Status);%radiobutton down
set(handles.menu_turbo,'Enable',Status);%radiobutton both
VisibleFrameMovie(Status);%uipanel movie
VisibleFrameMarkerFile(Status);%uipanel marker file
EnableAllButton(Status);
switch Status
    case 'off'
        VisibleFrameCC('on');%uipanel cross correlation
    case 'on'
        VisibleFrameCC('off');
end
% ----------- Function EnableAllButton -----------
%set the properties 'Enable' (on or off) of all Button (except Quit button)
function EnableAllButton(Status)
h_menu=findobj(0,'Tag','menu');
handles=guidata(h_menu);
set(handles.new_conf,'Enable',Status);%button New project
set(handles.load_conf,'Enable',Status);%button Load project
set(handles.save_conf,'Enable',Status);%button Save project
set(handles.menu_prev,'Enable',Status);%button Previous image
set(handles.menu_next,'Enable',Status);%button Next image
set(handles.menu_zoom,'Enable',Status);%button Zoom
set(handles.menu_showworkima,'Enable',Status);%button Show Working image
set(handles.menu_showqwim,'Enable',Status);%button Show Information image
set(handles.menu_cont_doit,'Enable',Status);%button constrast do it
set(handles.menu_filt_doit,'Enable',Status);%button filter do it
set(handles.menu_active,'Enable',Status);%button Active marker
set(handles.menu_rename,'Enable',Status);%button New marker
set(handles.marker_add,'Enable',Status);%button Add marker
set(handles.marker_del,'Enable',Status);%button Delete marker
set(handles.menu_autostart,'Enable',Status);%button Start automatic procedure
set(handles.menu_autostop,'Enable',Status);%button Stop automatic procedure
set(handles.menu_showmark,'Enable',Status);%button Show marker
set(handles.menu_movie,'Enable',Status);%button Movie
set(handles.menu_loadmarkfile,'Enable',Status);%button Load MF
set(handles.menu_savemarkfile,'Enable',Status);%button Save MF
set(handles.menu_view,'Enable',Status);%button View MF
set(handles.menu_tilt,'Enable',Status);%button Tilt Lines
set(handles.menu_align,'Enable',Status);%button Alignment
% ----------- Function VisibleFrameMovie -----------
% set the properties 'Enable' (on or off) of the uipanel movie
function VisibleFrameMovie(Status)
h_menu=findobj(0,'Tag','menu');
handles=guidata(h_menu);
set(handles.frame_movie,'Visible',Status);%frame
set(handles.mov_title,'Visible',Status);%text
set(handles.mov1,'Visible',Status);%text
set(handles.mov2,'Visible',Status);%text
set(handles.mov3,'Visible',Status);%text
set(handles.menu_no_mark,'Visible',Status);%edit box
set(handles.menu_no_ima,'Visible',Status);%text
set(handles.menu_movie,'Visible',Status);%button movie
set(handles.menu_showmark,'Visible',Status);%button show
% ----------- VisibleFrameMarkerFile -----------
% set the properties 'Enable' (on or off) of the uipanel marker file
function VisibleFrameMarkerFile(Status)
h_menu=findobj(0,'Tag','menu');
handles=guidata(h_menu);
set(handles.frame_markerfile,'Visible',Status);%frame
set(handles.mf_title,'Visible',Status);%text
set(handles.menu_name_mf,'Visible',Status);%text
set(handles.menu_loadmarkfile,'Visible',Status);%button load MF
set(handles.menu_savemarkfile,'Visible',Status);%button save MF
set(handles.menu_view,'Visible',Status);%button view MF
% ----------- VisibleFrameMarkerFile -----------
% set the properties 'Enable' (on or off) of the uipanel cross correlation
function VisibleFrameCC(Status)
h_menu=findobj(0,'Tag','menu');
handles=guidata(h_menu);
set(handles.frame_cc,'Visible',Status);%frame CC parameter
set(handles.frame_dir,'Visible',Status);%frame dir
set(handles.frame_filter,'Visible',Status);%frame filter
set(handles.frame_size,'Visible',Status);%frame size of cc
set(handles.paramcc1,'Visible',Status);%text limit ...
set(handles.paramcc2,'Visible',Status);%text Gold ...
set(handles.paramcc3,'Visible',Status);%text Error ...
set(handles.cc_title,'Visible',Status);%text
set(handles.dir_title,'Visible',Status);%text
set(handles.filter_title,'Visible',Status);%text
set(handles.paramcc_title,'Visible',Status);%text
set(handles.filterlow,'Visible',Status);%text low
set(handles.filterhi,'Visible',Status);%text hi
set(handles.sccx,'Visible',Status);%text x
set(handles.sccy,'Visible',Status);%text y
set(handles.dir_up,'Visible',Status);%radio up
set(handles.dir_down,'Visible',Status);%radio down
set(handles.dir_both,'Visible',Status);%radio both
set(handles.menu_filterlow,'Visible',Status);%editbox filter low
set(handles.menu_filterhi,'Visible',Status);%editbox filter hi
set(handles.menu_crx,'Visible',Status);%editbox size CC x
set(handles.menu_cry,'Visible',Status);%editbox size CC y
set(handles.menu_limccf,'Visible',Status);%editbox limit CCF
set(handles.menu_gmp,'Visible',Status);%editbox gold marker radius
set(handles.menu_threshold,'Visible',Status);%editbox error
set(handles.menu_cc_doit,'Visible',Status);%button do it

% ----------- Function MouseDown_WI -----------
% Called by ButtonDownFcn in WI 
% Active when click on WI. 
% Change square position 
function MouseDown_WI(hObject,Status)
h_menu=findobj(0,'Tag','menu');
handles=guidata(h_menu);
am=findobj('Tag','menu');param=get(am,'Userdata');
showima=get(handles.menu_wi,'String');ImageNr=str2num(showima);
point1 = get(gca,'CurrentPoint');    % button down detected
pt = point1(1,1:2);
rect1=handles.rectangle;
Xz= pt(1)-(rect1(3)/2);Yz=pt(2)-(rect1(4)/2);
m=findobj('Tag','r1');set(m,'Position',[Xz Yz rect1(3) rect1(4)]);
handles.rectangle=get(m,'Position');
handles=tom_limz('work',handles);rect1=handles.rectangle;
m=findobj('Tag','Figwo_z');set(0,'CurrentFigure',m);
handles=tom_limz('work',handles);
temp=handles.WorkingImage.Value(round(rect1(1)):round((rect1(1)+rect1(3))),round(rect1(2)):round((rect1(2)+rect1(4))));
if get(handles.menu_contrast,'Value')%contrast
    if  get(handles.cont_auto,'Value') %automatic
        av=mean2(temp);standart=std2(temp);
        maxi=max(max(temp));mini=min(min(temp));
        low=av-(3*standart);hi=av+(3*standart);
        set(handles.cont_mean,'String',av);
        set(handles.cont_std,'String',standart);
        set(handles.scale1,'String',low);
        set(handles.scale2,'String',hi);
    else%manual
        low=str2num(get(handles.scale1,'String'));
        hi=str2num(get(handles.scale2,'String'));
    end
    tom_imagesc(temp,'noinfo','range',[low hi]);
elseif get(handles.menu_filter,'Value')%filter
    low=str2num(get(handles.filter_val_low,'String'));
    hi=str2num(get(handles.filter_val_hi,'String'));
    switch handles.FilterName
        case 'lowpass filter'
            rect1=round(get(findobj(0,'Tag','r1'),'Position'));
            temp=handles.WorkingImage.Value(rect1(1):(rect1(1)+rect1(3)),rect1(2):(rect1(2)+rect1(4)));
            im=tom_bandpass(double(temp),0,low);
            tom_imagesc(im,'noinfo');
        case 'highpass filter'
            rect1=round(get(findobj(0,'Tag','r1'),'Position'));
            temp=handles.WorkingImage.Value(rect1(1):(rect1(1)+rect1(3)),rect1(2):(rect1(2)+rect1(4)));
            im=tom_bandpass(double(temp),hi,500);
            tom_imagesc(im,'noinfo');
        case 'bandpass filter'
            rect1=round(get(findobj(0,'Tag','r1'),'Position'));
            temp=handles.WorkingImage.Value(rect1(1):(rect1(1)+rect1(3)),rect1(2):(rect1(2)+rect1(4)));
            im=tom_bandpass(double(temp),low,hi);
            tom_imagesc(im,'noinfo');
    end
else
    tom_setmark_imagesc(temp)%display the file
end
title('WORKING IMAGE ZOOMED','FontWeight','bold','Units','normalized');
switch Status %when refresh, 'ButtonDownFcn' is reset. So restore 'ButtonDownFcn'
    case 'S' %for Single mark, when button Add is pressed
        %handle_image=findall(findobj(0,'Tag','Figwo_z'),'Type','image');
        %set(handle_image,'ButtonDownFcn','tom_setmark(''MouseDown_AddWIZ'',gcbo)');
        %handle_image=findall(findobj(0,'Tag','Figwo'),'Type','image');
        %set(handle_image,'ButtonDownFcn','tom_setmark(''MouseDown_WI'',gcbo,''S'')');
        SetButtonDownFcn_wiz('tom_setmark(''MouseDown_AddWIZ'',gcbo)');
        SetButtonDownFcn_wi('tom_setmark(''MouseDown_WI'',gcbo,''S'')');
    case 'A' %for Automatic, when button Start is pressed
        handle_image=findall(findobj(0,'Tag','Figwo_z'),'Type','image');
        set(handle_image,'ButtonDownFcn','tom_setmark(''MouseDown_WIZ'',gcbo)');
        handle_image=findall(findobj(0,'Tag','Figwo'),'Type','image');
        set(handle_image,'ButtonDownFcn','tom_setmark(''MouseDown_WI'',gcbo,''A'')');
end
guidata(h_menu,handles);

% ----------- Function MouseDown_WIZ -----------
% Called by ButtonDownFcn in WI zoomed 
% This function is used in the automatic procedure, to put markers.
% It is called just after a click on the Working Image Zoomed
function MouseDown_WIZ(hObject)
h_menu=findobj(0,'Tag','menu');handles=guidata(h_menu);
am=findobj('Tag','menu');param=get(am,'Userdata');
showima=get(handles.menu_wi,'String');ImageNr=str2num(showima);
if get(handles.menu_radioup,'Value')==1 %up
    point1 = get(gca,'CurrentPoint');    % button down detected
    pt = point1(1,1:2); pt=round(pt);
    if get(handles.menu_turbo,'Value')==0 %no mode turbo
        m=findobj('Tag','r1');handles.rectangle=get(m,'Position');
    end
    i=str2num(get(handles.menu_activemk,'String'));
    param.Matrixmark(2,ImageNr,i)=pt(1)+handles.rectangle(1); % in order to conform with the EM-Programm the x and y-Coordinate must be flipped
    param.Matrixmark(3,ImageNr,i)=pt(2)+handles.rectangle(2); % this is done in tom_3d, change by SN, 111002 -> FF to check !!!
    param.currentmark=i;
    set(am,'UserData',param);%write param.Matrixmark
    temp=param.Matrixmark;
    temp=tom_emheader(temp);
    tom_emwrite(param.myfilemarker_default,temp);%write in markerfile by default
    if ImageNr+1>str2num(param.mylastnb)
        handle_image=findall(findobj(0,'Tag','Figwo_z'),'Type','image');
        set(handle_image,'ButtonDownFcn','');
        handle_image=findall(findobj(0,'Tag','Figwo'),'Type','image');
        set(handle_image,'ButtonDownFcn','');
        set(findobj('Tag','figw_z_boxinfo2'),'String','STOP');
        set(findobj('Tag','Figwo'),'Pointer','arrow');
        set(findobj('Tag','Figwo_z'),'Pointer','arrow');        
        %%%%%%%%%%  working image  %%%%%%%%%%%%%
        if get(handles.menu_turbo,'Value')==0 %turbo is off
            handles=tom_refresh_ima(num2str(ImageNr),'work',handles);
        else %turbo is on
            set(handles.menu_turbo,'Value',0);
            handles=update(param,handles);
            guidata(h_menu,handles);
        end
        EnableAllButton('on');
        menu_contrast_Callback(handles.menu_contrast,'',handles);
        menu_filter_Callback(handles.menu_filter, '', handles);
        return;
    end
    %display the next image
    if param.currentmark>1%More than one mark
        pt=prediction(param,ImageNr+1,'up');
    else%less than one mark
        pt(1)=param.Matrixmark(2,ImageNr,param.currentmark);% extract x from Matrimark.So no prediction
        pt(2)=param.Matrixmark(3,ImageNr,param.currentmark);% extract y from Matrimark.So no prediction
    end
    Xz= pt(1)-(handles.rectangle(3)/2);Yz=pt(2)-(handles.rectangle(4)/2);%pt is the coordinate of the prediction, Xz Yz coordinates of rectangle with pt in the center
    aa=[param.mypathname param.myfilename num2str(ImageNr+1) param.myext];%read next picture
    handles.PresentationImage=handles.WorkingImage; %direction up
    if get(handles.menu_turbo,'Value')==0 %turbo is off
        %%%%%%%%%%  working image  %%%%%%%%%%%%%
        handles.rectangle=[Xz Yz handles.rectangle(3) handles.rectangle(4)];
        handles.WorkingImage=tom_emreadc(aa);
        handles=tom_refresh_ima(num2str(ImageNr+1),'work',handles);
        %%%%%%%%%%  information image  %%%%%%%%%%%%%        
        handles=tom_refresh_ima(num2str(ImageNr),'prev',handles);        
    else%turbo is on
        %%%%%%%%%%  information image  %%%%%%%%%%%%%
        m=findobj('Tag','Fig_previ');set(0,'CurrentFigure',m);%Let II & WI in this order otherwise, there is a shift
        tit=get(m,'Name');                                    %on II due to the rectangle
        handles=tom_turbo(num2str(ImageNr),'prev',handles);
        %%%%%%%%%%  working image  %%%%%%%%%%%%%
        handles.rectangle=[Xz Yz handles.rectangle(3) handles.rectangle(4)];
        r1=round(handles.rectangle(1));r2=round(handles.rectangle(2));
        w=round(handles.rectangle(3));h=round(handles.rectangle(4));
        [handles.WorkingImage,handles.rectangle]=fillbymean(aa,r1,r2,w,h);
        handles=tom_turbo(num2str(ImageNr+1),'work',handles);
    end
    %Set ButtonDownFcn after displaying image
    %handle_image=findall(findobj(0,'Tag','Figwo_z'),'Type','image');
    %set(handle_image,'ButtonDownFcn','tom_setmark(''MouseDown_WIZ'',gcbo)');
    %handle_image=findall(findobj(0,'Tag','Figwo'),'Type','image');
    %set(handle_image,'ButtonDownFcn','tom_setmark(''MouseDown_WI'',gcbo,''A'')');
    SetButtonDownFcn_wiz('tom_setmark(''MouseDown_WIZ'',gcbo)');
    SetButtonDownFcn_wi('tom_setmark(''MouseDown_WI'',gcbo,''A'')');

    guidata(h_menu,handles);
else%down
    point1 = get(gca,'CurrentPoint');    % button down detected
    pt = point1(1,1:2); pt=round(pt);
    if get(handles.menu_turbo,'Value')==0 %no mode turbo
        m=findobj('Tag','r1');handles.rectangle=get(m,'Position');
    end
    i=str2num(get(handles.menu_activemk,'String'));
    param.Matrixmark(2,ImageNr,i)=pt(1)+handles.rectangle(1); % in order to conform with the EM-Programm the x and y-Coordinate must be flipped
    param.Matrixmark(3,ImageNr,i)=pt(2)+handles.rectangle(2); % this is done in tom_3d, change by SN, 111002 -> FF to check !!!
    param.currentmark=i;
    set(am,'UserData',param);%write param.Matrixmark
    temp=param.Matrixmark;
    temp=tom_emheader(temp);
    tom_emwrite(param.myfilemarker_default,temp);%write in markerfile by default
    if ImageNr-1<str2num(param.myfirstnb)
        handle_image=findall(findobj(0,'Tag','Figwo_z'),'Type','image');
        set(handle_image,'ButtonDownFcn','');
        handle_image=findall(findobj(0,'Tag','Figwo'),'Type','image');
        set(handle_image,'ButtonDownFcn','');
        set(findobj('Tag','figw_z_boxinfo2'),'String','STOP');
        set(findobj('Tag','Figwo'),'Pointer','arrow');
        set(findobj('Tag','Figwo_z'),'Pointer','arrow');
        %%%%%%%%%%  working image  %%%%%%%%%%%%%
        if get(handles.menu_turbo,'Value')==0 %turbo is off
            handles=tom_refresh_ima(num2str(ImageNr),'work',handles);
        else %turbo is on
            set(handles.menu_turbo,'Value',0);
            handles=update(param,handles);
            guidata(h_menu,handles);
        end
        EnableAllButton('on');
        menu_contrast_Callback(handles.menu_contrast,'',handles);
        menu_filter_Callback(handles.menu_filter, '', handles);                
        return;
    end
   %display the next image
    if param.currentmark>1%More than one mark
        pt=prediction(param,ImageNr-1,'down');
    else%less than one mark
        pt(1)=param.Matrixmark(2,ImageNr,param.currentmark);% extract x from Matrimark.So no prediction
        pt(2)=param.Matrixmark(3,ImageNr,param.currentmark);% extract y from Matrimark.So no prediction
    end
    Xz= pt(1)-(handles.rectangle(3)/2);Yz=pt(2)-(handles.rectangle(4)/2);%pt is the coordinate of the prediction, Xz Yz coordinates of rectangle with pt in the center
    aa=[param.mypathname param.myfilename num2str(ImageNr-1) param.myext];%read next picture
    handles.PresentationImage=handles.WorkingImage; %direction up
    if get(handles.menu_turbo,'Value')==0 %turbo is off
        %%%%%%%%%%  working image  %%%%%%%%%%%%%
        handles.rectangle=[Xz Yz handles.rectangle(3) handles.rectangle(4)];
        handles.WorkingImage=tom_emreadc(aa);
        handles=tom_refresh_ima(num2str(ImageNr-1),'work',handles);
        %%%%%%%%%%  information image  %%%%%%%%%%%%%
        handles=tom_refresh_ima(num2str(ImageNr),'prev',handles);
    else%turbo is on
        %%%%%%%%%%  information image  %%%%%%%%%%%%%
        m=findobj('Tag','Fig_previ');set(0,'CurrentFigure',m);%Let II & WI in this order otherwise, there is a shift
        tit=get(m,'Name');                                    %on II due to the rectangle
        handles=tom_turbo(num2str(ImageNr),'prev',handles);
        %%%%%%%%%%  working image  %%%%%%%%%%%%%
        handles.rectangle=[Xz Yz handles.rectangle(3) handles.rectangle(4)];
        r1=round(handles.rectangle(1));r2=round(handles.rectangle(2));
        w=round(handles.rectangle(3));h=round(handles.rectangle(4));
        [handles.WorkingImage,handles.rectangle]=fillbymean(aa,r1,r2,w,h);
        handles=tom_turbo(num2str(ImageNr-1),'work',handles);
    end
    %Set ButtonDownFcn after displaying image
    %handle_image=findall(findobj(0,'Tag','Figwo_z'),'Type','image');
    %set(handle_image,'ButtonDownFcn','tom_setmark(''MouseDown_WIZ'',gcbo)');
    %handle_image=findall(findobj(0,'Tag','Figwo'),'Type','image');
    %set(handle_image,'ButtonDownFcn','tom_setmark(''MouseDown_WI'',gcbo,''A'')');
    SetButtonDownFcn_wiz('tom_setmark(''MouseDown_WIZ'',gcbo)');
    SetButtonDownFcn_wi('tom_setmark(''MouseDown_WI'',gcbo,''A'')');    
    guidata(h_menu,handles);
end

% ----------- Function MouseDown_AddWIZ -----------
% Called by ButtonDownFcn in WI zoomed 
% Active when click on WIZ 
% Add a marker 
function MouseDown_AddWIZ(hObject)
h_menu=findobj(0,'Tag','menu');handles=guidata(h_menu);
am=findobj('Tag','menu');param=get(am,'Userdata');
%%%%%%%%%%  working image  %%%%%%%%%%%%%
point1 = get(gca,'CurrentPoint');    % button down detected
pt = point1(1,1:2);pt=round(pt);
x=pt(1);y=pt(2);
m=findobj('Tag','r1');handles.rectangle=get(m,'Position');
showima=get(handles.menu_wi,'String');ImageNr=str2num(showima);
Currmark=str2num(get(handles.menu_activemk,'String'));
if param.Matrixmark(2,ImageNr,Currmark) & param.Matrixmark(2,ImageNr,Currmark) ~=-1
    message=['The mark ' num2str(Currmark) ' in the image ' showima ' already exist. Would you like to move it?'];
    Question=questdlg(message,'Delete marker','Yes','No','Yes');
    if strcmp(Question,'Yes')
        %tom_drawmark(x,y,num2str(Currmark));
    elseif strcmp(Question,'No')
        EnableAllButton('on');
        menu_contrast_Callback(handles.menu_contrast,'',handles);
        menu_filter_Callback(handles.menu_filter, '', handles);        
        return;
    end
end
param.Matrixmark(2,ImageNr,Currmark)=x+handles.rectangle(1); % in order to conform with the EM-Programm the x and y-Coordinate must be flipped
param.Matrixmark(3,ImageNr,Currmark)=y+handles.rectangle(2); % this is done in tom_3d, change by SN, 111002 -> FF to check !!!
param.currentmark=Currmark;
param.stop_auto_proc=0; %No action when param.stop_auto_proc=0;
set(am,'UserData',param);
handles=tom_refresh_ima(showima,'work',handles);
%%%%%%%%%%  information image   %%%%%%%%%%%%%
m=findobj('Tag','Fig_previ');asd=get(m,'Visible');
switch asd
    case 'off' %INFORMATION IMAGE  is NOT visible
        if (ImageNr-1)<str2num(param.myfirstnb)
            handles.PresentationImage=0;
        else
            bb=[param.mypathname param.myfilename num2str(ImageNr-1) param.myext];
            handles.PresentationImage=tom_emreadc(bb);
        end
        set(m,'Visible','on');
        m=findobj('Tag','Fig_previ_z');
        set(m,'Visible','on');
        handles=tom_refresh_ima(num2str(ImageNr-1),'prev',handles);
end
param.stop_auto_proc=1; %Stop the automatic procedure
set(am,'UserData',param);
temp=param.Matrixmark;
temp=tom_emheader(temp);
tom_emwrite(param.myfilemarker_default,temp);
EnableAllButton('on');
menu_contrast_Callback(handles.menu_contrast,'',handles);
menu_filter_Callback(handles.menu_filter, '', handles);
set(findobj('Tag','figw_z_boxinfo2'),'String','STOP');
set(findobj('Tag','Figwo'),'Pointer','arrow');
set(findobj('Tag','Figwo_z'),'Pointer','arrow');
SetButtonDownFcn_wi(''); %clear ButtonDownFcn
SetButtonDownFcn_wiz('');
guidata(h_menu,handles);
% ----------- Function MouseDown_zoom -----------
% Called by ButtonDownFcn in WI zoomed 
% change the square position 
function MouseDown_zoom(hObject,Status)
h_menu=findobj(0,'Tag','menu');handles=guidata(h_menu);
am=findobj('Tag','menu');param=get(am,'Userdata');
mm=findobj(0,'Tag','r1');
point1 = get(gca,'CurrentPoint');
pt = point1(1,1:2); pt=round(pt);% extract x and y
rect1=handles.rectangle;
Xz= pt(1)-(rect1(3)/2);
Yz=pt(2)-(rect1(4)/2);
set(mm,'Position',[Xz Yz rect1(3) rect1(4)]);
rect1=get(mm,'Position');
handles.rectangle=rect1;
m=findobj('Tag','Figwo_z');set(0,'CurrentFigure',m);
a=handles.WorkingImage;
handles=tom_limz('work',handles);rect1=handles.rectangle;
temp=a.Value(round(rect1(1)):(round(rect1(1)+rect1(3))),round(rect1(2)):round((rect1(2)+rect1(4))));
if get(handles.menu_contrast,'Value')%contrast
    if  get(handles.cont_auto,'Value') %automatic
        av=mean2(temp);standart=std2(temp);
        maxi=max(max(temp));mini=min(min(temp));
        low=av-(3*standart);hi=av+(3*standart);
        set(handles.cont_mean,'String',av);
        set(handles.cont_std,'String',standart);
        set(handles.scale1,'String',low);
        set(handles.scale2,'String',hi);
    else%manual
        low=str2num(get(handles.scale1,'String'));
        hi=str2num(get(handles.scale2,'String'));
    end
    tom_imagesc(temp,'noinfo','range',[low hi]);
elseif get(handles.menu_filter,'Value')%filter
    low=str2num(get(handles.filter_val_low,'String'));
    hi=str2num(get(handles.filter_val_hi,'String'));
    switch handles.FilterName
        case 'lowpass filter'
            rect1=round(get(findobj(0,'Tag','r1'),'Position'));
            temp=handles.WorkingImage.Value(rect1(1):(rect1(1)+rect1(3)),rect1(2):(rect1(2)+rect1(4)));
            im=tom_bandpass(double(temp),0,low);
            tom_imagesc(im,'noinfo');
        case 'highpass filter'
            rect1=round(get(findobj(0,'Tag','r1'),'Position'));
            temp=handles.WorkingImage.Value(rect1(1):(rect1(1)+rect1(3)),rect1(2):(rect1(2)+rect1(4)));
            im=tom_bandpass(double(temp),hi,500);
            tom_imagesc(im,'noinfo');
        case 'bandpass filter'
            rect1=round(get(findobj(0,'Tag','r1'),'Position'));
            temp=handles.WorkingImage.Value(rect1(1):(rect1(1)+rect1(3)),rect1(2):(rect1(2)+rect1(4)));
            im=tom_bandpass(double(temp),low,hi);
            tom_imagesc(im,'noinfo');
    end
else
    tom_setmark_imagesc(temp)%display the file
end
title('WORKING IMAGE ZOOMED','FontWeight','bold','Units','normalized');
showima=get(handles.menu_wi,'String');
handles=tom_refresh_mark(showima,param.Matrixmark,'workz',handles);
param.rectangle=rect1;
handles.Param=param;
switch Status
    case 'D'
        SetButtonDownFcn_wiz('tom_setmark(''MouseDown_DelMK'',gcbo)');
end       
guidata(h_menu,handles);
set(am,'UserData',param);

% ----------- Function MouseDown_DelMK -----------
% Called by ButtonDownFcn in WI zoomed 
% delete the marker clicked 
function MouseDown_DelMK(hObject)
h_menu=findobj(0,'Tag','menu');handles=guidata(h_menu);
am=findobj('Tag','menu');param=get(am,'Userdata');
flag=0;
showima=get(handles.menu_wi,'String');
ImageNr=str2num(showima);
handles.rectangle=get(findobj(0,'Tag','r1'),'Position'); 
point1 = get(gca,'CurrentPoint');    % button down detected
pt = point1(1,1:2); pt=round(pt);% extract x and y
Xm= pt(1)+handles.rectangle(1);
Ym=pt(2)+handles.rectangle(2);
Nbm=size(param.Matrixmark,3);
for i=1:Nbm    
    D=sqrt((param.Matrixmark(2,ImageNr,i)-Xm).^2 + (param.Matrixmark(3,ImageNr,i)-Ym).^2);
    if D<=10
        flag=1;
        x=param.Matrixmark(2,ImageNr,i)-handles.rectangle(1);y=param.Matrixmark(3,ImageNr,i)-handles.rectangle(2);
        Center= x + y*sqrt(-1);
        Radius = 10;Gridpt = 100;        
        [u,v]=setmark_circle(Center,Radius,Gridpt);
        line(u,v,'LineWidth',1,'Color',[0 1 0]);%green light
        uu = [x x x x-Radius x+Radius];vv = [y-Radius y+Radius y y y];        
        line(uu,vv,'LineWidth',1,'Color',[0 1 0]);%green light    
        text(x+25,y+3,num2str(i),'FontWeight','bold','Color',[0 1 0],'Fontsize',20);%green light
        message=['Do you really want to delete the marker number ' num2str(i) ' in the current image n???' num2str(ImageNr) '  or in all image?'];
        Question=questdlg(message,'Delete marker','Current','All','No','Current');
        if strcmp(Question,'Current') 
            tom_delmark(ImageNr,Question,param.Matrixmark,i,handles);
        elseif strcmp(Question,'All')
            message=['Are you sure you would like to delete all marker number ' num2str(i) '?'];
            Question2=questdlg(message,'Delete marker','Yes','No','No');
            if strcmp(Question2,'Yes')
                tom_delmark(ImageNr,Question,param.Matrixmark,i,handles);
            else              
               tom_drawmark(x,y,num2str(i)); 
            end
        else
            tom_drawmark(x,y,num2str(i));
        end 
    end
end
if flag==0
    message=['No marker selected !'];
    uiwait(msgbox(message,'Delete marker','help'));                
end
menu_autostop_Callback(handles.menu_autostop, '', handles);
% ----------- Function SetButtonDownFcn_wi -----------
% set the callback in ButtonDownFcn of WI
function SetButtonDownFcn_wi(fcn)
%all the text
T=findall(findobj(0,'Tag','Figwo'),'Type','text');
for i=1:size(T,1)
    set(T(i),'ButtonDownFcn',fcn);%S for single mark
end
%all the line
L=findall(findobj(0,'Tag','Figwo'),'Type','line');
for i=1:size(L,1)
    set(L(i),'ButtonDownFcn',fcn);
end
%the working image
handle_image=findall(findobj(0,'Tag','Figwo'),'Type','image');
set(handle_image,'ButtonDownFcn',fcn);
%the rectangle
set(findobj(0,'Tag','r1'),'ButtonDownFcn',fcn);

% ----------- Function SetButtonDownFcn_wiz -----------
% set the callback in ButtonDownFcn of WIZ
function SetButtonDownFcn_wiz(fcn)
%All the text
T=findall(findobj(0,'Tag','Figwo_z'),'Type','text');
for i=1:size(T,1)
    set(T(i),'ButtonDownFcn',fcn);
end
%All the line
L=findall(findobj(0,'Tag','Figwo_z'),'Type','line');
for i=1:size(L,1)
    set(L(i),'ButtonDownFcn',fcn);
end
%the working image zoomed
handle_image=findall(findobj(0,'Tag','Figwo_z'),'Type','image');
set(handle_image,'ButtonDownFcn',fcn);

% ----------- Function resetcc ------------
% reinitialize after cross correlation 
function resetcc(hObject, eventdata, handles)
set(findobj(0,'Tag','figw_z_boxinfo2'),'String','STOP');
set(findobj(0,'Tag','menu_cross'),'Value',0);
menu_cross_Callback(hObject, '', handles);
% ----------- Function cross_corr -----------
function [status,param]=cross_corr(param,handles,i,direction)
% CROSS_CORR is used to search gold marker by cross correlation
% Input
%    param: variable param
%    i: the current image
%    direction: 'up' or 'down'
% Output
%    status: error message
%    param: variable param is back
squasizex=str2num(get(handles.menu_crx,'String'));
squasizey=str2num(get(handles.menu_cry,'String'));
%get the gold mk to search
pt(1)=param.Matrixmark(2,i,param.currentmark);%
pt(2)=param.Matrixmark(3,i,param.currentmark);
Xz= round(pt(1)-(squasizex/2));Yz=round(pt(2)-(squasizey/2));
aa=[param.mypathname param.myfilename num2str(i) param.myext];
[act_im,h]=fillbymean(aa,Xz,Yz,squasizex,squasizey);%fill by the mean the 1st image
act_im.Value=tom_smooth(act_im.Value,round((squasizex*10)/100));%smooth 10% of the images's edge
if i==4
    t='pause';
end
switch direction
    case 'up'
        ptnew=prediction(param,i+1,direction);%prediction. Search where the gold maker will be
        bb=[param.mypathname param.myfilename num2str(i+1) param.myext];%read next image
    case 'down'
        ptnew=prediction(param,i-1,direction);%prediction. Search where the gold maker will be
        bb=[param.mypathname param.myfilename num2str(i-1) param.myext];%read previous image
end
%Xz=(ptnew(1)-(squasizex/2));Yz=(ptnew(2)-(squasizey/2));
Xz=ptnew(1);Yz=ptnew(2);
[oth_im,h]=fillbymean(bb,Xz-(squasizex/2),Yz-(squasizey/2),squasizex,squasizey);%fill by the mean the 2nd image
oth_im.Value=tom_smooth(oth_im.Value,round((squasizex*10)/100));%smooth 10% of the images's edge
hi=str2num(get(handles.menu_filterhi,'String')); %get filter hi
low=str2num(get(handles.menu_filterlow,'String'));%get filter hi
a=tom_bandpass(act_im.Value,low,hi);
b=tom_bandpass(oth_im.Value,low,hi);
ccf=tom_corr(a,b,'norm');%cross correlate the 2 images
p=0;
while 1
    [pic,val,M]=tom_peak(ccf,str2num(get(handles.menu_gmp,'String')));pic=pic-1;
    if val<str2num(get(handles.menu_limccf,'String'))
        status=['Cross correlation failed around image ' num2str(i) '. Limit CCF peak is too big. Please reduce it. Atual peak value is: ' num2str(val)];
        return;
    end
    switch direction
        case 'up'
            param.Matrixmark(2,i+1,param.currentmark)=(Xz)+(pic(1)-squasizex/2);
            param.Matrixmark(3,i+1,param.currentmark)=(Yz)+(pic(2)-squasizey/2);
        case 'down'
            param.Matrixmark(2,i-1,param.currentmark)=(Xz)+(pic(1)-squasizex/2);
            param.Matrixmark(3,i-1,param.currentmark)=(Yz)+(pic(2)-squasizey/2);
    end
    [pm,psi,error,x,y,z]=tom_alignment3d(param.Matrixmark,param.currentmark);
    if error<str2num(get(handles.menu_threshold,'String')) | isnan(error) | isinf(error)
        set(handles.menu_error,'String',num2str(error));
        status=[''];
        break;
    else
        ccf=M;
    end
    p=p+1;
    if p>15        
        status=['Cross correlation failed in image ' num2str(i) ' or in the next or previous one. Increase the value of the threshold or the size of CC'];
        return;
    end
end
%figure;tom_imagesc(oth_im.Value);figure;tom_imagesc(act_im.Value);figure;tom_imagesc(ccf);
%hold on;plot(squasizex/2,squasizey/2,'r+');plot(pic(1),pic(2),'g+');hold off;

%%% REFINE MARK %%%,important otherwise the result is shifted
%sq=100;%str2num(get(handles.menu_gmp,'String'))+10;%size of the marker
%switch direction
%    case 'up'
%        Xz=param.Matrixmark(2,i+1,param.currentmark)-sq/2;
%        Yz=param.Matrixmark(3,i+1,param.currentmark)-sq/2;
%    case 'down'
%        Xz=param.Matrixmark(2,i-1,param.currentmark)-sq/2;
%        Yz=param.Matrixmark(3,i-1,param.currentmark)-sq/2;
%end
%%if i==27
%%%    r='';
%%end
%[rc_im,h]=fillbymean(bb,Xz,Yz,sq,sq);
%m=mean(mean(rc_im.Value));
%mask=tom_sphere([sq sq],str2num(get(handles.menu_gmp,'String')),1,[sq/2+1 sq/2+1]);
%for ii=1:sq
%    for j=1:sq
%        if mask(ii,j)==0
%            mask(ii,j)=m;
%        end
%    end
%end
%ccf=tom_corr(mask,rc_im.Value);
%[pic,val]=tom_peak(ccf);pic=pic-1;
%if isempty(pic)
%    status=['Cross correlation failed in image ' num2str(i) ' or before. Please continue manually or change the parameter'];
%    return;
%end
%switch direction
%    case 'up'
%        param.Matrixmark(2,i+1,param.currentmark)=(Xz)+(pic(1));
%        param.Matrixmark(3,i+1,param.currentmark)=(Yz)+(pic(2));
%    case 'down'
%        param.Matrixmark(2,i-1,param.currentmark)=(Xz)+(pic(1));
%        param.Matrixmark(3,i-1,param.currentmark)=(Yz)+(pic(2));
%end
set(findobj('Tag','menu'),'UserData',param);





