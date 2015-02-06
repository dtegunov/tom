function varargout = tom_embrowse(varargin)
%TOM_EMBROWSE is an interactive tool to view EM-Format and MRC-Format files
%
%   varargout = tom_embrowse(varargin)
%
%   With tom_embrowse, the user can easily have a preview visualization of
%   EM-Format and MRC-Format files. The user just have to select in the list
%   boxes the file he wants to visualize. 
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
%   tom_embrowse(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_ISEMFILE, TOM_ISMRCFILE, tom_emreadc, TOM_MRCREAD
%
%   created by William Del Net 01/15/03
%   updated by William Del Net 07/09/05
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
                   'gui_OpeningFcn', @tom_embrowse_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_embrowse_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before tom_embrowse is made visible.
function tom_embrowse_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_embrowse (see VARARGIN)

% Choose default command line output for tom_embrowse
handles.output = hObject;
a=['byte    ';'short   ';'        ';'long int';'float   ';,...
   '        ';'        ';'complex '; 'double  '];
handles.DataType=cellstr(a);
handles.InitialPath=pwd;
handles.Path=pwd;
handles.Filename='';
handles.rnb1=1;
handles.rnb2=1;
% Update handles structure
guidata(hObject, handles);

if nargin == 3,
    initial_dir = pwd;
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
% Populate the listbox
load_listbox(initial_dir,handles)
% Return figure handle as first output argument

% UIWAIT makes tom_embrowse wait for user response (see UIRESUME)
% uiwait(handles.tom_embrowse);


% --- Outputs from this function are returned to the command line.
function varargout = tom_embrowse_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ----------------------------------------------------------
%
%----------------------
% --- BUTTON BROWSE ---
%----------------------
function browse_Callback(hObject, eventdata, handles)
p = uigetdir;
if isequal(p,0)
    %nothing because Cancel is ckicked
else
	load_listbox(p,handles);
	handles.Path=p;
	handles.Filename='';
	guidata(hObject,handles);
end
%--------------------------
% --- BUTTON NEW FOLDER ---
%--------------------------
function NewFolder_Callback(hObject, eventdata, handles)
prompt='Enter a folder name';
Title='Folder name';
answer = inputdlg(prompt);
if ~isempty(answer)
    answer = [handles.Path '/' answer{1}];
    ex=exist(answer,'dir');
    if ex==7
        message='Error! This folder name exist already.';
        uiwait(msgbox(message,'Error','error'));
    else
        mkdir(answer);
        load_listbox(handles.Path,handles);
    end
end
%---------------------------
% --- LIST BOX DIRECTORY ---
%---------------------------
function listbox_dir_Callback(hObject, eventdata, handles)
InitField(handles);
index_selected = get(handles.listbox_dir,'Value');
file_list = get(handles.listbox_dir,'String');
filename = [handles.Path '/' file_list{index_selected}];
file= file_list{index_selected};
[path,name,ext] = fileparts(filename);
ver='';

switch file
    case '.'
        filename = handles.Path;
        load_listbox(filename,handles);
    case '..'
        cd (handles.Path);
        cd ('..');
        filename=pwd;
        load_listbox(filename,handles);
    otherwise
        load_listbox(filename,handles);
end
handles=guidata(handles.tom_embrowse);
handles.Filename='';
disp_im(0);

%----------------------
% --- BUTTON DELETE ---
%----------------------
function delete_Callback(hObject, eventdata, handles)
index_selected = get(handles.listbox_filename,'Value');
file_list = get(handles.listbox_filename,'String');
message='Do you really want to delete the selected file?';
Question2=questdlg(message,'Delete file','Yes','No','No');
if strcmp(Question2,'Yes')
    for i=1:size(index_selected,2)
        filename = [handles.Path '/' file_list{index_selected(i)}];
        delete(filename);
    end
    load_listbox(handles.Path,handles);
end
%--------------------
% --- BUTTON COPY ---
%--------------------
function copy_file_Callback(hObject, eventdata, handles)
cd (handles.Path);
index_selected = get(handles.listbox_filename,'Value');
file_list = get(handles.listbox_filename,'String');
setappdata(0,'UseNativeSystemDialogs',0);
dirname = uigetdir(handles.Path, 'Copy file');
if dirname~=0 % when 'cancel' is clicked
    for i=1:size(index_selected,2)
        filename = [handles.Path '/' file_list{index_selected(i)}];
        copyfile(filename,dirname,'f')
    end
    load_listbox(handles.Path,handles);
end
cd (handles.InitialPath);
%--------------------
% --- BUTTON MOVE ---
%--------------------
function move_file_Callback(hObject, eventdata, handles)
cd (handles.Path);
index_selected = get(handles.listbox_filename,'Value');
file_list = get(handles.listbox_filename,'String');
setappdata(0,'UseNativeSystemDialogs',0);
dirname = uigetdir(handles.Path, 'Move file');
if dirname~=0 % when 'cancel' is clicked
    for i=1:size(index_selected,2)
        filename = [handles.Path '/' file_list{index_selected(i)}];
        movefile(filename,dirname,'f')
    end
    load_listbox(handles.Path,handles);
end
cd (handles.InitialPath);
%--------------------------
% --- LIST BOX FILENAME ---
%--------------------------
function varargout = listbox_filename_Callback(hObject, eventdata, handles)
if strcmp(get(handles.tom_embrowse,'SelectionType'),'normal')%'open'
    index_selected = get(handles.listbox_filename,'Value');
    if size(index_selected,2)>1
        return
    end
    file_list = get(handles.listbox_filename,'String');
    filename = [handles.Path '/' file_list{index_selected}];
    file= file_list{index_selected};
    [path,name,ext] = fileparts(filename);
    ver='';
    
    %---- EM Format file ----
    if tom_isemfile(filename)==1
        info=tom_reademheader(filename);
        set(handles.name,'String',file);
        set(handles.size,'String',[num2str(info.Header.Size(1)) ' x ' num2str(info.Header.Size(2)) ' x ' num2str(info.Header.Size(3))]);
        set(handles.angle,'String',[info.Header.Parameter(19)./1000.0]);
        set(handles.tiltaxis,'String',info.Header.Tiltaxis);
        set(handles.defocus,'String',info.Header.Defocus);
        set(handles.objpixsize,'String',info.Header.Objectpixelsize);
        set(handles.comment,'String',char(info.Header.Comment));
        if info.Header.Size(3)>1%for 3D images
            i=tom_emreadc(filename,'subregion',[1 1 round(info.Header.Size(3)/2)],[info.Header.Size(1)-1 info.Header.Size(2)-1 0]);
        else%for 2D images
            if info.Header.Magic(1)==5 | info.Header.Magic(1)==0 | info.Header.Magic(1)==3 %OS9, mac, SGI
                i=tom_emreadc(filename);
                s=size(i.Value);
                f=s(1)/256;
                if f>1%resize if image's size >256
                    i.Value=i.Value(1:f:i.Header.Size(1),1:f:i.Header.Size(2),1);
                end
            else %option 'resample' doesn't work with images from a mac
                resx=round(info.Header.Size(1)/256);
                resy=round(info.Header.Size(2)/256);
                if resx==0 |resy==0
                    resx=1;resy=1;
                end                 
                i=tom_emreadc(filename,'resample',[resx resy 1]);
            end
        end
        i_magic=i.Header.Magic(4);
        set(handles.type,'String',handles.DataType(i_magic));
        s=size(i.Value);
        if info.Header.Size(1)==12%markerfile
            disp_im(0);
            set(handles.comment,'String',['It can not be displayed  because ' file ' is a Marker file']);
        else
            disp_im(i.Value);
        end
        handles.Filename=file;
    %---- MRC Format file ----    
    elseif tom_ismrcfile(filename)==1
        i=tom_mrcread(filename);
        set(handles.name,'String',file);
        set(handles.size,'String',[num2str(i.Header.Size(1)) ' x ' num2str(i.Header.Size(2)) ' x ' num2str(i.Header.Size(3))]);
        set(handles.angle,'String',i.Header.Tiltangle); 
        set(handles.tiltaxis,'String',i.Header.Tiltaxis); 
        set(handles.defocus,'String',i.Header.Defocus);
        set(handles.objpixsize,'String',i.Header.Objectpixelsize);
        ff=i.Value;qq=whos('ff');
        set(handles.type,'String',qq.class);
        if i.Header.Size(3)>1%for volume
            i.Value=i.Value(:,:,round(i.Header.Size(3)/2));
        end
        s=size(i.Value);
        f=round(s(1)/256);
        if f>=1%resize if image's size >256
            i.Value=i.Value(1:f:i.Header.Size(1),1:f:i.Header.Size(2),1);
            disp_im(i.Value);
        else
            disp_im(i.Value);
        end

        handles.Filename=file;


        %---- Spider Format file ----    added by SN
    elseif tom_isspiderfile(filename)==1
        i=tom_spiderread(filename);
        set(handles.name,'String',file);
        set(handles.size,'String',[num2str(i.Header.Size(1)) ' x ' num2str(i.Header.Size(2)) ' x ' num2str(i.Header.Size(3))]);
        set(handles.angle,'String',i.Header.Spider.phi);
        set(handles.tiltaxis,'String',i.Header.Spider.theta);
        %        set(handles.defocus,'String',i.Header.Defocus);
        %        set(handles.objpixsize,'String',i.Header.Objectpixelsize);
        ff=i.Value;qq=whos('ff');
        set(handles.type,'String',qq.class);
        if i.Header.Size(3)>1%for volume
            i.Value=i.Value(:,:,round(i.Header.Size(3)/2));
        end
        s=size(i.Value);
        f=round(s(1)/256);
        if f>=1%resize if image's size >256
            i.Value=i.Value(1:f:i.Header.Size(1),1:f:i.Header.Size(2),1);
            disp_im(i.Value);
        else
            disp_im(i.Value);
        end

        handles.Filename=file;

    %----imagic format ------ added by fb 
    elseif tom_isimagicfile(filename)
        i=tom_imagicread(filename);
        if i.Header.Size(3)>1%for volume
            i.Value=i.Value(:,:,round(i.Header.Size(3)/2));
        end
        set(handles.name,'String',file);
        set(handles.size,'String',[num2str(i.Header.Size(1)) ' x ' num2str(i.Header.Size(2)) ' x ' num2str(i.Header.Size(3))]);
        set(handles.objpixsize,'String',i.Header.Objectpixelsize);
        ff=i.Value;qq=whos('ff');
        set(handles.type,'String',qq.class);
        if i.Header.Size(3)>1%for volume
            i.Value=i.Value(:,:,round(i.Header.Size(3)/2));
        end
        s=size(i.Value);
        f=round(s(1)/256);
        if f>=1%resize if image's size >256
            i.Value=i.Value(1:f:i.Header.Size(1),1:f:i.Header.Size(2),1);
            disp_im(i.Value);
        else
            disp_im(i.Value);
        end
        
        handles.Filename=file;
    
    elseif tom_isxmippsell(filename) 
        i=tom_xmippsellread(filename);
        set(handles.name,'String',file);
        set(handles.size,'String',[num2str(i.Header.Size(1)) ' x ' num2str(i.Header.Size(2)) ' x ' num2str(i.Header.Size(3))]);
         ff=i.Value;qq=whos('ff');
        set(handles.type,'String',qq.class);
        if i.Header.Size(3)>1%for volume
            i.Value=i.Value(:,:,round(i.Header.Size(3)/2));
        end
        s=size(i.Value);
        f=round(s(1)/256);
        if f>=1%resize if image's size >256
            i.Value=i.Value(1:f:i.Header.Size(1),1:f:i.Header.Size(2),1);
            disp_im(i.Value);
        else
            disp_im(i.Value);
        end
        
        handles.Filename=file;
        
        %---- Directory ---- 
    elseif isdir(filename)
        InitField(handles);
        switch file
            case '.'
                filename = handles.Path;
                load_listbox(filename,handles);
            case '..'
                cd (handles.Path);
                cd ('..');
                filename=pwd;
                load_listbox(filename,handles);
            otherwise
                load_listbox(filename,handles);                
        end
        handles=guidata(handles.tom_embrowse);
        handles.Filename='';
        disp_im(0);
    else
        if ~isempty(findstr(filename,'.bmp'))|~isempty(findstr(filename,'.tif'))|~isempty(findstr(filename,'.tiff'))|~isempty(findstr(filename,'.png'))|~isempty(findstr(filename,'.jpeg'))|~isempty(findstr(filename,'.jpg'))
            i.Value= imread(filename);
            i.Header=imfinfo(filename);
            s=size(i.Value);
            if s(1)>256 | s(2)>256
                if s(1)>=s(2)
                    fact=s(1)/256;
                    mrows=round(s(1)/fact);
                    ncols=round(s(2)/fact);
                else
                    fact=s(2)/256;
                    mrows=round(s(1)/fact);
                    ncols=round(s(2)/fact);
                end
                i.Value=imresize(i.Value,[mrows ncols],'nearest');
            end
            imagesc(i.Value(:,:,1:3));
            axis image;
            set(handles.name,'String',file);
            set(handles.size,'String',[num2str(i.Header.Height) ' x ' num2str(i.Header.Width) ' x 1']);
            set(handles.angle,'String',''); 
            set(handles.tiltaxis,'String',''); 
            set(handles.defocus,'String','');
            set(handles.objpixsize,'String','');            
            axis image;
            handles.Filename=file;
        else 
            handles.Filename=file;
            InitField(handles);
            disp_im(0);      
        end
    end    
    dd=dir(filename);
    set(handles.LastFileMod,'String',dd.date); 
    guidata(hObject,handles);
%elseif strcmp(get(handles.tom_embrowse,'SelectionType'),'extend')%'open'
%    index_selected = get(handles.listbox_filename,'Value');
%    file_list = get(handles.listbox_filename,'String');
%    filename = [file_list{index_selected}];
end

%--------------------
% --- BUTTON VIEW ---
%--------------------
function view_Callback(hObject, eventdata, handles)
filename=[handles.Path '/' handles.Filename];
%---- EM Format file ----
if tom_isemfile(filename)==1
    i=tom_emreadc(filename);
    %    i.Value=tom_norm(double(i.Value),1);
    i.Header.Format='em';
    if i.Header.Size(1)==12%markerfile
        return;
    end
    %---- MRC Format file ----
elseif tom_ismrcfile(filename)==1
    i=tom_mrcread(filename);
    %    i.Value=tom_norm(double(i.Value),1);
    i.Header.Format='mrc';
    %---- SPIDER Format file ----
elseif tom_isspiderfile(filename)==1
    i=tom_spiderread(filename);
    %    i.Value=tom_norm(double(i.Value),1);
    i.Header.Format='spi';
elseif tom_isimagicfile(filename)==1
    i=tom_imagicread(filename);
    %    i.Value=tom_norm(double(i.Value),1);
    i.Header.Format='ima';
    i.Header.Filename=filename;
elseif tom_isxmippsell(filename)==1  
    i=tom_xmippsellread(filename);
    %    i.Value=tom_norm(double(i.Value),1);
    i.Header.Format='xse';
    i.Header.Filename=filename;

else
    i.Value=imread(filename);
    i.Header=imfinfo(filename);
    i.Header.Size=size(i.Value);
    if size(i.Header.Size,2)<3;i.Header.Size(3)=1;end;
end
if i.Header.Size(3)==1 %2D image
    [mean, max, min, std, variance] = tom_devinternal(i.Value);
    v=version;
    if findstr(v,'6.5.0')==1%Matlab version 6.5.0
        message='Can not display the image. Update to Matlab R14 (version 7.0.4)';
        msgbox(message,'Update Matlab','warn');
    elseif findstr(v,'6.5.1')==1%Matlab version 6.5.1
        message='Can not display the image. Update to Matlab R14 (version 7.0.4)';
        msgbox(message,'Update Matlab','warn');
    else %Matlab version 7.0.1 or higher
        if (mean-3*std)>=(mean+3*std)
            vv=imtool(i.Value','InitialMagnification','fit');
        else
            vv=imtool(i.Value',[mean-(3*std) mean+(3*std)],'InitialMagnification','fit');
        end
        set(vv,'NAME',[handles.Filename '  ' 'mean:' num2str(mean,3) ' std:' num2str(std,3) ]);
        r=get(0);
        if (isunix)
            if r.MonitorPositions(1,4) > 1024
                set(vv,'Position',[r.MonitorPositions(1,3)-1024 r.MonitorPositions(1,4)-1024 1024 1024]);
            end
        else
            if r.MonitorPositions(1,4) > 1024
                set(vv,'Position',[r.MonitorPositions(1,3)-1024 r.MonitorPositions(1,4)-1024 768 768]);
            end
        end
    end
elseif i.Header.Size(3)>1
    switch i.Header.Format
        case {'png', 'jpg', 'bmp', 'tif'}%2D image
            vv=imtool(i.Value,'InitialMagnification','fit');
        otherwise
            if get(handles.dspcub_view,'Value')
                direction = str2num(get(handles.edit_dspcub_direction,'String'));
                if tom_isxmippsell(filename)==1  
                    imtool(tom_norm(tom_gallery(i.Value)',1));
                else
                    figure;tom_dspcub(i.Value,direction);%3D image
                end
            else
                tom_volxyz(i);%3D image
            end
    end   
end
%----------------------------------
% --- BUTTON OPEN A TEXT EDITOR ---
%----------------------------------
function text_editor_Callback(hObject, eventdata, handles)
filename=[handles.Path '/' handles.Filename];
edit(filename);
%---------------------------------
% --- BUTTON MOVE_FILE TO WORKSPACE ---
%---------------------------------
function workspace_Callback(hObject, eventdata, handles)
filename=[handles.Path '/' handles.Filename];
rnb_vol=1; rnb_img=1;
%---- EM Format file ----
if tom_isemfile(filename)==1
    i=tom_emreadc(filename);i.Value=double(i.Value);
    if i.Header.Size(1)==12%markerfile
        return;
    end    
%---- MRC Format file ----
elseif tom_ismrcfile(filename)==1
    i=tom_mrcread(filename);
%---- Spider Format file ----
elseif tom_isspiderfile(filename)==1
    i=tom_spiderread(filename);
 
  %---- imagic Format file ----
elseif tom_isimagicfile(filename)==1
    i=tom_imagicread(filename);  
    %---- Image ----    
elseif tom_isxmippsell(filename)==1
    i=tom_xmippsellread(filename);  


else
    i.Value=imread(filename);
    i.Header=imfinfo(filename);
    i.Header.Size=size(i.Value);
    if size(i.Header.Size)<3;i.Header.Size(3)=1;end;
end
var_base = cellstr(evalin('base','who'));
if i.Header.Size(3)>1;var_exp='vol_';%volume
else;var_exp='img_';%image
end
if ~isempty(var_base)
    while 1
        a=mean(strcmp(var_base,[var_exp num2str(rnb_vol)]));
        if a==0;break;
        else;rnb_vol=rnb_vol+1;
        end
    end
end
var1=[var_exp num2str(rnb_vol)];
assignin('base',var1,i);
guidata(hObject, handles);

%---------------------
% --- BUTTON PRINT ---
%---------------------
function print_Callback(hObject, eventdata, handles)
filename=[handles.Path '/' handles.Filename];
if tom_isemfile(filename)
    i=tom_emreadc(filename);i.Value=double(i.Value);
    if i.Header.Size(1)==12%markerfile
        return;
    end    
elseif tom_ismrcfile(filename)
    i=tom_mrcread(filename);i.Value=double(i.Value);
else
    i.Value=imread(filename);i.Value=double(i.Value);
    i.Header=imfinfo(filename);
    i.Header.Size=size(i.Value);
    if size(i.Header.Size)<3;i.Header.Size(3)=1;end;
end
if i.Header.Size(3)>1
    i.Value=i.Value(:,:,round(i.Header.Size(3)/2));
end
[mean, max, min, std, variance] = tom_devinternal(i.Value);
m=figure('Papertype','A4','Visible','off');set(0,'CurrentFigure',m);
set(m,'NumberTitle','off','Name',filename);
if (mean-3*std)>=(mean+3*std)               
    imshow(i.Value',[mean-(3*std) mean+(3*std)],'InitialMagnification','fit');              
else
    %imagesc(i.Value',[mean-(3*std) mean+(3*std)]);
    imshow(i.Value',[mean-(3*std) mean+(3*std)],'InitialMagnification','fit');
    %imtool(i.Value',[mean-(3*std) mean+(3*std)],'InitialMagnification','fit');
end  
colormap gray; axis image;
t=title('matlab');
set(t,'Interpreter','none',...
'String',[ filename ', Info: (' num2str(size(i.Value)) '), mean:' num2str(mean,3) ', std:' num2str(std,3) ]);
print (m) ;
close(m);disp('document printed');

%----------------------
% --- BUTTON EXPORT ---
%----------------------
function export_Callback(hObject, eventdata, handles)
filename = fullfile(handles.Path,handles.Filename);
ext = {'*.tif','Tagged Image File Format (*.tif)';,...
       '*.jpg', 'Joint Photographic Experts Group (*.jpg)';...
       '*.bmp','Windows Bitmap (*.bmp)';,...
       '*.png','Portable Network Graphics (*.png)'};
varp=findstr(handles.Filename,'.');
if size(varp,2)>1
    varp=varp(size(varp,2));
end
a=size(handles.Filename,2)-varp;
if a<=3;
    weme=handles.Filename(1:varp-1);
else
    weme=handles.Filename;
end
cd(handles.Path);
[myfile,mypathname,index] = uiputfile(ext,'Save File Name As',weme);
if isequal(myfile,0) | isequal(mypathname,0)
    return;%nothing because Cancel is ckicked
else
    copyto=fullfile(mypathname,myfile);
    %---- EM Format file ----
    if tom_isemfile(filename)==1
        i=tom_emreadc(filename);    
        if i.Header.Size(1)==12%markerfile
            return;
        end
        %---- MRC Format file ----
    elseif tom_ismrcfile(filename)==1
        i=tom_mrcread(filename);i.Value=double(i.Value);
    else
        i.Value=imread(filename);i.Value=double(i.Value);
        i.Header=imfinfo(filename);
        i.Header.Size=size(i.Value);
        if size(i.Header.Size)<3;i.Header.Size(3)=1;end;
    end
end
if i.Header.Size(3)>1
    i.Value=i.Value(:,:,round(i.Header.Size(3)/2));
end
if findstr(copyto,'.tif')
    a=tom_norm(i.Value,1);
    imwrite(a',copyto,'Compression','none');%no compression
elseif findstr(copyto,'.jpg')
    a=tom_norm(i.Value,1);
    imwrite(a',copyto,'Quality',50);%best quality 100
elseif findstr(copyto,'.bmp')
    a=tom_norm(i.Value,1);
    imwrite(a',copyto);
elseif findstr(copyto,'.png')
    a=tom_norm(i.Value,1);
    imwrite(double(a)',copyto,'String',filename,'BitDepth',16);%better quality 16
end
cd(handles.InitialPath);

%---------------------------
% --- BUTTON EDIT HEADER ---
%---------------------------
function edit_Callback(hObject, eventdata, handles)
filename = fullfile(handles.Path,handles.Filename);
if tom_isemfile(filename)
    i=tom_emreadc(filename);i.Value=double(i.Value);
    if i.Header.Size(1)==12%markerfile
        return;
    end
    tom_editemheader(i);
end

%----------------------------------------
% --- BUTTON CHANGE CURRENT DIRECTORY ---
%----------------------------------------
function change_dir_Callback(hObject, eventdata, handles)
cd (handles.Path);
directory=pwd
ls;

%******************************************************
%            fonction need by tom_embrowse
%******************************************************            

% ----------- Function disp_im -----------
function disp_im(in)%disp_im(in,parameter)
in_red=imresize(double(in),.1);
[meanv max min std]=tom_devinternal(in_red);
%changed 070509 by RK (preview in tom_embrowse looks terrible)
%if (meanv-4*std)>=(meanv+4*std)
%    imagesc(in');
%else
%    imagesc(in',[meanv-4*std meanv+4*std]);colormap gray;axis image;
%end;
imagesc (in');
colormap gray;   
axis image; axis ij; 

% ----------- Function tom_devinternal -----------
function [a,b,c,d,e]=tom_devinternal(A);
A=double(A); % change by SN for Linux
[s1,s2,s3]=size(A);
a=sum(sum(sum(A)))/(s1*s2*s3);
b=max(max(max(A)));
c=min(min(min(A)));
d=std2(A);
e=d^2;

% ----------- Function load_listbox -----------
function load_listbox(dir_path,handles)
cc=1;dd=1;aa=cell(1);bb=cell(1);
cd (dir_path);
dir_struct = dir(dir_path);
a={dir_struct.name};
b={dir_struct.isdir};
for i=1:size(dir_struct,1)
    if b{i}==1
        aa{cc}=a{i};
        cc=cc+1;
    else
        bb{dd}=a{i};
        dd=dd+1;
    end
end
handles.dir=aa';%char(aa);
handles.file_names=bb';%char(bb);
handles.Path=pwd;
guidata(handles.tom_embrowse,handles);
set(handles.listbox_dir,'String',handles.dir,...
                        'Value',1,...
                        'FontWeight','bold');
set(handles.listbox_filename,'String',handles.file_names,'Value',1);
set(handles.text1,'String',pwd);
cd (handles.InitialPath);

function InitField(handles);
set(handles.name,'String','');
set(handles.size,'String','');
set(handles.type,'String','');
set(handles.angle,'String','');
set(handles.defocus,'String','');
set(handles.tiltaxis,'String','');
set(handles.objpixsize,'String','');
set(handles.LastFileMod,'String','');
set(handles.comment,'String','');





function edit_dspcub_direction_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dspcub_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dspcub_direction as text
%        str2double(get(hObject,'String')) returns contents of edit_dspcub_direction as a double


% --- Executes during object creation, after setting all properties.
function edit_dspcub_direction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dspcub_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


