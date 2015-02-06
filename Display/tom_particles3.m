function varargout = tom_particles3(varargin)
%TOM_PARTICLES3 is a tool for interactive 3D particles picking 
%
%   varargout = tom_particles3(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout   		...
%
%   tom_particles3(volume) This function can be used to pick particles from
%   tomograms interactively. The volumes is not kept in memory and he is read 
%   slice by slice from the files.
%   It is a GUI to scan through the volume in z-direction and lets you
%   adjust the contrast interactively. Additionally, a running average in
%   z-direction can be calculated to increase contrast.
%   The List Marker Point is the following:
%    Each ligne correspond to a marker number
%           ligne 1 -> marker n???1
%           ligne 2 -> marker n???2
%               .           .
%           ligne n -> marker n???n
%    The row is defined as the following
%       -row 1: particle's number
%       -row 2: picture size (x)
%       -row 3: picture size (y)
%       -row 4: picture size (z)
%       -row 5: binning of the picture
%       -row 6: x coordinate of the marker point
%       -row 7: y coordinate of the marker point
%       -row 8: z coordinate of the marker point
%       -row 9: 2*radius
%       -row 10: x coordinate of the particle's position with regard to the
%               center of the image
%       -row 11: y coordinate of the particle's position with regard to the
%               center of the image
%       -row 12: z coordinate of the particle's position with regard to the
%               center of the image
%       -row 13: 2*radius*binning
%  
%   Syntaxe: tom_particles3 or tom_particles3(1,'fixed')
%            tom_particles3(1,Filename) or tom_particles3(1,Filename,'fixed') 
%       Input :
%           - 1: You have to put 1 if you put parameter as string.
%           - Filename: name of EM-File. If it's omitted, a browser is
%           opened.
%           -'fixed': to see the real size of the image.
%                    (1 pixel of the image correspond to 1 pixel 
%                    on the display).
%           -'fixed noinfo': to see the real size of the image without
%                            image informations.
%
%EXAMPLE
%    You don't have the name of the file and want to open a browser:
%       tom_particles3;
%       tom_particles3 or tom_particles3(1,'fixed')
%       tom_particles3(1,'fixed noinfo');
%   
%    You want to give as parameter the file's name
%       tom_particles3(1,'C:\MATLAB6p5\work\mem_2bin1_vol.em');
%       tom_particles3(1,'C:\MATLAB6p5\work\mem_2bin1_vol.em','fixed');
%
%REFERENCES
%
%SEE ALSO
%   TOM_VOLXY, TOM_INTERVOL
%
%   created by WDN 11/13/03
%   updated by WDN 02/04/04
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

% error(nargchk(0, 1, nargin, 'struct'))


% Edit the above text to modify the response to help tom_particles3

% Last Modified by GUIDE v2.5 02-Nov-2009 13:21:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_particles3_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_particles3_OutputFcn, ...
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


% --- Executes just before tom_particles3 is made visible.
function tom_particles3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_particles3 (see VARARGIN)

% Choose default command line output for tom_particles3
handles.output = hObject;
% Update handles structure
handles.Option='';
if size(varargin,1) <1 
[filename, pathname] = uigetfile({'*.em;*.vol;*.norm;*.ccf';'*.*'}, 'Pick an EM-file');
    if isequal(filename,0) | isequal(pathname,0) 
        error('Cancel button pressed. No data loaded.');
        return; 
    end;
    handles.Filename= [pathname,filename];
    
else
    if size(varargin,2)==2
        if findstr(varargin{2},'fixed')
            handles.Option='fixed';
            if findstr(varargin{2},'noinfo')
                handles.Option='fixed noinfo';
            end
            [filename, pathname] = uigetfile({'*.em;*.vol;*.norm;*.ccf';'*.*'}, 'Pick an EM-file');
            if isequal(filename,0) | isequal(pathname,0) 
                error('Cancel button pressed. No data loaded.');
                return; 
            end;
            handles.Filename= [pathname,filename];            
        else   
            handles.Filename=cell2mat(varargin(2));
        end
    elseif size(varargin,2)==3
        handles.Filename=cell2mat(varargin(2));
        handles.Option=varargin{3};
    end
end

%if size(varargin,1) <1 
%[filename, pathname] = uigetfile({'*.em;*.vol';'*.*'}, 'Pick an EM-file');
%    if isequal(filename,0) | isequal(pathname,0) error('No data loaded.');
%        return; 
%    end;
%   handles.Filename= [pathname,filename];
%else
%    handles.Filename=cell2mat(varargin(2));
%end;
clear varargin
guidata(hObject,handles);
handles.Header=tom_reademheader(handles.Filename);

%%
tmpp = zeros(8,8,128);
for lauf=1:128
    rand_pos=rand(3,1);
    rand_pos=rand_pos.*handles.Header.Header.Size;
    rand_pos=rand_pos+1;
    rand_pos=floor(rand_pos);
    rand_pos(1)=abs(rand_pos(1)-10);
    rand_pos(2)=abs(rand_pos(2)-10);
    tmp=tom_emreadc(handles.Filename,'subregion', rand_pos' ,[7 7 0]);
    tmpp(:,:,lauf)=tmp.Value;
end

%%

[mean max min std]=tom_dev(tmpp);
handles.DataScale=[mean-2*std mean+2*std];
set(handles.limit_down,'String',handles.DataScale(1));
set(handles.limit_up,'String',handles.DataScale(2));
handles.HistogramData=tmpp;
dim_x=handles.Header.Header.Size(1);
dim_y=handles.Header.Header.Size(2);
dim_z=handles.Header.Header.Size(3);
handles.Image=tmp;
handles.actualaxis=[1 dim_x 1 dim_y];
handles.radius = 2;
handles.ListMarkerPoint=zeros(4,1);
handles.ListMarkerPoint(:,:)=-1;
handles.ListMarkerPoint(4,:)=2;
handles.ParticleNumber='yes';
set(handles.volume_name,'String',[handles.Filename ' ( ' num2str(dim_x) ' x ' num2str(dim_y) ' x ' num2str(dim_z) ' )']); 
guidata(hObject,handles);
clear tmpp;
% UIWAIT makes tom_particles3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = tom_particles3_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% generate all Min and Max
set(handles.XY_slider,'Min',1);
set(handles.XY_slider,'Max',handles.Header.Header.Size(3))
set(handles.XY_slider,'SliderStep',[1./handles.Header.Header.Size(3) 5./handles.Header.Header.Size(3)]);
set(handles.AVG_slider,'Min',1);
set(handles.AVG_slider,'Max',handles.Header.Header.Size(3))
set(handles.AVG_slider,'SliderStep',[1./handles.Header.Header.Size(3) 1./handles.Header.Header.Size(3)]);
updatehistogram_Callback(hObject, eventdata, handles);
set(handles.figure1,'Toolbar','figure');
prompt = {'Enter the binning relative to the full reconstruction'};
dlg_title = 'Full reconstruction''s binning';
num_lines= 1;
def     = {'2'};
answer  = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(cell2mat(answer))   
    handles.Binning=2;
else
    handles.Binning=str2num(cell2mat(answer));
end
prompt = {'Enter the z-offset used for the actual reconstruction.'};
dlg_title = 'z-offset';
num_lines= 1;
def     = {'0'};
answer_zoffset  = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(cell2mat(answer_zoffset))   
    handles.z_offset=0;
else
    handles.z_offset=str2num(cell2mat(answer_zoffset)).*(2.^(handles.Binning));
end
guidata(hObject,handles);

%***********************************************************
% --- MENU LIST ---
function file_Callback(hObject, eventdata, handles)

% --- MENU RECONSTRUCTION ---
function reconstruction_Callback(hObject, eventdata, handles)

% --- MENU LOAD LIST ---
function load_Callback(hObject, eventdata, handles)
[myname, mypathname] = uigetfile('*.em', 'LOAD LIST');
myfile=[mypathname myname]
if myfile(1)~=0 & myfile(2)~=0 %myfile= 0 0 when 'cancel' is clicked
    temp2=tom_emread(myfile);
    temp=temp2.Value';
    Nbm=size(temp,1);
    for i=1:Nbm
        handles.ListMarkerPoint(1,i)=temp(i,6);
        handles.ListMarkerPoint(2,i)=temp(i,7);
        handles.ListMarkerPoint(3,i)=temp(i,8);
        handles.ListMarkerPoint(4,i)=round(temp(i,9)/2);
    end
    RefreshImage_Callback(hObject, eventdata, handles);
    guidata(hObject,handles);
end

% --- MENU SAVE LIST ---
function save_Callback(hObject, eventdata, handles)
[myname, mypathname] = uiputfile('*.em', 'SAVE YOUR LIST AS');
myfile=[mypathname myname];
if myfile(1)~=0 & myfile(2)~=0 %myfile= 0 0 when 'cancel' is clicked
    if isempty(findstr(myfile,'.em'))
        myfile=strcat(myfile,'.em');
    end    
    for i=1:size(handles.ListMarkerPoint,2)
        temp(i,1)=i;
        temp(i,2)=handles.Header.Header.Size(1);
        temp(i,3)=handles.Header.Header.Size(2);
        temp(i,4)=handles.Header.Header.Size(3);
        temp(i,5)=handles.Binning;
        temp(i,6)=handles.ListMarkerPoint(1,i);
        temp(i,7)=handles.ListMarkerPoint(2,i);
        temp(i,8)=handles.ListMarkerPoint(3,i);
        temp(i,9)=2*handles.ListMarkerPoint(4,i);
        temp(i,10)=(temp(i,6)-((temp(i,2)/2)+1))*2^temp(i,5);
        temp(i,11)=(temp(i,7)-((temp(i,3)/2)+1))*2^temp(i,5);
        temp(i,12)=(temp(i,8)-((temp(i,4)/2)+1))*2^temp(i,5);
        temp(i,13)=temp(i,9)*2^temp(i,5);
    end
    temp2=temp';
    temp2=tom_emheader(temp2);  
    tom_emwrite(myfile,temp2);
end

% --- MENU RECONSTRUCTION PARTICULES ---
function rec_particles_Callback(hObject, eventdata, handles)
for i=1:size(handles.ListMarkerPoint,2)
    temp(i,1)=i;
    temp(i,2)=handles.Header.Header.Size(1);
    temp(i,3)=handles.Header.Header.Size(2);
    temp(i,4)=handles.Header.Header.Size(3);
    temp(i,5)=handles.Binning;
    temp(i,6)=handles.ListMarkerPoint(1,i);
    temp(i,7)=handles.ListMarkerPoint(2,i);
    temp(i,8)=handles.ListMarkerPoint(3,i);
    temp(i,9)=2*handles.ListMarkerPoint(4,i);
    temp(i,10)=(temp(i,6)-((temp(i,2)/2)+1))*2^temp(i,5);
    temp(i,11)=(temp(i,7)-((temp(i,3)/2)+1))*2^temp(i,5);
    temp(i,12)=(temp(i,8)-((temp(i,4)/2)+1))*2^temp(i,5)+handles.z_offset;
    temp(i,13)=temp(i,9)*2^temp(i,5);
    temp(i,14)=handles.z_offset;
end
tom_recparticles(temp);

% --- BUTTON SET HISTOGRAM ---
function sethistogram_Callback(hObject, eventdata, handles)
%tmp_obj=findobj('Tag','Histogram');
axes(handles.Histogram);
k = waitforbuttonpress;
     point1 = get(gca,'CurrentPoint');    % button down detected
     finalRect = rbbox;                   % return figure units
     point2 = get(gca,'CurrentPoint');    % button up detected
     point1 = point1(1,1:2);              % extract x and y
     point2 = point2(1,1:2);
     p1 = min(point1,point2);             % calculate locations
     offset = abs(point1-point2);         % and dimensions
     x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
     y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
     hold on;
     axis manual;    
handles.DataScale=[x(1) x(2)];
set(handles.limit_down,'String',handles.DataScale(1));
set(handles.limit_up,'String',handles.DataScale(2));
guidata(hObject,handles);
set(gca,'Xlim',[x(1) x(2)]);
UPDATE_ALL_Callback(hObject, eventdata, handles);
set(handles.Histogram,'Tag','Histogram');
RefreshImage_Callback(hObject, eventdata, handles);

% --- BUTTON RESET HISTOGRAM ---
function updatehistogram_Callback(hObject, eventdata, handles)
axes(handles.Histogram);
[h,n]=tom_hist3d(handles.HistogramData);
handles.DataScale=[n(1)  n(size(n,2))];
set(handles.limit_down,'String',handles.DataScale(1));
set(handles.limit_up,'String',handles.DataScale(2));
h=200.*h./(100.*handles.Header.Header.Size(1).*handles.Header.Header.Size(2).*handles.Header.Header.Size(3));
bar(n,h);axis auto;
guidata(hObject,handles);
UPDATE_ALL_Callback(hObject, eventdata, handles);
set(handles.Histogram,'Tag','Histogram');
RefreshImage_Callback(hObject, eventdata, handles);

% --- BUTTON SET MANUALLY HISTOGRAM ---
function setmanual_histogramm_Callback(hObject, eventdata, handles)
tmp_obj=findobj('Tag','Histogram');
min=str2num(get(handles.limit_down,'String'));
max=str2num(get(handles.limit_up,'String'));
handles.DataScale=[min max];
guidata(hObject,handles);
set(tmp_obj,'Xlim',[min max]);
UPDATE_ALL_Callback(hObject, eventdata, handles);
RefreshImage_Callback(hObject, eventdata, handles);

% --- EDIT BOX LIMIT MIN ---
function limit_down_Callback(hObject, eventdata, handles)

% --- EDIT BOX LIMIT MAX ---
function limit_up_Callback(hObject, eventdata, handles)

% --- BUTTON ZOOM IN ---
function zoom_in_Callback(hObject, eventdata, handles)
axes(handles.XY_slice);
k = waitforbuttonpress;
     point1 = get(gca,'CurrentPoint');    % button down detected
     finalRect = rbbox;                   % return figure units
     point2 = get(gca,'CurrentPoint');    % button up detected
     point1 = point1(1,1:2);              % extract x and y
     point2 = point2(1,1:2);
     p1 = min(point1,point2);             % calculate locations
     offset = abs(point1-point2);         % and dimensions
     x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
     y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
x=round(x);
y=round(y);
handles.actualaxis=[x(1) x(2) y(1) y(3)];
axis([x(1) x(2) y(1) y(3)]);
guidata(hObject,handles);
UPDATE_ALL_Callback(hObject, eventdata, handles);
RefreshImage_Callback(hObject, eventdata, handles);

% --- BUTTON ZOOM RESET ---
function zoo_reset_Callback(hObject, eventdata, handles)
axes(handles.XY_slice);
dim_z=handles.Header.Header.Size(3);
dim_x=handles.Header.Header.Size(1);
dim_y=handles.Header.Header.Size(2);
slicexy=round(get(handles.XY_slider,'Value'));
sliceavg=round(get(handles.AVG_slider,'Value'));
sliceavgud=round(sliceavg/2);
if slicexy >= sliceavgud; 
    slicexyud=slicexy-sliceavgud+1;
else slicexyud=1; 
    sliceavg=sliceavg-(sliceavgud-slicexy); 
end;
if slicexy+sliceavgud>dim_z; sliceavg=2*(dim_z-slicexy);end;
tmp=double(tom_emreadinc(handles.Filename,[1 1 slicexyud],[dim_x-1 dim_y-1 sliceavg-1]));
tmp=mean(tmp,3);
axes(handles.XY_slice);
if findstr(handles.Option,'fixed')
    handles.Image=tmp;
    display_real(hObject, eventdata, handles);    
else
    imagesc(tmp',[handles.DataScale]);
end
set(handles.XY_slice,'Tag','XY_slice');
handles.actualaxis=[1 dim_x 1 dim_y];
axis([1 dim_x 1 dim_y]);
guidata(hObject,handles);
RefreshImage_Callback(hObject, eventdata, handles);

% --- BUTTON PICK ---
function pickit_Callback(hObject, eventdata, handles)
m=handles.output;
set(m,'Pointer','crosshair');
k = waitforbuttonpress;   
if k==0 %mouse button press
    point1 = get(gca,'CurrentPoint');    % button down detected
    pt = round(point1(1,1:2));
    Rad = 10;
    x=pt(1);y=pt(2);
    uu = [x x x x-Rad x+Rad];
    vv = [y-Rad y+Rad y y y];
    line(uu,vv,'LineWidth',2,'color',[1 0.75 0]);%[1 0 0]red dark 
    if handles.ListMarkerPoint(1,1) & handles.ListMarkerPoint(2,1)==-1 %non existant marker
        handles.ListMarkerPoint(1,1)=x;%x coordinate
        handles.ListMarkerPoint(2,1)=y;%y coordinate
        handles.ListMarkerPoint(3,1)=str2num(get(handles.TEXT_XY,'String'));%z coordinate
        handles.ListMarkerPoint(4,1)=round(str2double(get(handles.editradius,'String')));
        switch handles.ParticleNumber
            case 'yes'
                text(x+5,y+10,'1','FontWeight','bold','Color',[1 0.75 0],'Fontsize',20);%orange
        end                
    else    %already existant marker         
        ss=size(handles.ListMarkerPoint,2);
        handles.ListMarkerPoint(1,ss+1)=x;%x coordinate
        handles.ListMarkerPoint(2,ss+1)=y;%y coordinate
        handles.ListMarkerPoint(3,ss+1)=str2num(get(handles.TEXT_XY,'String'));%z coordinate
        handles.ListMarkerPoint(4,ss+1)=round(str2double(get(handles.editradius,'String')));
        switch handles.ParticleNumber
            case 'yes'
                text(x+5,y+10,num2str(ss+1),'FontWeight','bold','Color',[1 0.75 0],'Fontsize',20);%orange
        end                
    end
    set(m,'Pointer','arrow');
end
guidata(hObject,handles);

% --- BUTTON DELETE MARK ---
function delete_Callback(hObject, eventdata, handles)
uiwait(msgbox('Please, select the marker you want to delete','Delete marker','help'));
[Xz,Yz]=ginput(1);
Nbm=size(handles.ListMarkerPoint,2);
for i=1:Nbm
    D=sqrt((handles.ListMarkerPoint(1,i)-Xz).^2 + (handles.ListMarkerPoint(2,i)-Yz).^2);
    if D<=5        
        x=handles.ListMarkerPoint(1,i);y=handles.ListMarkerPoint(2,i);
        Rad = 10;     
        uu = [x x x x-Rad x+Rad];vv = [y-Rad y+Rad y y y];        
        line(uu,vv,'LineWidth',2,'Color',[1 0 0]);%red dark
        switch handles.ParticleNumber
            case 'yes'
                text(x+5,y+10,num2str(i),'FontWeight','bold','Color',[1 0 0],'Fontsize',20);%red dark
        end                
        message=['Do you really want to delete the Current marker in red or All of them?'];
        asd=0;
        Question=questdlg(message,'Delete marker','Current','All','No','Current');
        if strcmp(Question,'Current') 
            handles.ListMarkerPoint(:,i)=-1;
            t=size(handles.ListMarkerPoint,2);
            if t~=1
                handles.ListMarkerPoint(:,i)=[];
            end
            if findstr(handles.Option,'fixed')
                handles.Image=handles.Image;
                display_real(hObject, eventdata, handles);
            else
                imagesc(handles.Image',[handles.DataScale]);
            end           
            RefreshImage_Callback(hObject, eventdata, handles);
            guidata(hObject,handles);
        elseif strcmp(Question,'All')
            message=['Are you sure you want to delete all markers number ?'];
            Question2=questdlg(message,'Delete marker','Yes','No','No');
            if strcmp(Question2,'Yes')
                handles.ListMarkerPoint=zeros(4,1);
                handles.ListMarkerPoint(:,:)=-1;
                if findstr(handles.Option,'fixed')
                    handles.Image=handles.Image;
                    display_real(hObject, eventdata, handles);
                else
                    imagesc(handles.Image',[handles.DataScale]);
                end 
            else
                if findstr(handles.Option,'fixed')
                    handles.Image=handles.Image;
                    display_real(hObject, eventdata, handles);
                else
                    imagesc(handles.Image',[handles.DataScale]);
                end                                       
                RefreshImage_Callback(hObject, eventdata, handles);
            end
        else
            if findstr(handles.Option,'fixed')
                handles.Image=handles.Image;
                display_real(hObject, eventdata, handles);
            else
                imagesc(handles.Image',[handles.DataScale]);
            end                       
            RefreshImage_Callback(hObject, eventdata, handles);
        end 
        break;
    end
end
guidata(hObject,handles);

% --- BUTTON SHOW MARK ---
function show_Callback(hObject, eventdata, handles)
    for i=1:size(handles.ListMarkerPoint,2)
        temp(i,1)=i;
        temp(i,2)=handles.Header.Header.Size(1);
        temp(i,3)=handles.Header.Header.Size(2);
        temp(i,4)=handles.Header.Header.Size(3);
        temp(i,5)=handles.Binning;
        temp(i,6)=handles.ListMarkerPoint(1,i);
        temp(i,7)=handles.ListMarkerPoint(2,i);
        temp(i,8)=handles.ListMarkerPoint(3,i);
        temp(i,9)=2*handles.ListMarkerPoint(4,i);
        temp(i,10)=(temp(i,6)-((temp(i,2)/2)+1))*2^temp(i,5);
        temp(i,11)=(temp(i,7)-((temp(i,3)/2)+1))*2^temp(i,5);
        temp(i,12)=(temp(i,8)-((temp(i,4)/2)+1))*2^temp(i,5);
        temp(i,13)=temp(i,9)*2^temp(i,5);
    end    
ListMarkerPoint=temp

% --- PARTICLES'S NUMBER ---
function checkbox1_Callback(hObject, eventdata, handles)
a=get(handles.checkbox1,'Value');
if a==0;
    handles.ParticleNumber='no';
else;
    handles.ParticleNumber='yes';
end
guidata(hObject,handles);
if findstr(handles.Option,'fixed')
    handles.Image=handles.Image;
    display_real(hObject, eventdata, handles);
else
    imagesc(handles.Image',[handles.DataScale]);
end           
RefreshImage_Callback(hObject, eventdata, handles);

% --- RADIUS ---
function editradius_Callback(hObject, eventdata, handles)
handles.radius = round(str2double(get(handles.editradius,'String')));
guidata(hObject,handles);

% --- SLIDER Z DIRECTION ---
function XY_slider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject,'SliderStep',[1./16 1./16]);
set(hObject,'Value',[1]);

% --- SLIDER Z DIRECTION ---
function XY_slider_Callback(hObject, eventdata, handles)
dim_x=handles.Header.Header.Size(1);
dim_y=handles.Header.Header.Size(2);
dim_z=handles.Header.Header.Size(3);
slicexy=round(get(hObject,'Value'));
set(handles.TEXT_XY,'String',num2str(slicexy));
sliceavg=round(get(handles.AVG_slider,'Value'));
sliceavgud=round(sliceavg/2);
if slicexy >= sliceavgud; 
    slicexyud=slicexy-sliceavgud+1;
else slicexyud=1; 
    sliceavg=sliceavg-(sliceavgud-slicexy); 
end;
rangex_min=handles.actualaxis(1);
rangex=handles.actualaxis(2)-rangex_min;
rangey_min=handles.actualaxis(3);
rangey=handles.actualaxis(4)-rangey_min;
if slicexy+sliceavgud>dim_z; sliceavg=2*(dim_z-slicexy);end;
%tmp=double(tom_emreadinc(handles.Filename,[rangex_min rangey_min slicexy],[rangex rangey sliceavg-1]));
tmp=tom_emreadc(handles.Filename,'subregion',[1 1 slicexyud],[dim_x-1 dim_y-1 sliceavg-1]);
tmp=tmp.Value;
tmp=mean(tmp,3);
axes(handles.XY_slice);
if findstr(handles.Option,'fixed')
    handles.Image=tmp;
    display_real(hObject, eventdata, handles);    
else
    imagesc(tmp',[handles.DataScale]);
end
axis(handles.actualaxis);
set(gcf,'DoubleBuffer','on');
handles.Image=tmp;
set(handles.XY_slice,'Tag','XY_slice');
RefreshImage_Callback(hObject, eventdata, handles);
guidata(hObject,handles);

% --- IMAGE NUMBER ---
function TEXT_XY_Callback(hObject, eventdata, handles)
dim_z=handles.Header.Header.Size(3);
dim_x=handles.Header.Header.Size(1);
dim_y=handles.Header.Header.Size(2);
slicexy=str2num(get(hObject,'String'));
set(handles.XY_slider,'Value',slicexy);
sliceavg=round(get(handles.AVG_slider,'Value'));
sliceavgud=round(sliceavg/2);
if slicexy >= sliceavgud; 
    slicexyud=slicexy-sliceavgud+1;
else slicexyud=1; 
    sliceavg=sliceavg-(sliceavgud-slicexy); 
end;
rangex_min=handles.actualaxis(1);
rangex=handles.actualaxis(2)-rangex_min;
rangey_min=handles.actualaxis(3);
rangey=handles.actualaxis(4)-rangey_min;
if slicexy+sliceavgud>dim_z; sliceavg=2*(dim_z-slicexy);end;
%tmp=double(tom_emreadinc(handles.Filename,[rangex_min rangey_min slicexy],[rangex rangey sliceavg-1]));
tmp=double(tom_emreadinc(handles.Filename,[1 1 slicexyud],[dim_x-1 dim_y-1 sliceavg-1]));
tmp=mean(tmp,3);
axes(handles.XY_slice);
if findstr(handles.Option,'fixed')
    handles.Image=tmp;
    display_real(hObject, eventdata, handles);    
else
    imagesc(tmp',[handles.DataScale]);
end
axis(handles.actualaxis);
set(gcf,'DoubleBuffer','on');
handles.Image=tmp;
set(handles.XY_slice,'Tag','XY_slice');
RefreshImage_Callback(hObject, eventdata, handles);
guidata(hObject,handles);

% --- SLIDER AVERAGE ---
function AVG_slider_Callback(hObject, eventdata, handles)
sliceavg=round(get(hObject,'Value'));
set(handles.TEXT_AVG,'String',num2str(sliceavg));

% --- AVERAGE NUMBER ---
function TEXT_AVG_Callback(hObject, eventdata, handles)
dim_z=handles.Header.Header.Size(3);
sliceavg=round(eval(get(handles.TEXT_AVG,'String')));
if sliceavg<1 sliceavg=1; set(handles.TEXT_AVG,'String',sliceavg); end;
if sliceavg>dim_z sliceavg=dim_z; set(handles.TEXT_AVG,'String',sliceavg); end;
set(handles.AVG_slider,'Value',sliceavg);

% --- BUTTON QUIT ---
function byebye_Callback(hObject, eventdata, handles)
delete(handles.figure1);

%********************************************************************
%*****   Other function  ********************************************
%********************************************************************

% ---------- UPDATE_ALL ----------
function UPDATE_ALL_Callback(hObject, eventdata, handles)
dim_x=handles.Header.Header.Size(1);
dim_y=handles.Header.Header.Size(2);
dim_z=handles.Header.Header.Size(3);
get(handles.TEXT_XY,'Value');
slicexy=round(eval(get(handles.TEXT_XY,'String')));
sliceavg=round(get(handles.AVG_slider,'Value'));
sliceavgud=round(sliceavg/2);
if slicexy >= sliceavgud; 
    slicexyud=slicexy-sliceavgud+1;
else slicexyud=1; 
    sliceavg=sliceavg-(sliceavgud-slicexy); 
end;
if slicexy+sliceavgud>dim_z; sliceavg=2*(dim_z-slicexy);end;
tmp=tom_emreadc(handles.Filename,'subregion',[1 1 slicexyud],[dim_x-1 dim_y-1 sliceavg-1]);
tmp=tmp.Value;
tmp=mean(tmp,3);
axes(handles.XY_slice);
if findstr(handles.Option,'fixed')
    handles.Image=tmp;
    display_real(hObject, eventdata, handles);    
else
    imagesc(tmp',[handles.DataScale]);
end
axis(handles.actualaxis);
colormap gray;
set(handles.XY_slice,'Tag','XY_slice');
slicexy=round(get(handles.XY_slider,'Value'));
axes(handles.XY_slice);

% ---------- RefreshImage ----------
function RefreshImage_Callback(hObject, eventdata, handles)
if handles.ListMarkerPoint(1:2,:)==-1
    %no action
else
    ss=size(handles.ListMarkerPoint);
    for i=1:ss(2)
        refmark=handles.ListMarkerPoint(3,i);
        radius=handles.ListMarkerPoint(4,i);
        currimage = str2num(get(handles.TEXT_XY,'String'));
        if currimage-refmark>0
            if refmark+radius>=currimage
                x=handles.ListMarkerPoint(1,i);
                y=handles.ListMarkerPoint(2,i);
                Rad = 10;
                uu = [x x x x-Rad x+Rad];
                vv = [y-Rad y+Rad y y y];
                line(uu,vv,'LineWidth',1,'color',[0 1 0]);%green dark
                switch handles.ParticleNumber
                    case 'yes'
                        text(x+5,y+10,num2str(i),'FontWeight','bold','Color',[0 1 0],'Fontsize',20);%green dark
                end
            end
        elseif currimage-refmark<0
            if refmark-radius<=currimage
                x=handles.ListMarkerPoint(1,i);
                y=handles.ListMarkerPoint(2,i);
                Rad = 10;
                uu = [x x x x-Rad x+Rad];
                vv = [y-Rad y+Rad y y y];
                line(uu,vv,'LineWidth',1,'color',[0 1 0]);%green dark
                switch handles.ParticleNumber
                    case 'yes'
                        text(x+5,y+10,num2str(i),'FontWeight','bold','Color',[0 1 0],'Fontsize',20);%green dark
                end
            end
        elseif currimage-refmark==0
            x=handles.ListMarkerPoint(1,i);
            y=handles.ListMarkerPoint(2,i);
            Rad = 10;
            uu = [x x x x-Rad x+Rad];
            vv = [y-Rad y+Rad y y y];
            line(uu,vv,'LineWidth',2,'color',[1 0.75 0]);%orange
            switch handles.ParticleNumber
                case 'yes'
                    text(x+5,y+10,num2str(i),'FontWeight','bold','Color',[1 0.75 0],'Fontsize',20);%orange
            end
        end
    end
end

% ---------- Display the image with 'fixed'option ----------
function display_real(hObject, eventdata, handles);
in=handles.Image';
in_red=imresize(in,.25);
imagesc(in,[handles.DataScale]);
colormap gray(256);
param1='';param2='';
if findstr(handles.Option,'fixed')
    param1='fixed';
end
if findstr(handles.Option,'noinfo')
    param2='noinfo';
end
switch param1        
     case 'fixed'
        set(gca,'Units','pixels');
        pp=get(gca,'Position');sf=size(in);            
        set(gca,'Position',[pp(1) pp(2) sf(1) sf(2)]);
        if isempty(param2)
            t=title(['Info: (' num2str(size(in)) ') pixel']); %, mean:' num2str(meanv,3) ' std:' num2str(std,3)
        end
    otherwise
        axis image; axis ij; %colormap gray; %nothing changed, as nargin=1
        if isempty(param2)
            t=title(['Info: (' num2str(size(in)) ') pixel']);
        end
end


% --------------------------------------------------------------------
function Motl_Callback(hObject, eventdata, handles)
% hObject    handle to Motl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ListMarkerPoint=[];
[fname path]=uigetfile( ...
       {'*.em'}, ...
        'Pick a MOTL file');;
list=tom_emread([path fname]);
 h = waitbar(0,'LI LA LOADING...');

for i=1:size(list.Value,2)
     handles.ListMarkerPoint(1,i)=list.Value(8,i);%x coordinate
     handles.ListMarkerPoint(2,i)=list.Value(9,i);;%y coordinate
     handles.ListMarkerPoint(3,i)=list.Value(10,i);%z coordinate
     handles.ListMarkerPoint(4,i)=round(str2double(get(handles.editradius,'String')));
    waitbar(i/size(list.Value,2),h);
end;
close(h);

guidata(hObject,handles);

