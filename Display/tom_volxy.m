function varargout = tom_volxy(varargin)
%TOM_VOLXY is a 3D visualization tool for tomograms. 
%   This function can be
%   used to view a volume directly from the file. It is a GUI to scan
%   through the volume in z-direction and lets you adjust the contrast 
%   interactively. Additionally, a running average in z-direction can be 
%   calculated to increase the contrast.
%
%   varargout = tom_volxy(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout   		...
%
%   Syntax: tom_volxy or tom_volxy(1,'fixed')
%           tom_volxy(1,Filename) or tom_volxy(1,Filename,'fixed') 
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
%       tom_volxy;
%       tom_volxy or tom_volxy(1,'fixed')
%       tom_volxy(1,'fixed noinfo');
%   
%    You want to give as parameter the file's name
%       tom_volxy(1,'C:\MATLAB6p5\work\mem_2bin1_vol.em');
%       tom_volxy(1,'C:\MATLAB6p5\work\mem_2bin1_vol.em','fixed');
%
%REFERENCES
%
%SEE ALSO
%   TOM_INTERVOL, TOM_PARTICLES
%
%   created by SN 12/18/2002
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


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_volxy_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_volxy_OutputFcn, ...
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


% --- Executes just before tom_volxy is made visible.
function tom_volxy_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_volxy (see VARARGIN)

% Choose default command line output for tom_volxy
handles.output = hObject;
% Update handles structure
handles.Option='';
if size(varargin,1) <1 
[filename, pathname] = uigetfile({'*.em;*.vol';'*.*'}, 'Pick an EM-file');
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
            [filename, pathname] = uigetfile({'*.em;*.vol';'*.*'}, 'Pick an EM-file');
            if isequal(filename,0) | isequal(pathname,0) 
                error('Cancel button pressed. No data loaded.');
                return; 
            end;
            handles.Filename= [pathname,filename];            
        else   
            handles.Filename=cell2mat(varargin(2));
            filename=handles.Filename;
        end
    elseif size(varargin,2)==3        
        handles.Filename=cell2mat(varargin(2));
        filename=handles.Filename;
        handles.Option=varargin{3};
    end
end
clear varargin
guidata(hObject,handles);
handles.Header=tom_reademheader(handles.Filename);
for lauf=1:500
    rand_pos=rand(3,1);
    rand_pos=rand_pos.*handles.Header.Header.Size;
    rand_pos=rand_pos+1;
    rand_pos=floor(rand_pos);
    rand_pos(1)=abs(rand_pos(1)-10);
    rand_pos(2)=abs(rand_pos(2)-10);
    tmp=double(tom_emreadinc(handles.Filename,[rand_pos],[8 8 0]));
    tmpp(:,:,lauf)=tmp;
end;
[mean max min std]=tom_dev(tmpp);
handles.DataScale=[mean-2*std mean+2*std];
set(handles.limit_down,'String',handles.DataScale(1));
set(handles.limit_up,'String',handles.DataScale(2));
handles.Histogram=tmpp;
dim_x=handles.Header.Header.Size(1);
dim_y=handles.Header.Header.Size(2);
dim_z=handles.Header.Header.Size(3);
mdimz=abs(dim_z/2);
%set(handles.TEXT_XY,'String','128');
%set(handles.XY_slider,'Value',mdimz);
%set(handles.volume_name,'String',filename);
set(handles.volume_name,'String',[filename ' ( ' num2str(dim_x) ' x ' num2str(dim_y) ' x ' num2str(dim_z) ' )']); 
set(gcf,'Toolbar','figure');
if dim_x~=dim_y
    set(handles.XY_slice,'Units','Pixel');
    dim=get(handles.XY_slice,'Position');    
    if dim_x>dim_y
        pr=dim(3)/dim_x;
        dim(4)=pr*dim_y;
        set(handles.XY_slice,'Position',dim);
    else
        pr=dim(4)/dim_y;
        dim(3)=pr*dim_x;
        set(handles.XY_slice,'Position',dim);        
    end  
    set(handles.XY_slice,'Units','Normalized');
end
handles.dimz=mdimz;
handles.actualaxis=[1 dim_x 1 dim_y];
handles.Option_Meas='Line';
handles.rnb=1;
guidata(hObject,handles);
clear tmpp;
% UIWAIT makes tom_volxy wait for user response (see UIRESUME)
% uiwait(handles.volxy);

% --- Outputs from this function are returned to the command line.
function varargout = tom_volxy_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% generate all Min and Max
tmp_obj=findobj('Tag','XY_slider');
set(tmp_obj,'Min',1);
set(tmp_obj,'Max',handles.Header.Header.Size(3))
set(tmp_obj,'SliderStep',[1./handles.Header.Header.Size(3) 5./handles.Header.Header.Size(3)]);

set(tmp_obj,'Value',handles.dimz);
set(handles.TEXT_XY,'String',num2str(handles.dimz));

tmp_obj=findobj('Tag','AVG_slider');
set(tmp_obj,'Min',1);
set(tmp_obj,'Max',handles.Header.Header.Size(3))
set(tmp_obj,'SliderStep',[1./handles.Header.Header.Size(3) 1./handles.Header.Header.Size(3)]);
set(handles.TEXT_AVG,'String','1');
set(tmp_obj,'Value',1);
updatehistogram_Callback(hObject, eventdata, handles);


% --- BUTTON SET HISTOGRAM ---
function sethistogram_Callback(hObject, eventdata, handles)
tmp_obj=findobj('Tag','Histogram');
axes(tmp_obj);
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
%set(tmp_obj,'Tag','Histogram');

% --- BUTTON RESET HISTOGRAM ---
function updatehistogram_Callback(hObject, eventdata, handles)
tmp_obj=findobj('Tag','Histogram');
axes(tmp_obj);
[h,n]=tom_hist3d(handles.Histogram);
handles.DataScale=[n(1)  n(size(n,2))];
set(handles.limit_down,'String',handles.DataScale(1));
set(handles.limit_up,'String',handles.DataScale(2));
h=200.*h./(100.*handles.Header.Header.Size(1).*handles.Header.Header.Size(2).*handles.Header.Header.Size(3));
bar(n,h);axis auto;
guidata(hObject,handles);
UPDATE_ALL_Callback(hObject, eventdata, handles);
set(tmp_obj,'Tag','Histogram');

% --- BUTTON SET HISTOGRAM MANUALLY ---
function setmanual_histogramm_Callback(hObject, eventdata, handles)
tmp_obj=findobj('Tag','Histogram');
min=str2num(get(handles.limit_down,'String'));
max=str2num(get(handles.limit_up,'String'));
handles.DataScale=[min max];
guidata(hObject,handles);
set(tmp_obj,'Xlim',[min max]);
UPDATE_ALL_Callback(hObject, eventdata, handles);

% --- EDIT BOX LIMIT MIN ---
function limit_down_Callback(hObject, eventdata, handles)

% --- EDIT BOX LIMIT MAX ---
function limit_up_Callback(hObject, eventdata, handles)

% --- BUTTON ZOOM IN ---
function zoom_in_Callback(hObject, eventdata, handles)
tmp_obj=findobj('Tag','XY_slice');
axes(tmp_obj);
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

% --- BUTTON ZOOM RESET ---
function zoom_reset_Callback(hObject, eventdata, handles)
tmp_obj=findobj('Tag','XY_slice');
axes(tmp_obj);
dim_z=handles.Header.Header.Size(3);
dim_x=handles.Header.Header.Size(1);
dim_y=handles.Header.Header.Size(2);
tmp_obj=findobj('Tag','XY_slider');
slicexy=round(get(tmp_obj,'Value'));
tmp_obj=findobj('Tag','AVG_slider');
sliceavg=round(get(tmp_obj,'Value'));
if slicexy+sliceavg>dim_z; sliceavg=dim_z-slicexy;end;
tmp=double(tom_emreadinc(handles.Filename,[1 1 slicexy],[dim_x-1 dim_y-1 sliceavg-1]));
tmp=mean(tmp,3);
tmp_obj=findobj('Tag','XY_slice');
axes(tmp_obj);
if findstr(handles.Option,'fixed')
    handles.Image=tmp;
    display_real(hObject, eventdata, handles);    
else
    imagesc(tmp',[handles.DataScale]);
end
set(gcf,'DoubleBuffer','on');
set(tmp_obj,'Tag','XY_slice');
handles.actualaxis=[1 dim_x 1 dim_y];
axis([1 dim_x 1 dim_y]);
guidata(hObject,handles);

% --- BUTTON LINE ---
function line_Callback(hObject, eventdata, handles)
set(handles.x1,'String','');set(handles.y1,'String','');
set(handles.x2,'String','');set(handles.y2,'String','');
set(handles.x3,'String','');set(handles.y3,'String','');
set(handles.x4,'String','');set(handles.y4,'String','');
set(handles.radio_line,'Value',1);
set(handles.radio_circle,'Value',0);
set(handles.radio_rectangle,'Value',0);
set(handles.tl1,'String','1st point');
set(handles.tl2,'String','2nd point');
set(handles.tl3,'String','length');
set(handles.tl4,'String','');
message='Click 2 times on the figure to determine the length of the line';
uiwait(msgbox(message));
k = waitforbuttonpress; %waiting for the 1st point
if k==0 %mouse button press
    point1=get(handles.XY_slice,'Currentpoint');
    pt = point1(1,1:2);
    x1=pt(1)*100;y1=pt(2)*100;
    x1=round(x1)/100;y1=round(y1)/100;
    set(handles.x1,'String',x1);
    set(handles.y1,'String',y1);
    h=drawmark(x1,y1);
end
k = waitforbuttonpress;%waiting for the 2nd point
if k==0 %mouse button press
    point1=get(handles.XY_slice,'Currentpoint');
    pt = point1(1,1:2);
    x2=pt(1)*100;y2=pt(2)*100;
    x2=round(x2)/100;y2=round(y2)/100;
    set(handles.x2,'String',x2);
    set(handles.y2,'String',y2);
    %drawmark(x2,y2);
end            
delete(h);            
h=line([x1,x2],[y1,y2],'Color',[1 0 0],'Linewidth',2);
len=sqrt(((x2-x1).^2) + ((y2-y1).^2))*100;
len=round(len)/100;
set(handles.x3,'String',num2str(len));
guidata(hObject,handles);

% --- BUTTON CIRCLE ---
function circle_Callback(hObject, eventdata, handles)
set(handles.x1,'String','');set(handles.y1,'String','');
set(handles.x2,'String','');set(handles.y2,'String','');
set(handles.x3,'String','');set(handles.y3,'String','');
set(handles.x4,'String','');set(handles.y4,'String','');
set(handles.radio_line,'Value',0);
set(handles.radio_circle,'Value',1);
set(handles.radio_rectangle,'Value',0);
set(handles.tl1,'String','center');
set(handles.tl2,'String','radius');
set(handles.tl3,'String','');
set(handles.tl4,'String','');
message='Click 2 times on the figure to determine the center and the radius of a circle';
uiwait(msgbox(message));
k = waitforbuttonpress;%waiting for the 1st point
if k==0 %mouse button press
    point1=get(handles.XY_slice,'Currentpoint');
    pt = point1(1,1:2);
    x1=pt(1)*100;y1=pt(2)*100;
    x1=round(x1)/100;y1=round(y1)/100;
    set(handles.x1,'String',x1);
    set(handles.y1,'String',y1);
    h=drawmark(x1,y1);
end
k = waitforbuttonpress;%waiting for the 2nd point
if k==0 %mouse button press
    point1=get(handles.XY_slice,'Currentpoint');
    pt = point1(1,1:2);
    x2=pt(1)*100;y2=pt(2)*100;
    x2=round(x2)/100;y2=round(y2)/100;
    len=sqrt(((x2-x1).^2) + ((y2-y1).^2))*100;
    len=round(len)/100;
    set(handles.x2,'String',num2str(len));                
end
center= x1 + y1*sqrt(-1);
[u v]=circle(center,len,100);
h=line(u,v,'LineWidth',2,'Color',[1 0 0]);%red dark
guidata(hObject,handles);


% --- BUTTON RECTANGLE ---
function rectangle_Callback(hObject, eventdata, handles)
set(handles.x1,'String','');set(handles.y1,'String','');
set(handles.x2,'String','');set(handles.y2,'String','');
set(handles.x3,'String','');set(handles.y3,'String','');
set(handles.x4,'String','');set(handles.y4,'String','');
set(handles.radio_line,'Value',0);
set(handles.radio_circle,'Value',0);
set(handles.radio_rectangle,'Value',1);
set(handles.tl1,'String','corner 1');
set(handles.tl2,'String','corner 2');
set(handles.tl3,'String','width');
set(handles.tl4,'String','height');
message='Click 2 times on the figure to determine the limit of the square';
uiwait(msgbox(message));
k = waitforbuttonpress; %waiting for the 1st point
if k==0 %mouse button press
    point1=get(handles.XY_slice,'Currentpoint');
    pt = point1(1,1:2);
    x1=pt(1)*100;y1=pt(2)*100;
    x1=round(x1)/100;y1=round(y1)/100;
    set(handles.x1,'String',x1);
    set(handles.y1,'String',y1);
    h=drawmark(x1,y1);
end
k = waitforbuttonpress;%waiting for the 2nd point
if k==0 %mouse button press
    point1=get(handles.XY_slice,'Currentpoint');
    pt = point1(1,1:2);
    x2=pt(1)*100;y2=pt(2)*100;
    x2=round(x2)/100;y2=round(y2)/100;
    set(handles.x2,'String',x2);
    set(handles.y2,'String',y2);               
end
delete(h);
h=rectangle('Position',[x1 y1 (x2-x1) (y2-y1)],'EdgeColor',[1 0 0],'Linewidth',2);
set(handles.x3,'String',(x2-x1));
set(handles.x4,'String',(y2-y1));              
guidata(hObject,handles);

% --- BUTTON COPY TO WORKSPACE ---
function copytoworkspace_Callback(hObject, eventdata, handles)
%filename=[handles.Path '\' handles.Filename];
sss=tom_isemfile(handles.Filename);
if sss
    i=tom_emreadc(handles.Filename);
    %i.Value=double(i.Value);    
else
    i.Value=imread(handles.Filename);
    i.Header=imfinfo(handles.Filename);
end
var1=['vol_' num2str(handles.rnb)];
handles.rnb=handles.rnb+1;
assignin('base',var1,i);
guidata(hObject, handles);


% --- SLIDER Z ---
function XY_slider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject,'SliderStep',[1./16 1./16]);
set(hObject,'Value',[1]);

% --- SLIDER Z. Executes on slider movement ---
function XY_slider_Callback(hObject, eventdata, handles)
dim_z=handles.Header.Header.Size(3);
dim_x=handles.Header.Header.Size(1);
dim_y=handles.Header.Header.Size(2);
slicexy=round(get(hObject,'Value'));
tmp_obj=findobj('Tag','TEXT_XY');
set(tmp_obj,'String',num2str(slicexy));
tmp_obj=findobj('Tag','AVG_slider');
sliceavg=round(get(tmp_obj,'Value'));
rangex_min=handles.actualaxis(1);
rangex=handles.actualaxis(2)-rangex_min;
rangey_min=handles.actualaxis(3);
rangey=handles.actualaxis(4)-rangey_min;
if slicexy+sliceavg>dim_z
    sliceavg=dim_z-slicexy;
    if sliceavg==0
        sliceavg=1;
    end
end
tmp=double(tom_emreadinc(handles.Filename,[rangex_min rangey_min slicexy],[rangex rangey sliceavg-1]));
tmp=mean(tmp,3);
tmp_obj=findobj('Tag','XY_slice');
axes(tmp_obj);
if findstr(handles.Option,'fixed')
    handles.Image=tmp;
    display_real(hObject, eventdata, handles);    
else
    imagesc(tmp',[handles.DataScale]);
end
set(gcf,'DoubleBuffer','on');
set(tmp_obj,'Tag','XY_slice');
guidata(hObject,handles);

% --- EDIT BOX IMAGE NUMBER IN Z ---
function TEXT_XY_Callback(hObject, eventdata, handles)
dim_z=handles.Header.Header.Size(3);
dim_x=handles.Header.Header.Size(1);
dim_y=handles.Header.Header.Size(2);
slicexy=str2num(get(hObject,'String'));
set(handles.XY_slider,'Value',slicexy);
sliceavg=round(get(handles.AVG_slider,'Value'));
rangex_min=handles.actualaxis(1);
rangex=handles.actualaxis(2)-rangex_min;
rangey_min=handles.actualaxis(3);
rangey=handles.actualaxis(4)-rangey_min;
if slicexy+sliceavg>dim_z
    sliceavg=dim_z-slicexy;
    if sliceavg==0
        sliceavg=1;
    end
end
tmp=double(tom_emreadinc(handles.Filename,[rangex_min rangey_min slicexy],[rangex rangey sliceavg-1]));
tmp=mean(tmp,3);
axes(handles.XY_slice);
if findstr(handles.Option,'fixed')
    handles.Image=tmp;
    display_real(hObject, eventdata, handles);    
else
    imagesc(tmp',[handles.DataScale]);
end
set(gcf,'DoubleBuffer','on');
handles.Image=tmp;
set(handles.XY_slice,'Tag','XY_slice');
guidata(hObject,handles);

% --- SLIDER AVERAGE. Reduce noise ---
function AVG_slider_Callback(hObject, eventdata, handles)
sliceavg=round(get(hObject,'Value'));
tmp_obj=findobj('Tag','TEXT_AVG');
set(tmp_obj,'String',num2str(sliceavg));

% --- EDIT BOX AVERAGE NUMBER ---
function TEXT_AVG_Callback(hObject, eventdata, handles)
dim_z=handles.Header.Header.Size(3);
tmp_obj=findobj('Tag','TEXT_AVG');
get(tmp_obj,'Value');
sliceavg=round(eval(get(tmp_obj,'String')));
if sliceavg<1 sliceavg=1; set(tmp_obj,'String',sliceavg); end;
if sliceavg>dim_z sliceavg=dim_z; set(tmp_obj,'String',sliceavg); end;
tmp_obj=findobj('Tag','AVG_slider');
set(tmp_obj,'Value',sliceavg);

%********************************************************************
%*****   Other function  ********************************************
%********************************************************************
function UPDATE_ALL_Callback(hObject, eventdata, handles)
dim_x=handles.Header.Header.Size(1);
dim_y=handles.Header.Header.Size(2);
dim_z=handles.Header.Header.Size(3);
tmp_obj=findobj('Tag','TEXT_XY');
get(tmp_obj,'Value');
slicexy=round(eval(get(tmp_obj,'String')));
tmp_obj=findobj('Tag','AVG_slider');
sliceavg=round(get(tmp_obj,'Value'));
if slicexy+sliceavg>dim_z; sliceavg=dim_z-slicexy;end;
tmp=double(tom_emreadinc(handles.Filename,[1 1 slicexy],[dim_x-1 dim_y-1 sliceavg-1]));
tmp=mean(tmp,3);
tmp_obj=findobj('Tag','XY_slice');
axes(tmp_obj);
if findstr(handles.Option,'fixed')
    handles.Image=tmp;
    display_real(hObject, eventdata, handles);
else
    imagesc(tmp',[handles.DataScale]);
    %t=title(['Info: (' num2str(size(tmp)) ') pixel']); %, mean:' num2str(meanv,3) ' std:' num2str(std,3)
    t=title(['Info: ( ' num2str(dim_x) ' x ' num2str(dim_y) ' x ' num2str(dim_z) ' )']);
end
axis(handles.actualaxis);
colormap gray;
%axis image;
set(tmp_obj,'Tag','XY_slice');
tmp2_obj=findobj('Tag','XY_slider');
slicexy=round(get(tmp2_obj,'Value'));
axes(tmp_obj);

% --- Display the image with 'fixed'option ---
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
        set(gca,'Position',[pp(1) pp(2) sf(2) sf(1)]);
        if isempty(param2)
            t=title(['Info: (' num2str(size(in)) ') pixel']); %, mean:' num2str(meanv,3) ' std:' num2str(std,3)
        end
    otherwise
        axis image; axis ij; %colormap gray; %nothing changed, as nargin=1
        if isempty(param2)
            %t=title(['Info: (' num2str(size(in)) ') pixel, mean:' num2str(meanv,3) ' std:' num2str(std,3) ]);
            t=title(['Info: (' num2str(size(in)) ') pixel']);
        end
end
% ----------- Function setmark_circle -----------
function [X, Y] = circle(w,r,n)
%CIRCLE is used to calculate the coordinate of a circle.
%
%   Syntax: [X,Y]=circle(w,r,n)
%       Input:
%           w: Center of the circle. W must be a complex number as w=x + yi 
%              (x and y are the coordinate)
%           r: Radius of the circle
%           n: The width of the line         
%       Output:
%           X: it is a matrix of coordinate to draw the circle
%           Y: it is a matrix of coordinate to draw the circle

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


% ----------- Function drawmark -----------
function h=drawmark(x,y)
% DRAWMARK draw a circle and a cross to represent the marker 
%   
%   Syntax: drawmark
%       Input:
%           -x: coordinate x
%           -y: coordinate y
%           -i: text to print (must be a string)
%       Output:
%           - 
%   Date: 09/04/03 WDN  

hold on;
Center= x + y*sqrt(-1);
Radius = 1;Gridpt = 100;
[u,v]=circle(Center,Radius,Gridpt);
%line(u,v,'LineWidth',1,'Color',[1 0 0]);%red dark
uu = [x x x x-Radius x+Radius];
vv = [y-Radius y+Radius y y y];
h=line(uu,vv,'LineWidth',1,'color',[1 0 0]);%red dark

hold off;






