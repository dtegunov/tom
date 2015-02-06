function varargout = tom_intervol(varargin)
%TOM_INTERVOL is a tool for interactive visualization of 3D volumes.
%   tom_intervol(volume) is used to view volumes which are in memory, and
%   tom_volxy is used to view volumes which are in a file
%
%   varargout = tom_intervol(varargin)
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
%   a= tom_sphere([64 64 64],4,10,[16 16 1]);
%   tom_intervol(a);
%
%REFERENCES
%
%SEE ALSO
%   TOM_VOLXY, TOM_PARTICLES
%
%   created by FF, WDN 02/20/03
%   updated by WDN 06/28/04
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

error(nargchk(0, 1, nargin, 'struct'))

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_intervol_new_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_intervol_new_OutputFcn, ...
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
%End initialization code - DO NOT EDIT


% --- Executes just before tom_intervol_new is made visible.
function tom_intervol_new_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
%Update handles structure
if isstruct(varargin{1})
    handles.Volume=varargin{1};
    set(handles.volume_name,'String',handles.Volume.Header.Filename);
else
    handles.Volume=tom_emheader(varargin{1});
end
calc_hist(hObject, eventdata, handles);
handles=guidata(handles.figure1);

tmp_obj=findobj('Tag','AVG_slider');
set(tmp_obj,'Min',1);
set(tmp_obj,'Max',handles.Volume.Header.Size(3))
set(tmp_obj,'SliderStep',[1./handles.Volume.Header.Size(3) 1./handles.Volume.Header.Size(3)]);
dim_x=handles.Volume.Header.Size(1);
dim_y=handles.Volume.Header.Size(2);
dim_z=handles.Volume.Header.Size(3);
mdimz=abs(dim_z/2);
%set(handles.volume_name,'String',['Volume   ( ' num2str(dim_x) ' x ' num2str(dim_y) ' x ' num2str(dim_z) ' )']);
%set(handles.volume_name,'String',filename);
set(gcf,'Toolbar','figure');
if dim_x~=dim_y %case of the image has differente height and width
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
tmp_obj=findobj('Tag','XY_slider');
set(tmp_obj,'Min',1);
set(tmp_obj,'Max',handles.Volume.Header.Size(3))
set(tmp_obj,'SliderStep',[1./handles.Volume.Header.Size(3) 5./handles.Volume.Header.Size(3)]);
set(tmp_obj,'Value',handles.dimz);
set(handles.TEXT_XY,'String',num2str(handles.dimz));

handles.actualaxis=[1 dim_x 1 dim_y];
handles.Option_Meas='Line';

handles.actualaxis=[1 dim_x 1 dim_y];
UPDATE_ALL_Callback(hObject, eventdata, handles);
guidata(hObject,handles);

% --- Outputs from this function are returned to the command line.
function varargout = tom_intervol_new_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- BUTTON SET HISTOGRAM ---
function sethistogram_Callback(hObject, eventdata, handles)
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
handles.DataScale=[x(1) x(2)];
set(handles.limit_down,'String',handles.DataScale(1));
set(handles.limit_up,'String',handles.DataScale(2));
guidata(hObject,handles);
set(gca,'Xlim',[x(1) x(2)]);
UPDATE_ALL_Callback(hObject, eventdata, handles);

% --- BUTTON RESET HISTOGRAM ---
function updatehistogram_Callback(hObject, eventdata, handles)
calc_hist(hObject, eventdata, handles);
handles=guidata(handles.figure1);
UPDATE_ALL_Callback(hObject, eventdata, handles);

% --- BUTTON SET HISTOGRAM MANUALLY ---
function setmanual_histogramm_Callback(hObject, eventdata, handles)
min=str2num(get(handles.limit_down,'String'));
max=str2num(get(handles.limit_up,'String'));
handles.DataScale=[min max];
guidata(hObject,handles);
set(handles.Histogram,'Xlim',[min max]);
UPDATE_ALL_Callback(hObject, eventdata, handles);

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

% --- BUTTON ZOOM RESET ---
function zoom_reset_Callback(hObject, eventdata, handles)
axes(handles.XY_slice);
dim_x=handles.Volume.Header.Size(1);
dim_y=handles.Volume.Header.Size(2);
handles.actualaxis=[1 dim_x 1 dim_y];
axis([1 dim_x 1 dim_y]);
guidata(hObject,handles);
UPDATE_ALL_Callback(hObject, eventdata, handles);

% --- BUTTON PIXEL VALUE ---
function pix_Callback(hObject, eventdata, handles)
set(gcf,'Units','pixel');
pixval('on');
guidata(hObject,handles);

% --- BUTTON LINE ---
function line_Callback(hObject, eventdata, handles)
set(handles.x1,'String','');set(handles.y1,'String','');
set(handles.x2,'String','');set(handles.y2,'String','');
set(handles.x3,'String','');set(handles.y3,'String','');
set(handles.x4,'String','');set(handles.y4,'String','');
set(handles.line,'Value',1);
set(handles.circle,'Value',0);
set(handles.rectangle,'Value',0);
set(handles.tl1,'String','1st point');
set(handles.tl2,'String','2nd point');
set(handles.tl3,'String','length');set(handles.y3,'String','pix');
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
set(handles.line,'Value',0);
set(handles.circle,'Value',1);
set(handles.rectangle,'Value',0);
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
set(handles.line,'Value',0);
set(handles.circle,'Value',0);
set(handles.rectangle,'Value',1);
set(handles.tl1,'String','corner 1');
set(handles.tl2,'String','corner 2');
set(handles.tl3,'String','width');set(handles.y3,'String','pix');
set(handles.tl4,'String','height');set(handles.y4,'String','pix');
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

% --- SLIDER Z. Executes on slider movement ---
function XY_slider_Callback(hObject, eventdata, handles)
slicexy=round(get(hObject,'Value'));
set(handles.TEXT_XY,'String',num2str(slicexy));
UPDATE_ALL_Callback(hObject, eventdata, handles)

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

% --- EDIT BOX IMAGE NUMBER IN Z ---
function TEXT_XY_Callback(hObject, eventdata, handles)
editz=round(eval(get(handles.TEXT_XY,'String')));
set(handles.XY_slider,'Value',editz);
UPDATE_ALL_Callback(hObject, eventdata, handles)

% --- SLIDER AVERAGE. Reduce noise ---
function AVG_slider_Callback(hObject, eventdata, handles)
sliceavg=round(get(hObject,'Value'));
set(handles.TEXT_AVG,'String',num2str(sliceavg));

% --- EDIT BOX AVERAGE NUMBER ---
function TEXT_AVG_Callback(hObject, eventdata, handles)
dim_z=handles.Volume.Header.Size(3);
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
dim_x=handles.Volume.Header.Size(1);
dim_y=handles.Volume.Header.Size(2);
dim_z=handles.Volume.Header.Size(3);
%set(handles.volume_name,'String',['Volume   ( ' num2str(dim_x) ' x ' num2str(dim_y) ' x ' num2str(dim_z) ' )']);

slicexy=round(eval(get(handles.TEXT_XY,'String')));
if slicexy>handles.Volume.Header.Size(3)
    slicexy=handles.Volume.Header.Size(3);
    set(handles.XY_slider,'Value',slicexy);
    set(handles.TEXT_XY,'String',slicexy);
elseif slicexy<1
    slicexy=1;
    set(handles.XY_slider,'Value',slicexy);
    set(handles.TEXT_XY,'String',slicexy);
end

sliceavg=round(str2num(get(handles.TEXT_AVG,'String')));
%if slicexy+sliceavg>dim_z;
%    sliceavg=dim_z-slicexy;
%end
rangex_min=handles.actualaxis(1);
rangex=handles.actualaxis(2)-rangex_min;
rangey_min=handles.actualaxis(3);
rangey=handles.actualaxis(4)-rangey_min;
if slicexy+sliceavg>dim_z; 
    sliceavg=dim_z-slicexy;
    if sliceavg==0
        sliceavg=1;
    end    
end;
tmp=handles.Volume.Value(rangex_min:(rangex_min+rangex), rangey_min:(rangey_min+rangey), slicexy:(slicexy+sliceavg-1));
tmp=mean(tmp,3);
axes(handles.XY_slice);
imagesc(tmp',[handles.DataScale]); 
t=title(['Info: ( ' num2str(dim_x) ' x ' num2str(dim_y) ' x ' num2str(dim_z) ' )']);
%axis(handles.actualaxis);
colormap gray;


function calc_hist(hObject, eventdata, handles)
[h,n]=tom_hist3d(handles.Volume.Value(:,:,1));
handles.DataScale=[n(1)  n(size(n,2))];
h=200.*h./(100.*handles.Volume.Header.Size(1).*handles.Volume.Header.Size(2).*handles.Volume.Header.Size(3));
axes(handles.Histogram);bar(n,h);axis auto;
set(handles.limit_down,'String',handles.DataScale(1));
set(handles.limit_up,'String',handles.DataScale(2));
guidata(handles.figure1,handles);

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




