function varargout = tom_ctffit(varargin)
%TOM_CTFFIT creates ...
%
%   varargout = tom_ctffit(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout   		...
%
%  EXAMPLE
%   tom_ctffit
%   or
%   XXX=tom_emreadc('C:\MATLAB6p5\work\bild_1.em');
%   tom_ctffit(XXX,-4500);
%
%REFERENCES
%
%SEE ALSO
%   TOM_EMREAD, TOM_EMWRITE, TOM_CTFZERO, TOM_PS, IMAGESC
%
%   created by William Del Net 04/30/02
%   updated by William Del Net 07/11/03
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

%error(nargchk(0, 1, nargin, 'struct'))

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_ctffit_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_ctffit_OutputFcn, ...
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

% --- Executes just before tom_ctffit is made visible.
function tom_ctffit_OpeningFcn(hObject, eventdata, handles, varargin)
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_ctffit (see VARARGIN)
handles.output = hObject;
if nargin>4
    im=varargin{1};
    defocus=varargin{2}*10^-3;
elseif nargin>3
    im=varargin{1};
    defocus=im.Header.Defocus*0.1*10^-3;%units in nm
end
if nargin>3
    set(handles.name_image,'String','Name: ?');
    
    if (isstruct(im)==0)
        im=tom_emheader(im);
    end
    im.Value=double(im.Value);

    axes(handles.axes1);
    tom_imagesc(im.Value,'noinfo');%tom_setmark_imagesc(im.Value);
    drawnow;
    if isfield(im,'Header')==1 %case of existing Header
        if im.Header.Voltage==0 %case of voltage 0
            uiwait(tom_crea_header(handles, im));
            qw=guidata(handles.load);
            im.Header=qw.Header;
        end
    else 
        %case of non existing Header
    end
    if ~isequal(varargin{3},'ps')
        ps = tom_ps(im.Value);
        ps =log(ps);
    else
        ps=im.Value;
    end;

    ps(size(ps,1)/2+1,size(ps,1)/2+1)=0; %set middle to zero
    
    axes(handles.axes2);
    %tom_setmark_imagesc(ps);
    [mean, max, min, std, variance] = tom_dev(ps,'noinfo');
    imagesc(ps',[(1*mean)-(1*std) (1*mean)+(1*std)]);%colormap gray;
    slider_max=get(handles.slider1,'Max');slider_min=get(handles.slider1,'Min');
    slider_step(1)=1000/(slider_max-slider_min);slider_step(2)=0.1; 
    set(handles.slider1,'Sliderstep',slider_step)
    set(handles.slider1,'Value',defocus*10^3);
    if defocus<0.001 & defocus>-0.001
        defocus=0;
    end
    set(handles.Defocus,'String',defocus*10^3);
    handles.Deviation.Mean=mean;
    handles.Deviation.Max=max;
    handles.Deviation.Min=min;
    handles.Deviation.Std=std;
    handles.Deviation.Variance=variance;
    handles.Deviation.Cont=1;
    handles.Deviation.Bright=1;
    handles.NameImage=[];
    handles.Image=im;
    handles.PowerSpektrum=ps;
    if im.Header.Defocus~=0
        defocus=im.Header.Defocus*0.1*10^-3;
        set(handles.slider1,'Value',defocus*10^3);        
        set(handles.Defocus,'String',defocus*10^3);
        disp_circle(defocus,handles);
    end

    %disp_circle(defocus,handles);
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_ctffit wait for user response (see UIRESUME)
% uiwait(handles.ctffit);


% --- Outputs from this function are returned to the command line.
function varargout = tom_ctffit_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%*************************************************************************%
% --- LOAD IMAGE ---
function load_Callback(hObject, eventdata, handles)
set(handles.Defocus,'String','');
set(handles.step_10,'Value',0);
set(handles.step_100,'Value',0);
set(handles.step_1000,'Value',1);
set(handles.step_10000,'Value',0);
[myfile, mypathname] = uigetfile('*.*', 'LOAD IMAGE (EM format)');
myfile=num2str(myfile);
mypathname=num2str(mypathname);
if (strcmp(myfile,'0') & strcmp(mypathname,'0'))
    %no action. User has pressed cancel button
else
    im=tom_emreadc([mypathname myfile]);
    set(handles.name_image,'String',['Name: ' myfile]);
    im.Value=double(im.Value);
    axes(handles.axes1);
    %tom_imagesc(im.Value,'noinfo');
    tom_setmark_imagesc(im.Value);%imagesc(im.Value);colormap gray;%
    drawnow;
    if isfield(im,'Header')==1 %case of existing Header
        if im.Header.Voltage==0 %case of voltage 0
            uiwait(tom_crea_header(handles, im));
            qw=guidata(handles.load);
            im.Header=qw.Header;
        end
    else %case of non existing Header
    end
    ps = tom_ps(im.Value);
    ps(size(ps,1)/2+1,size(ps,1)/2+1)=0; %set mean to zero
    axes(handles.axes2);
    [mean, max, min, std, variance] = tom_dev(ps,'noinfo');
    imagesc(ps',[(1*mean)-(1*std) (1*mean)+(1*std)]);%colormap gray;
    %imagesc(ps,[(mean-3) (mean+3)]);
    %tom_setmark_imagesc(ps);
    slider_max=get(handles.slider1,'Max');slider_min=get(handles.slider1,'Min');
    slider_step(1)=1000/(slider_max-slider_min);slider_step(2)=0.1; 
    set(handles.slider1,'Sliderstep',slider_step,'Value',0);
    set(handles.Defocus,'String',0);
    handles.Deviation.Mean=mean;
    handles.Deviation.Max=max;
    handles.Deviation.Min=min;
    handles.Deviation.Std=std;
    handles.Deviation.Variance=variance;
    handles.Deviation.Cont=1;
    handles.Deviation.Bright=1;
    handles.NameImage=[mypathname myfile];
    handles.Image=im;  
    handles.PowerSpektrum=ps;
    guidata(hObject, handles);
    if im.Header.Defocus~=0
        defocus=im.Header.Defocus*0.1*10^-3;
        set(handles.slider1,'Value',defocus*10^3);        
        set(handles.Defocus,'String',defocus*10^3);
        disp_circle(defocus,handles);
    end
    guidata(hObject, handles);
end
% --- EXIT TOM_CTFFIT ---
function exit_Callback(hObject, eventdata, handles)
close(handles.ctffit);
% --- Executes on button press in semicircle.
function semicircle_Callback(hObject, eventdata, handles)
set(handles.semicircle,'Value',1);
set(handles.circle,'Value',0);
defocus=str2num(get(handles.Defocus,'String'))*10^-3;
set(handles.slider1,'Value',defocus*10^3);
axes(handles.axes2);
K=handles.Deviation.Cont;
B=handles.Deviation.Bright;
imagesc(handles.PowerSpektrum',[(B*handles.Deviation.Mean)-(K*handles.Deviation.Std) (B*handles.Deviation.Mean)+(K*handles.Deviation.Std)]);
%imagesc(handles.PowerSpektrum',[handles.Deviation.Mean-(handles.Deviation.Cont*handles.Deviation.Std) handles.Deviation.Mean+(handles.Deviation.Cont*handles.Deviation.Std)]);
if defocus<0.001 & defocus>-0.001
    defocus=0;
end
disp_circle(defocus,handles);

% --- Executes on button press in circle.
function circle_Callback(hObject, eventdata, handles)
set(handles.semicircle,'Value',0);
set(handles.circle,'Value',1);
defocus=str2num(get(handles.Defocus,'String'))*10^-3;
set(handles.slider1,'Value',defocus*10^3);
axes(handles.axes2);
K=handles.Deviation.Cont;
B=handles.Deviation.Bright;
imagesc(handles.PowerSpektrum',[(B*handles.Deviation.Mean)-(K*handles.Deviation.Std) (B*handles.Deviation.Mean)+(K*handles.Deviation.Std)]);
%imagesc(handles.PowerSpektrum',[handles.Deviation.Mean-(handles.Deviation.Cont*handles.Deviation.Std) handles.Deviation.Mean+(handles.Deviation.Cont*handles.Deviation.Std)]);
if defocus<0.001 & defocus>-0.001
    defocus=0;
end
disp_circle(defocus,handles);

% --- DEFOCUS RANGE 10nm. ---
function step_10_Callback(hObject, eventdata, handles)
set(handles.step_10,'Value',1);
set(handles.step_100,'Value',0);
set(handles.step_1000,'Value',0);
set(handles.step_10000,'Value',0);
slider_max=get(handles.slider1,'Max');slider_min=get(handles.slider1,'Min');
slider_step(1)=10/(slider_max-slider_min);slider_step(2)=0.1; 
set(handles.slider1,'Sliderstep',slider_step);

% --- DEFOCUS RANGE 100nm. ---
function step_100_Callback(hObject, eventdata, handles)
set(handles.step_10,'Value',0);
set(handles.step_100,'Value',1);
set(handles.step_1000,'Value',0);
set(handles.step_10000,'Value',0);
slider_max=get(handles.slider1,'Max');slider_min=get(handles.slider1,'Min');
slider_step(1)=100/(slider_max-slider_min);slider_step(2)=0.1; 
set(handles.slider1,'Sliderstep',slider_step);

% --- DEFOCUS RANGE 1000nm. ---
function step_1000_Callback(hObject, eventdata, handles)
set(handles.step_10,'Value',0);
set(handles.step_100,'Value',0);
set(handles.step_1000,'Value',1);
set(handles.step_10000,'Value',0);
slider_max=get(handles.slider1,'Max');slider_min=get(handles.slider1,'Min');
slider_step(1)=1000/(slider_max-slider_min);slider_step(2)=0.1; 
set(handles.slider1,'Sliderstep',slider_step);

% --- DEFOCUS RANGE 10000nm. ---
function step_10000_Callback(hObject, eventdata, handles)
set(handles.step_10,'Value',0);
set(handles.step_100,'Value',0);
set(handles.step_1000,'Value',0);
set(handles.step_10000,'Value',1);
slider_max=get(handles.slider1,'Max');slider_min=get(handles.slider1,'Min');
slider_step(1)=10000/(slider_max-slider_min);slider_step(2)=0.1; 
set(handles.slider1,'Sliderstep',slider_step);

% --- DEFOCUS SLIDER ---
function slider1_Callback(hObject, eventdata, handles)
axes(handles.axes2);
K=handles.Deviation.Cont;
B=handles.Deviation.Bright;
imagesc(handles.PowerSpektrum',[(B*handles.Deviation.Mean)-(K*handles.Deviation.Std) (B*handles.Deviation.Mean)+(K*handles.Deviation.Std)]);
defocus=get(hObject,'Value')*10^-3;
if defocus<0.001 & defocus>-0.001
    defocus=0;
end
set(handles.Defocus,'String',defocus*10^3);
disp_circle(defocus,handles);

% --- DEFOCUS VALUE ---
function Defocus_Callback(hObject, eventdata, handles)
defocus=str2num(get(handles.Defocus,'String'))*10^-3;
set(handles.slider1,'Value',defocus*10^3);
axes(handles.axes2);
K=handles.Deviation.Cont;
B=handles.Deviation.Bright;
imagesc(handles.PowerSpektrum',[(B*handles.Deviation.Mean)-(K*handles.Deviation.Std) (B*handles.Deviation.Mean)+(K*handles.Deviation.Std)]);
if defocus<0.001 & defocus>-0.001
    defocus=0;
end
disp_circle(defocus,handles);

% --- SAVE DEFOCUS VALUE ---
function Save_Callback(hObject, eventdata, handles)
defocus=str2num(get(handles.Defocus,'String'))*10; %save in Angstroem
temp=handles.Image;
temp.Header.Defocus=defocus;
if isempty(handles.NameImage)
    [myfile, mypathname] = uiputfile('*.*', 'SAVE DEFOCUS AS');
    myfile=num2str(myfile);
    mypathname=num2str(mypathname);
    if (strcmp(myfile,'0') & strcmp(mypathname,'0'))
        %no action. User has pressed cancel button
    else
        handles.NameImage=[mypathname myfile];
        set(handles.name_image,'String',['Name: ' myfile]);
        guidata(hObject, handles);   
        tom_emwrite(handles.NameImage,temp);                             
    end
else
    tom_emwrite(handles.NameImage,temp);
end

% --- CLEAR THE CIRCLE ---
function clear_Callback(hObject, eventdata, handles)
set(handles.Defocus,'String',0);
set(handles.slider1,'Value',(0*10^3));
axes(handles.axes2);
imagesc(handles.PowerSpektrum',[handles.Deviation.Mean-handles.Deviation.Std handles.Deviation.Mean+handles.Deviation.Std]);
handles.Deviation.Cont=1;
handles.Deviation.Bright=1;
guidata(hObject, handles); 

% --- CONTRAST UP ---
function Contrast_up_Callback(hObject, eventdata, handles)
axes(handles.axes2);
K=(handles.Deviation.Cont-0.02);
B=handles.Deviation.Bright;
if K<0.00001
    K=0.00001;
end
imagesc(handles.PowerSpektrum',[(B*handles.Deviation.Mean)-(K*handles.Deviation.Std) (B*handles.Deviation.Mean)+(K*handles.Deviation.Std)]);
handles.Deviation.Cont=K;
guidata(hObject, handles); 
Defocus_Callback(hObject, eventdata, handles);

% --- CONTRAST DOWN ---
function Contrast_down_Callback(hObject, eventdata, handles)
axes(handles.axes2);
K=(handles.Deviation.Cont+0.02);
B=handles.Deviation.Bright;
if K>2
    K=2;
end
imagesc(handles.PowerSpektrum',[(B*handles.Deviation.Mean)-(K*handles.Deviation.Std) (B*handles.Deviation.Mean)+(K*handles.Deviation.Std)]);
handles.Deviation.Cont=K;
guidata(hObject, handles);  
Defocus_Callback(hObject, eventdata, handles);

% --- BRIGHTNESS UP ---
function Bright_up_Callback(hObject, eventdata, handles)
axes(handles.axes2);
B=(handles.Deviation.Bright-0.05);
K=handles.Deviation.Cont;
if B<0.00001
    B=0.00001;
end
imagesc(handles.PowerSpektrum',[(B*handles.Deviation.Mean)-(K*handles.Deviation.Std) (B*handles.Deviation.Mean)+(K*handles.Deviation.Std)]);
handles.Deviation.Bright=B;
guidata(hObject, handles); 
Defocus_Callback(hObject, eventdata, handles);

% --- BRIGHTNESS DOWN ---
function Bright_down_Callback(hObject, eventdata, handles)
axes(handles.axes2);
B=(handles.Deviation.Bright+0.05);
K=handles.Deviation.Cont;
if B>3
    B=3;
end
imagesc(handles.PowerSpektrum',[(B*handles.Deviation.Mean)-(K*handles.Deviation.Std) (B*handles.Deviation.Mean)+(K*handles.Deviation.Std)]);
handles.Deviation.Bright=B;
guidata(hObject, handles); 
Defocus_Callback(hObject, eventdata, handles);

function disp_circle(defocus,handles)
Pixelsize=handles.Image.Header.CCDArea./handles.Image.Header.Size(1);
if handles.Image.Header.Postmagnification==0
    handles.Image.Header.Postmagnification=1;
end
pix_size= Pixelsize/(handles.Image.Header.Magnification*handles.Image.Header.Postmagnification)*10^3;
voltage = handles.Image.Header.Voltage/1000;
Cs=handles.Image.Header.Cs;
out=tom_ctfzero(defocus,pix_size, voltage,size(handles.Image.Value,1),Cs);
ss=size(out);
axes(handles.axes2);
sxz=size(handles.Image.Value);
hold on;
Center= ((sxz(1)/2)+1) + ((sxz(2)/2)+1)*sqrt(-1);
Gridpt = 200;
for i=1:ss(2)       
    Radius = out(i);
    if get(handles.semicircle,'Value')==1
        %[u,v]=semi_circle(Center,Radius,Gridpt);
        [u,v]=circle(Center,Radius,Gridpt,'semi-circle');
    elseif get(handles.circle,'Value')==1
        %[u,v]=circle(Center,Radius,Gridpt);
        [u,v]=circle(Center,Radius,Gridpt,'circle');
    end
    line(u,v,'LineWidth',1,'Color',[1 0 0]);%red dark
end
hold off;

function [X, Y] = circle(w,r,n,What)
%   This function is used to calculate a matrix of a 1/2 circle (between pi/2 and -pi/2) or a circle.
%   Syntax: [X,Y]=circle(w,r,n,What)
%       Input:
%           w: Center of the 1/2 circle. W must be a complex number as w=x + yi 
%              (x and y are the coordinate)
%           r: Radius of the circle
%           n: nb of point. 
%           What: case 'semi-circle': matrix of a 1/2 circle (between pi/2 and -pi/2)
%                 case 'circle':matrix of a circle
%       Output:
%           X: it is a matrix of coordinate to draw the circle
%           Y: it is a matrix of coordinate to draw the circle
%
w1 = real(w);
w2 = imag(w);
for k = 1:n
    switch What
        case 'semi-circle'
            t = -k*(pi)/n;
            X(k) = w1 + r*sin(t);
            Y(k) = w2 + r*cos(t);
        case 'circle'                
            t = k*pi/n;
            X(k) = w1 + r*cos(t);
            Y(k) = w2 + r*sin(t);
            X(n+k) = w1 - r*cos(t);
            Y(n+k) = w2 - r*sin(t);
            X(2*n+1) = X(1);
            Y(2*n+1) = Y(1);
    end
end


function tom_setmark_imagesc(in,parameter)
% TOM_SETMARK_IMAGESC
%    Syntaxe: tom_setmark_imagesc(in);
%       Input:
%           in: Array of data of the image
%           parameter: 'fixed' to have one pixel on screen equal to one
%                       pixel on file
%       Output:
%           - 
%   Date: 12/08/02 SN   
%   Last modified: 25/04/03 WDN

in_red=imresize(in,.1);
[meanv max min std]=tom_devinternal(in_red);
if (meanv-4*std)>=(meanv+4*std)
    imagesc(in');
else
    imagesc(in',[meanv-4*std meanv+4*std]);colormap gray;axis image;
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
%TOM_DEVINTERNAL is used only by TOM_SETMARKSC
%   Date: 12/08/02 SN   

[s1,s2,s3]=size(A);
a=sum(sum(sum(A)))/(s1*s2*s3);
b=max(max(max(A)));
c=min(min(min(A)));
d=std2(A);
e=d^2;

