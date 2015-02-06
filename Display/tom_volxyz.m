function varargout = tom_volxyz(varargin)
%TOM_VOLXYZ is a 3D visualization tool for tomograms
%
%   varargout = tom_volxyz(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout   		...
%
%   This function is used to view a volume inp_vol directly from the memory. 
%   It is a GUI to scan through the volume in x,y and z direction and lets you adjust 
%   the contrast interactively. Additionally, a running average in every direction
%   can be calculated to increase the contrast. A bandpass filter can be applied
%   in each direction.
%
%EXAMPLE
%   a=tom_emread('pyrodictium.vol');
%   tom_volxyz(a.Value);
%
%REFERENCES
%
%SEE ALSO
%   TOM_INTERVOL, TOM_VOLXY, TOM_PARTICLES, TOM_DSPCUB
%
%   created by AK 03/20/05
%   completely rewritten AK 09/23/05
%   updated by AK 10/07/05
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_volxyz_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_volxyz_OutputFcn, ...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Opening Function                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tom_volxyz_OpeningFcn(hObject, eventdata, handles, varargin)

set(0,'RecursionLimit',1000);

%due to problems accessing and updating the handles structure in some
%callbacks this global is used for application data throughout this code.
global storage;
storage = struct();
setappdata(0,'UseNativeSystemDialogs',0);
disp('Loading data, please wait...');

%read in volume from file if no argument was passed to us
if nargin == 3
    varargin{1} = tom_emread;
end

%initialise volume and add header if needed
if isstruct(varargin{1})
    storage.Volume=varargin{1};
    set(handles.filename,'String',storage.Volume.Header.Filename);
else
    storage.Volume=tom_emheader(varargin{1});
end

%particle picking structure (only available if av3_new directory is in path)
try
    if ~isstruct(varargin{2})
        storage.Align =  tom_av3_createnewalignstruct();
    else
        storage.Align = varargin{2};
    end
    storage.particleradius = 10;
end



%initialise some storage variables
storage.dimensions.x = storage.Volume.Header.Size(1);
storage.dimensions.y = storage.Volume.Header.Size(2);
try
    storage.dimensions.z = storage.Volume.Header.Size(3);
catch
    storage.dimensions.z = 1;
end
storage.figures.vline.xy = [];
storage.figures.hline.xy = [];
storage.figures.hline.xz = [];
storage.figures.hline.yz = [];
storage.actualaxis.xy = [1 storage.dimensions.x 1 storage.dimensions.y];
storage.actualaxis.xz = [1 storage.dimensions.x 1 storage.dimensions.z];
storage.actualaxis.yz = [1 storage.dimensions.z 1 storage.dimensions.y];
storage.subvol.data = [];
storage.slices.xy = [];
storage.slices.xz = [];
storage.slices.yz = [];
storage.rotate.phi = 0;
storage.rotate.psi = 0;
storage.rotate.theta = 0;
storage.profile.plane = 'xy_axes';

storage.profile.valcell = {};
storage.profile.coordcell = {};

%measurement variables
storage.measurement.point1 = [];
storage.measurement.point2 = [];
storage.measurement.point3 = [];


%These values are the borders in the volume figure window
storage.display.outer_padding = 50;
storage.display.inner_padding = 20;
storage.display.outerdiff = 20;
storage.display.slicefontsize = 0;
storage.display.ticks = 1;
storage.display.scalebar = 0;
storage.display.pixelinfo = 0;
storage.display.slicenumbers = 0;
storage.display.crosshair = 1;
storage.display.box = 0;

%initialise bandpass filter values
set(handles.bandpass_low,'String', num2str(0));
set(handles.bandpass_high,'String',num2str(round(storage.dimensions.x ./ 2)));
storage.filter.low = 0;
storage.filter.high = round(storage.dimensions.x ./ 2);
storage.filter.xy = 0;
storage.filter.xz = 0;
storage.filter.yz = 0;
storage.mean = 0;

if storage.dimensions.z == 1
    vol = zeros(storage.dimensions.x, storage.dimensions.y, 2);
    vol(:,:,1) = storage.Volume.Value;
    vol(:,:,2) = storage.Volume.Value;
    storage.Volume.Value = vol;
    storage.dimensions.z = 2;
    storage.actualaxis.xz = [1 storage.dimensions.x 1 storage.dimensions.z];
    storage.actualaxis.yz = [1 storage.dimensions.z 1 storage.dimensions.y];
end


%Create image slice window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
handles.figures.image = figure();

if ~strcmp(computer,'PCWIN')
    set(handles.figures.image,'CloseRequestFcn',@(h,e)closefcn(hObject, eventdata, handles));
else
    set(handles.figures.image,'CloseRequestFcn',@dummy);
end
%set(handles.figures.image,'CloseRequestFcn',@dummy);
set(handles.figures.image,'Tag','tom_volxyz_imagewindow');
set(handles.figures.image,'Name','tom_volxyz','Color',[1 1 1]);
set(handles.figures.image,'Resize','on');
set(handles.figures.image,'Renderer','zbuffer');
set(handles.figures.image,'DoubleBuffer','on');
figurewidth = storage.display.outer_padding * 2 + storage.display.inner_padding + storage.dimensions.z + storage.dimensions.x;
figureheight = storage.display.outer_padding * 2 + storage.display.inner_padding + storage.dimensions.y + storage.dimensions.z;
storage.figureaspect = figureheight ./ figurewidth;

set(handles.figures.image,'Position',[100 100 figurewidth figureheight]);
set(handles.figures.image,'ToolBar','none');
set(handles.figures.image,'MenuBar','none');
set(handles.figures.image,'Resize','on');
set(handles.figures.image,'Units','Pixel');


%Position of xy Slice
handles.axes.xy = axes();
axis ij;
left = storage.display.outer_padding;
bottom = storage.display.outer_padding + storage.dimensions.z + storage.display.inner_padding;
width = storage.dimensions.x;
height = storage.dimensions.y;
set(handles.axes.xy,'Units','pixel');
set(handles.axes.xy,'Tag','xy_axes');
set(handles.axes.xy, 'Position', [left bottom width height],'visible','on');
set(handles.axes.xy,'Units','normalized');

%Position of xz slice
handles.axes.xz = axes();
axis ij;
left = storage.display.outer_padding;
bottom = storage.display.outer_padding;
width = storage.dimensions.x;
height = storage.dimensions.z;
set(handles.axes.xz,'Units','pixel');
set(handles.axes.xz, 'Tag', 'xz_axes');
set(handles.axes.xz, 'Position', [left bottom width height],'visible','on');
set(handles.axes.xz,'Units','normalized');

%Position of yz slice
handles.axes.yz = axes();
axis ij;
left = storage.display.outer_padding + storage.display.inner_padding + storage.dimensions.x;
bottom = storage.display.outer_padding + storage.dimensions.z + storage.display.inner_padding;
width = storage.dimensions.z;
height = storage.dimensions.y;
set(handles.axes.yz,'Units','pixel');
set(handles.axes.yz,'Tag','yz_axes');
set(handles.axes.yz, 'Position', [left bottom width height],'visible','on');
set(handles.axes.yz,'Units','normalized');
set(handles.axes.yz,'YAxisLocation','right','Color','none');

%Position of 3D box
handles.axes.box = axes();
axis ij;axis off;
left =storage.display. outer_padding + storage.display.inner_padding + storage.dimensions.x;
bottom = storage.display.outer_padding + 1;
width = storage.dimensions.z;
height = storage.dimensions.z;
set(handles.axes.box,'Units','pixel');
set(handles.axes.box,'Tag','box');
set(handles.axes.box,'Position',[left bottom width height],'visible','on');
set(handles.axes.box,'Units','normalized','Visible','off');

% center the view
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
center_position();


% adjust average sliders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(handles.averageslider_x,'Min',1);
set(handles.averageslider_x,'Max',storage.dimensions.x)
set(handles.averageslider_x, 'Value', 1);
set(handles.averageslider_x,'SliderStep',[1 ./ storage.dimensions.x, 1 ./ storage.dimensions.x]);
storage.average.x = 1;

set(handles.averageslider_y,'Min',1);
set(handles.averageslider_y,'Max',storage.dimensions.y)
set(handles.averageslider_y, 'Value', 1);
set(handles.averageslider_y,'SliderStep',[1 ./ storage.dimensions.y, 1 ./ storage.dimensions.y]);
storage.average.y = 1;


set(handles.averageslider_z,'Min',1);
set(handles.averageslider_z,'Max',storage.dimensions.z)
set(handles.averageslider_z, 'Value', 1);
set(handles.averageslider_z,'SliderStep',[1 ./ storage.dimensions.z, 1  ./ storage.dimensions.z]);
storage.average.z = 1;


%calculate histogram of the volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calculate_histogram();
render_histogram();
storage.DataScale_orig = storage.DataScale;

%render the slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
redraw_all_slices();

set(handles.figures.image,'ResizeFcn',@(h,e)(resize_window(hObject,eventdata,handles)));

%extract object pixel size from volume header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    if storage.Volume.Header.Objectpixelsize ~= 0
        storage.objpixelsize = storage.Volume.Header.Objectpixelsize;
        set(findobj('Tag','scalebar_factor'),'String',num2str(storage.objpixelsize .* 10));
    else
        storage.objpixelsize = 0;
    end
catch
    storage.objpixelsize = 0;
end;

%Adjust the gui size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set(0,'Units','pixel');
%screensize = get(0,'ScreenSize');
movegui(findobj('Tag','tom_volxyz'),'east');


handles.output = hObject;
guidata(hObject, handles);

disp('Finished.');

try;
    if isstruct(varargin{2})
        uiwait();
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Close Function                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function closefcn(hObject, eventdata, handles)

global storage;
uiresume;

tmpobj = findobj('Tag','tom_volxyz');
if tmpobj ~= 0
    delete(tmpobj);
end

tmpobj = findobj('Tag','tom_volxyz_imagewindow');
if tmpobj ~= 0
    delete(tmpobj);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_volxyz_OutputFcn(hObject, eventdata, handles) 

global storage;
try
    varargout{1} = storage.Align(1,2:end);
end
%varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Set Histogram                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_sethistogram_Callback(hObject, eventdata, handles)

global storage;

axes(findobj('Tag','histogram'));
waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
finalRect = rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions
x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
storage.DataScale=[x(1) x(2)];
set(findobj('Tag','histogram_low'),'String',storage.DataScale(1));
set(findobj('Tag','histogram_high'),'String',storage.DataScale(2));
guidata(hObject,handles);
set(gca,'Xlim',[x(1) x(2)]);

redraw_all_slices();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Set Histogram Manually                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_setmanually_Callback(hObject, eventdata, handles)

global storage;

min=str2num(get(findobj('Tag','histogram_low'),'String'));
max=str2num(get(findobj('Tag','histogram_high'),'String'));

if max>min
    storage.DataScale = [min max];
    set(findobj('Tag','histogram'),'Xlim',[min max]);
else
    set(findobj('Tag','histogram_low'),'String',num2str(storage.DataScale(1)));
    set(findobj('Tag','histogram_high'),'String',num2str(storage.DataScale(2)));
end

redraw_all_slices();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Reset Histogram                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_resethistogram_Callback(hObject, eventdata, handles)

global storage;

storage.DataScale = storage.DataScale_orig;

set(findobj('Tag','histogram'),'Xlim',storage.DataScale);

redraw_all_slices();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Pick lower limit for Histogram                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_histogram_picklow_Callback(hObject, eventdata, handles)

global storage;

figure(findobj('Tag','tom_volxyz_imagewindow'));
[x y] = ginput(1);

if x > 0 & y > 0 & strcmp(get(gca,'Tag'),'box') == 0 & ~isempty(get(gca,'Tag'))
    if strcmp(get(gca,'Tag'),'xy_axes') == 1
        val = storage.Volume.Value(round(x),round(y),storage.position.xy);
    elseif strcmp(get(gca,'Tag'),'xz_axes') == 1
        val = storage.Volume.Value(round(x),storage.position.xz,round(y));
    elseif strcmp(get(gca,'Tag'),'yz_axes') == 1
        val = storage.Volume.Value(storage.position.yz,round(y),round(x));
    end
    
    if val >= storage.DataScale(2)
        errordlg('lower limit must be smaller than upper limit!','Histrogram error');
    else
        storage.DataScale = [val storage.DataScale(2)];
        set(findobj('Tag','histogram_low'),'String',num2str(val));
        redraw_all_slices();
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Pick upper limit for Histogram                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_histogram_pickhigh_Callback(hObject, eventdata, handles)

global storage;

figure(findobj('Tag','tom_volxyz_imagewindow'));
[x y] = ginput(1);

if x > 0 & y > 0 & strcmp(get(gca,'Tag'),'box') == 0 & ~isempty(get(gca,'Tag'))
    if strcmp(get(gca,'Tag'),'xy_axes') == 1
        val = storage.Volume.Value(round(x),round(y),storage.position.xy);
    elseif strcmp(get(gca,'Tag'),'xz_axes') == 1
        val = storage.Volume.Value(round(x),storage.position.xz,round(y));
    elseif strcmp(get(gca,'Tag'),'yz_axes') == 1
        val = storage.Volume.Value(storage.position.yz,round(y),round(x));
    end

    if val <= storage.DataScale(1)
        errordlg('upper limit must be greater than lower limit!','Histrogram error');
    else
        storage.DataScale = [storage.DataScale(1) val];
        set(findobj('Tag','histogram_high'),'String',num2str(val));
        redraw_all_slices();
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Zoom                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_zoom_Callback(hObject, eventdata, handles)

global storage;

unset_callbacks();

waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
finalRect = rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions

offset(1) = offset(2);

xdim = offset(1);
ydim = offset(1);
zdim = offset(1);

orientation = get(gco,'Tag');

if strcmp(orientation,'xy_image') == 1
	x = [p1(1) p1(1)+offset(1)];
	y = [p1(2) p1(2)+offset(2)];
	x=round(x);
	y=round(y);
	
	storage.actualaxis.xy=[x(1) x(2) y(1) y(2)];
	storage.actualaxis.xz=[x(1) x(2) storage.position.xy-zdim./2 storage.position.xy+zdim./2];
	storage.actualaxis.yz=[storage.position.xy-zdim./2 storage.position.xy+zdim./2 y(1) y(2)];
elseif strcmp(orientation,'xz_image') == 1
	x = [p1(1) p1(1)+offset(1)];
	z = [p1(2) p1(2)+offset(2)];
	x=round(x);
	z=round(z);
	
	storage.actualaxis.xy=[x(1) x(2) storage.position.xy-ydim./2 storage.position.xy+ydim./2];
	storage.actualaxis.xz=[x(1) x(2) z(1) z(2)];
	storage.actualaxis.yz=[z(1) z(2) storage.position.xy-ydim./2 storage.position.xy+ydim./2];
elseif strcmp(orientation, 'yz_image') == 1
	z = [p1(1) p1(1)+offset(1)];
	y = [p1(2) p1(2)+offset(2)];
	z=round(z);
	y=round(y);
	
	storage.actualaxis.xy=[storage.position.yz-xdim./2 storage.position.yz+xdim./2 y(1) y(2)];
	storage.actualaxis.xz=[storage.position.yz-xdim./2 storage.position.yz+xdim./2 z(1) z(2)];
	storage.actualaxis.yz=[z(1) z(2) y(1) y(2)];
end

set_callbacks();
redraw_all_slices();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Zoom Reset                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_zoomreset_Callback(hObject, eventdata, handles)

global storage;

storage.actualaxis.xy = [1 storage.dimensions.x 1 storage.dimensions.y];
storage.actualaxis.xz = [1 storage.dimensions.x 1 storage.dimensions.z];
storage.actualaxis.yz = [1 storage.dimensions.z 1 storage.dimensions.y];

redraw_all_slices();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Bandpass xy checkbox                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bandpass_xy_Callback(hObject, eventdata, handles)

global storage;

storage.filter.xy = get(hObject,'Value');
storage.filter.low = str2num(get(findobj('Tag','bandpass_low'),'String'));
storage.filter.high = str2num(get(findobj('Tag','bandpass_high'),'String'));

render_slice('xy');
render_measuremarks();

if storage.display.pixelinfo == 1
    impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Bandpass xz checkbox                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bandpass_xz_Callback(hObject, eventdata, handles)

global storage;

storage.filter.xz = get(hObject,'Value');
storage.filter.low = str2num(get(findobj('Tag','bandpass_low'),'String'));
storage.filter.high = str2num(get(findobj('Tag','bandpass_high'),'String'));

render_slice('xz');
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Bandpass yz checkbox                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bandpass_yz_Callback(hObject, eventdata, handles)

global storage;

storage.filter.yz = get(hObject,'Value');
storage.filter.low = str2num(get(findobj('Tag','bandpass_low'),'String'));
storage.filter.high = str2num(get(findobj('Tag','bandpass_high'),'String'));

render_slice('yz');
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Bandpass filter low                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bandpass_low_Callback(hObject, eventdata, handles)

global storage;

storage.filter.low = str2num(get(hObject,'String'));

%redraw slices if needed (if corresponding checkbox is on)
if get(findobj('Tag','bandpass_xy'),'Value') == 1
    render_slice('xy');
    render_measuremarks();
end
if get(findobj('Tag','bandpass_xz'),'Value') == 1
    render_slice('xz');
    render_measuremarks();
end
if get(findobj('Tag','bandpass_yz'),'Value') == 1
    render_slice('yz');
    render_measuremarks();
end

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Bandpass filter high                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bandpass_high_Callback(hObject, eventdata, handles)

global storage;

storage.filter.high = str2num(get(hObject,'String'));

%redraw slices if needed (if corresponding checkbox is on)
if get(findobj('Tag','bandpass_xy'),'Value') == 1
    render_slice('xy');
    render_measuremarks();
end
if get(findobj('Tag','bandpass_xz'),'Value') == 1
    render_slice('xz');
    render_measuremarks();
end
if get(findobj('Tag','bandpass_yz'),'Value') == 1
    render_slice('yz');
    render_measuremarks();
end

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Subvolume                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_subvolume_Callback(hObject, eventdata, handles)

global storage;

unset_callbacks();

storage.subvol.data = [];

waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
finalRect = rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions

x = 0;
y = 0;
z = 0;
z2 = 0;
fy2 = 0;
fx2 = 0;

orientation = get(gco,'Tag');

if strcmp(orientation,'xy_image') == 1
	x = [p1(1) p1(1)+offset(1)];
	y = [p1(2) p1(2)+offset(2)];
	x=round(x);
	y=round(y);
	
	rectangle('Position',[x(1) y(1) abs(x(2) - x(1)) abs(y(2) - y(1))],'Edgecolor',[0 1 0], 'Linewidth',1,'Tag','subvol_rect','Parent',gca);
	axes(findobj('Tag','xz_axes'));
	line('Xdata',[x(1) x(2)],'Ydata',[storage.position.xy storage.position.xy],'LineWidth',1,'Color',[0 1 0]);
	axes(findobj('Tag','yz_axes'));
	line('Xdata',[storage.position.xy storage.position.xy],'Ydata',[y(1) y(2)],'LineWidth',1,'Color',[0 1 0]);

	while (strcmp(get(gco,'Tag'),'xz_image') == false & strcmp(get(gco,'Tag'),'yz_image') == false)
				
		k=waitforbuttonpress; %waiting for the 1st point
		if k == 0 %mouse button press
    			point1=get(gca,'Currentpoint');
    		
                pt = point1(1,1:2);
    			x1=pt(1)*100;y1=pt(2)*100;
    			x1=round(x1)/100;y1=round(y1)/100;
    			if strcmp(get(gco,'Tag'),'yz_image') == 1 
				y1 = ((y(2) - y(1)) / 2) + y(1);
				h=drawmark(x1,y1,[0 1 0]);
				z1 = round(x1);
			elseif strcmp(get(gco,'Tag'),'xz_image') == 1 
				x1 = ((x(2) - x(1)) / 2) + x(1);
				h=drawmark(x1,y1,[0,1,0]);
				z1 = round(y1);
			end
		
		end
		lastfig = get(gco, 'Tag');
	end
	
	while (z2 == 0)
		
		k=waitforbuttonpress;%waiting for the 2nd point
		if k==0 %mouse button press
    			point1=get(gca,'Currentpoint');
    			pt = point1(1,1:2);
    			x2=pt(1)*100;y2=pt(2)*100;
    			x2=round(x2)/100;y2=round(y2)/100;
    			%drawmark(x2,y2);
					
			if (strcmp(get(gco,'Tag'),'yz_image') == true & strcmp(lastfig,'yz_image') == true)
				z2 = round(x2);
			elseif (strcmp(get(gco,'Tag'),'xz_image') == true & strcmp(lastfig,'xz_image') == true)
				z2 = round(y2);
			end
			
		end
	
	end

	if z1 < z2
		z = [z1 z2];
	else
		z = [z2 z1];
	end
			
elseif strcmp(orientation,'xz_image') == 1
	x = [p1(1) p1(1)+offset(1)];
	z = [p1(2) p1(2)+offset(2)];
	x=round(x);
	z=round(z);
	
	rectangle('Position',[x(1) z(1) abs(x(2) - x(1)) abs(z(2) - z(1))],'Edgecolor',[0 1 0], 'Linewidth',1,'Tag','subvol_rect','Parent',gca);
	axes(findobj('Tag','xy_axes'));
	line('Xdata',[x(1) x(2)],'Ydata',[storage.position.xz storage.position.xz],'LineWidth',1,'Color',[0 1 0]);
	axes(findobj('Tag','yz_axes'));
	line('Xdata',[z(1) z(2)],'Ydata',[storage.position.xz storage.position.xz],'LineWidth',1,'Color',[0 1 0]);

	while (strcmp(get(gco,'Tag'),'xy_image') == false & strcmp(get(gco,'Tag'),'yz_image') == false)
				
		k=waitforbuttonpress; %waiting for the 1st point
		if k==0 %mouse button press
    			point1=get(gca,'Currentpoint');
    		
		    	pt = point1(1,1:2);
    			x1=pt(1)*100;y1=pt(2)*100;
    			x1=round(x1)/100;y1=round(y1)/100;
    			if strcmp(get(gco,'Tag'),'xy_image') == 1
				x1 = ((x(2) - x(1)) / 2) + x(1);
				h=drawmark(x1,y1,[0 1 0]);
				fy1 = round(y1);
			elseif strcmp(get(gco,'Tag'),'yz_image') == 1
				x1 = ((z(2) - z(1)) / 2) + z(1);
				h=drawmark(x1,y1,[0,1,0]);
				fy1 = round(y1);
			end
				
		
		end
		lastfig = get(gco, 'Tag');
	end
	
	while (fy2 == 0)
		
		k=waitforbuttonpress;%waiting for the 2nd point
		if k==0 %mouse button press
    			point1=get(gca,'Currentpoint');
    			pt = point1(1,1:2);
    			x2=pt(1)*100;y2=pt(2)*100;
    			x2=round(x2)/100;y2=round(y2)/100;
    			%drawmark(x2,y2);
					
			if (strcmp(get(gco,'Tag'),'xy_image') == true & strcmp(lastfig,'xy_image') == true)
				fy2 = round(y2);
			elseif (strcmp(get(gco,'Tag'),'yz_image') == true & strcmp(lastfig,'yz_image') == true)
				fy2 = round(y2);
			end
			
		end
	
	end
	
	if fy1 < fy2
		y = [fy1 fy2];
	else
		y = [fy2 fy1];
	end
	
	
elseif orientation == 'yz_image'
	z = [p1(1) p1(1)+offset(1)];
	y = [p1(2) p1(2)+offset(2)];
	z=round(z);
	y=round(y);
	
	rectangle('Position',[z(1) y(1) abs(z(2) - z(1)) abs(y(2) - y(1))],'Edgecolor',[0 1 0], 'Linewidth',1,'Tag','subvol_rect','Parent',gca);
	axes(findobj('Tag','xy_axes'));
	line('Xdata',[storage.position.yz storage.position.yz],'Ydata',[y(1) y(2)],'LineWidth',1,'Color',[0 1 0]);
	axes(findobj('Tag','xz_axes'));
	line('Xdata',[storage.position.yz storage.position.yz],'Ydata',[z(1) z(2)],'LineWidth',1,'Color',[0 1 0]);

	while (strcmp(get(gco,'Tag'),'xy_image') == false & strcmp(get(gco,'Tag'),'xz_image') == false)
				
		k=waitforbuttonpress; %waiting for the 1st point
		if k==0 %mouse button press
    			point1=get(gca,'Currentpoint');
    		
			pt = point1(1,1:2);
    			x1=pt(1)*100;y1=pt(2)*100;
    			x1=round(x1)/100;y1=round(y1)/100;
    			if strcmp(get(gco,'Tag'), 'xy_image') == 1
				y1 = ((y(2) - y(1)) / 2) + y(1);
				h=drawmark(x1,y1,[0,1,0]);
				fx1 = round(x1);
			elseif strcmp(get(gco,'Tag'),'xz_image') == 1
				y1 = ((y(2) - y(1)) / 2) + y(1);
				h=drawmark(x1,y1,[0,1,0]);
				fx1 = round(x1);
			end
				
		
		end
		lastfig = get(gco, 'Tag');
	end
	
	while (fx2 == 0)
		
		k = waitforbuttonpress;%waiting for the 2nd point
		if k==0 %mouse button press
    			point1=get(gca,'Currentpoint');
    			pt = point1(1,1:2);
    			x2=pt(1)*100;y2=pt(2)*100;
    			x2=round(x2)/100;y2=round(y2)/100;
    			%drawmark(x2,y2);
					
			if (strcmp(get(gco,'Tag'),'xy_image') == true & strcmp(lastfig,'xy_image') == true)
				fx2 = round(x2);
			elseif (strcmp(get(gco,'Tag'),'xz_image') == true & strcmp(lastfig,'xz_image') == true)
				fx2 = round(x2);
			end
			
		end
	
	end
	
	if fx1 < fx2
		x = [fx1 fx2];
	else
		x = [fx2 fx1];
	end
	
end


storage.subvol.data = [x y z];

redraw_all_slices();

%Write the subvol data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subvol = storage.Volume.Value(storage.subvol.data(1):storage.subvol.data(2),storage.subvol.data(3):storage.subvol.data(4),storage.subvol.data(5):storage.subvol.data(6));
tom_emwrite(subvol);
clear subvol;
storage.subvol.data = [];

redraw_all_slices();

set_callbacks();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Subvolume reset                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_subvolumereset_Callback(hObject, eventdata, handles)

global storage;

storage.subvol.data = [];
redraw_all_slices();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Navigation slider X                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function navslider_x_Callback(hObject, eventdata, handles)

global storage;

sliceyz=round(get(hObject,'Value'));

if sliceyz > storage.dimensions.x
    sliceyz = storage.dimensions.x;
    set(hObject,'Value',sliceyz);
elseif sliceyz < 1
    sliceyz = 1;
    set(hObject,'Value',sliceyz);
end

set(findobj('Tag','nav_x'),'String',num2str(sliceyz));
storage.position.yz = sliceyz;

render_slice('yz');
render_3dbox();
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Navigation slider Y                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function navslider_y_Callback(hObject, eventdata, handles)

global storage;

slicexz=round(get(hObject,'Value'));

if slicexz > storage.dimensions.y
    slicexz = storage.dimensions.y;
    set(hObject,'Value',slicexz);
elseif slicexz < 1
    slicexz = 1;
    set(hObject,'Value',slicexz);
end

set(findobj('Tag','nav_y'),'String',num2str(slicexz));
storage.position.xz = slicexz;

render_slice('xz');
render_3dbox();
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Navigation slider Z                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function navslider_z_Callback(hObject, eventdata, handles)

global storage;

slicexy=round(get(hObject,'Value'));

if slicexy > storage.dimensions.z
    slicexy = storage.dimensions.z;
    set(hObject,'Value',slicexy);
elseif slicexy < 1
    slicexy = 1;
    set(hObject,'Value',slicexy);
end

set(findobj('Tag','nav_z'),'String',num2str(slicexy));
storage.position.xy = slicexy;

render_slice('xy');
render_3dbox();
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Navigation number box X                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nav_x_Callback(hObject, eventdata, handles)

global storage;

editx=round(str2num(get(hObject,'String')));
set(hObject,'Value',editx);
set(findobj('Tag','navslider_x'),'Value',editx);

storage.position.yz = editx;
render_slice('yz');
render_3dbox();
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Navigation number box Y                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nav_y_Callback(hObject, eventdata, handles)

global storage;

edity=round(str2num(get(hObject,'String')));
set(hObject,'Value',edity);
set(findobj('Tag','navslider_y'),'Value',edity);

storage.position.xz = edity;
render_slice('xz');
render_3dbox();
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Navigation number box Z                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nav_z_Callback(hObject, eventdata, handles)

global storage;

editz=round(str2num(get(hObject,'String')));
set(hObject,'Value',editz);
set(findobj('Tag','navslider_z'),'Value',editz);

storage.position.xy = editz;
render_slice('xy');
render_3dbox();
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average slider X                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function averageslider_x_Callback(hObject, eventdata, handles)

global storage;

sliceavg=round(get(hObject,'Value'));
set(findobj('Tag','average_x'),'String',num2str(sliceavg));
storage.average.x = sliceavg;

render_slice('yz');
render_3dbox();
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average slider Y                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function averageslider_y_Callback(hObject, eventdata, handles)

global storage;

sliceavg=round(get(hObject,'Value'));
set(findobj('Tag','average_y'),'String',num2str(sliceavg));
storage.average.y = sliceavg;

render_slice('xz');
render_3dbox();
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average slider Z                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function averageslider_z_Callback(hObject, eventdata, handles)

global storage;

sliceavg=round(get(hObject,'Value'));
set(findobj('Tag','average_z'),'String',num2str(sliceavg));
storage.average.z = sliceavg;

render_slice('xy');
render_3dbox();
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average number box X                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function average_x_Callback(hObject, eventdata, handles)

global storage;

sliceavg=round(str2num(get(hObject,'String')));

if sliceavg < 1 
    sliceavg = 1;
elseif sliceavg > storage.dimensions.x 
    sliceavg = storage.dimensions.x;
end

set(hObject,'String',num2str(sliceavg));
set(findobj('Tag','averageslider_x'),'Value',sliceavg);
storage.average.x = sliceavg;

render_slice('yz');
render_3dbox();
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average number box Y                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function average_y_Callback(hObject, eventdata, handles)

global storage;

sliceavg=round(str2num(get(hObject,'String')));

if sliceavg < 1 
    sliceavg = 1;
elseif sliceavg > storage.dimensions.y 
    sliceavg = storage.dimensions.y;
end

set(hObject,'String',num2str(sliceavg));
set(findobj('Tag','averageslider_y'),'Value',sliceavg);
storage.average.y = sliceavg;

render_slice('xz');
render_3dbox();
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average number box Z                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function average_z_Callback(hObject, eventdata, handles)

global storage;

sliceavg=round(str2num(get(hObject,'String')));

if sliceavg < 1 
    sliceavg = 1;
elseif sliceavg > storage.dimensions.z 
    sliceavg = storage.dimensions.z;
end

set(hObject,'String',num2str(sliceavg));
set(findobj('Tag','averageslider_z'),'Value',sliceavg);
storage.average.z = sliceavg;

render_slice('xy');
render_3dbox();
render_measuremarks();

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Create Movie button                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_movie_Callback(hObject, eventdata, handles)

global storage;

tom_makemoviegui(storage);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Create Screenshot button                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_screenshot_Callback(hObject, eventdata, handles)

shot = frame2im(getframe(findobj('Tag','tom_volxyz_imagewindow')));
[filename, pathname] = uiputfile({'*.tif';'*.png';'*.jpg'}, 'Save as');
if ischar(filename)
    imwrite(shot,[pathname '/' filename]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Copy to Workspace button                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_copytows_Callback(hObject, eventdata, handles)

global storage;

var_base = cellstr(evalin('base','who'));
var_exp='vol_'; rnb_vol=1;
if ~isempty(var_base)
    while 1
        a=mean(strcmp(var_base,[var_exp num2str(rnb_vol)]));
        if a==0;break;
        else;rnb_vol=rnb_vol+1;
        end
    end
end
var1=[var_exp num2str(rnb_vol)];
assignin('base',var1,storage.Volume);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Exit                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_exit_Callback(hObject, eventdata, handles)

closefcn(hObject, eventdata, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Save                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_settings_save_Callback(hObject, eventdata, handles)

global storage;

[filename, pathname] = uiputfile('*.mat', 'Save Current Settings as');
if ischar(filename)
    storagesave = storage;
    storagesave.slices.xy = [];
    storagesave.slices.xz = [];
    storagesave.slices.yz = [];
    storagesave.Volume.Value = [];
    storagesave.Volume_Orig.Value = [];

    %Save window positions
    storagesave.figurepositions.main = get(findobj('Tag','tom_volxyz'),'Position');
    storagesave.figurepositions.imagewindow = get(findobj('Tag','tom_volxyz_imagewindow'),'Position');

    %Save rotation settings
    storagesave.rotation.phi = round(get(findobj('Tag','rotateslider_phi'),'Value'));
    storagesave.rotation.psi = round(get(findobj('Tag','rotateslider_psi'),'Value'));
    storagesave.rotation.theta = round(get(findobj('Tag','rotateslider_theta'),'Value'));
    storagesave.rotation.xaxis = get(findobj('Tag','rotate_xaxis'),'Value');
    storagesave.rotation.yaxis = get(findobj('Tag','rotate_yaxis'),'Value');
    storagesave.rotation.zaxis = get(findobj('Tag','rotate_zaxis'),'Value');

    save([pathname '/' filename], 'storagesave');
    disp('Settings Saved');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Load                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_settings_load_Callback(hObject, eventdata, handles)

global storage;
[filename, pathname] = uigetfile('*.mat', 'Load Current Settings');
if ischar(filename)
    s = load([pathname '/' filename]);

    if isfield(s,'storagesave') ~= 1
        errordlg('This is not a settings file.','File Error');
    elseif s.storagesave.dimensions.x ~= storage.dimensions.x | s.storagesave.dimensions.y ~= storage.dimensions.y | s.storagesave.dimensions.z ~= storage.dimensions.z
        errordlg('Settings File is not suited for current volume.','File Error');
    else

        s.storagesave.Volume.Value = storage.Volume.Value;
        if isfield(storage,'Volume_Orig') == 1
            s.storagesave.Volume_Orig = storage.Volume_Orig;
        end

        storage = s.storagesave;

        %Restore histogram controls
        set(findobj('Tag','histogram_low'),'String',num2str(storage.DataScale(1)));
        set(findobj('Tag','histogram_high'),'String',num2str(storage.DataScale(2)));

        %Restore bandpass controls
        set(findobj('Tag','bandpass_low'),'String',num2str(storage.filter.low));
        set(findobj('Tag','bandpass_high'),'String',num2str(storage.filter.high));
        set(findobj('Tag','bandpass_xy'),'Value', storage.filter.xy);
        set(findobj('Tag','bandpass_xz'),'Value', storage.filter.xz);
        set(findobj('Tag','bandpass_yz'),'Value', storage.filter.yz);

        %Restore display options
        set(findobj('Tag','display_scalebar'),'Value', storage.display.scalebar);
        if storage.objpixelsize ~= 0
            set(findobj('Tag','scalebar_factor'),'String', num2str(storage.objpixelsize));
        end
        set(findobj('Tag','display_pixelinfo'),'Value', storage.display.pixelinfo);
        set(findobj('Tag','display_ticks'),'Value', storage.display.ticks);
        set(findobj('Tag','display_slicenumbers'),'Value', storage.display.slicenumbers);
        set(findobj('Tag','display_3dbox'),'Value', storage.display.box);

        %Restore navigation
        set(findobj('Tag','navslider_x'),'Value',storage.position.yz);
        set(findobj('Tag','navslider_y'),'Value',storage.position.xz);
        set(findobj('Tag','navslider_z'),'Value',storage.position.xy);
        set(findobj('Tag','nav_x'),'String',num2str(storage.position.yz));
        set(findobj('Tag','nav_y'),'String',num2str(storage.position.xz));
        set(findobj('Tag','nav_z'),'String',num2str(storage.position.xy));

        %Restore average
        set(findobj('Tag','averageslider_x'),'Value',storage.average.x);
        set(findobj('Tag','averageslider_y'),'Value',storage.average.y);
        set(findobj('Tag','averageslider_z'),'Value',storage.average.z);
        set(findobj('Tag','average_x'),'String',num2str(storage.average.x));
        set(findobj('Tag','average_y'),'String',num2str(storage.average.y));
        set(findobj('Tag','average_z'),'String',num2str(storage.average.z));

        set(findobj('Tag','tom_volxyz'),'Position',storage.figurepositions.main);
        set(findobj('Tag','tom_volxyz_imagewindow'),'Position',storage.figurepositions.imagewindow);

        %Restore rotation
        set(findobj('Tag','rotateslider_phi'),'Value',storage.rotation.phi);
        set(findobj('Tag','rotateslider_psi'),'Value',storage.rotation.psi);    
        set(findobj('Tag','rotateslider_theta'),'Value',storage.rotation.theta);
        set(findobj('Tag','rotate_phi'),'String',num2str(storage.rotation.phi));
        set(findobj('Tag','rotate_psi'),'String',num2str(storage.rotation.psi));
        set(findobj('Tag','rotate_theta'),'String',num2str(storage.rotation.theta));
        set(findobj('Tag','rotate_xaxis'),'Value',storage.rotation.xaxis);
        set(findobj('Tag','rotate_yaxis'),'Value',storage.rotation.yaxis);
        set(findobj('Tag','rotate_zaxis'),'Value',storage.rotation.zaxis);        

        %Restore measurement
        if ~isempty(storage.measurement.point1)
            set(findobj('Tag','measuretable_x1'),'String',sprintf('%0.1f',storage.measurement.point1(1)));
            set(findobj('Tag','measuretable_y1'),'String',sprintf('%0.1f',storage.measurement.point1(2)));
            set(findobj('Tag','measuretable_z1'),'String',sprintf('%0.1f',storage.measurement.point1(3)));
        end
        if ~isempty(storage.measurement.point2)
            set(findobj('Tag','measuretable_x2'),'String',sprintf('%0.1f',storage.measurement.point2(1)));
            set(findobj('Tag','measuretable_y2'),'String',sprintf('%0.1f',storage.measurement.point2(2)));
            set(findobj('Tag','measuretable_z2'),'String',sprintf('%0.1f',storage.measurement.point2(3)));
        end
        if ~isempty(storage.measurement.point3)
            set(findobj('Tag','measuretable_x3'),'String',sprintf('%0.1f',storage.measurement.point3(1)));
            set(findobj('Tag','measuretable_y3'),'String',sprintf('%0.1f',storage.measurement.point3(2)));
            set(findobj('Tag','measuretable_z3'),'String',sprintf('%0.1f',storage.measurement.point3(3)));
        end
        %Deactivate or activate the rotation sliders as needed
        if storage.rotation.xaxis == 1

            set(findobj('Tag','rotate_yaxis'),'Value',0);
            set(findobj('Tag','rotate_zaxis'),'Value',0);

            set(findobj('Tag','rotateslider_phi'),'Value',0,'Enable','off');
            set(findobj('Tag','rotateslider_psi'),'Value',0,'Enable','off');
            set(findobj('Tag','rotateslider_theta'),'Value',0,'Enable','on');

            set(findobj('Tag','rotate_phi'),'String','0','Enable','off');
            set(findobj('Tag','rotate_psi'),'String','0','Enable','off');
            set(findobj('Tag','rotate_theta'),'String','0','Enable','on');

        elseif storage.rotation.yaxis == 1

            set(findobj('Tag','rotate_xaxis'),'Value',0);
            set(findobj('Tag','rotate_zaxis'),'Value',0);

            set(findobj('Tag','rotateslider_phi'),'Value',270,'Enable','off');
            set(findobj('Tag','rotateslider_psi'),'Value',90,'Enable','off');
            set(findobj('Tag','rotateslider_theta'),'Enable','on');

            set(findobj('Tag','rotate_phi'),'String','270','Enable','off');
            set(findobj('Tag','rotate_psi'),'String','90','Enable','off');
            set(findobj('Tag','rotate_theta'),'Enable','on');

        elseif storage.rotation.zaxis == 1

            set(findobj('Tag','rotate_xaxis'),'Value',0);
            set(findobj('Tag','rotate_yaxis'),'Value',0);

            set(findobj('Tag','rotateslider_phi'),'Enable','on');
            set(findobj('Tag','rotateslider_psi'),'Value',0,'Enable','off');
            set(findobj('Tag','rotateslider_theta'),'Value',0,'Enable','off');

            set(findobj('Tag','rotate_phi'),'Enable','on');
            set(findobj('Tag','rotate_psi'),'String','0','Enable','off');
            set(findobj('Tag','rotate_theta'),'String','0','Enable','off');

        else

            set(findobj('Tag','rotateslider_phi'),'Value',0,'Enable','on');
            set(findobj('Tag','rotateslider_psi'),'Value',0,'Enable','on');
            set(findobj('Tag','rotateslider_theta'),'Value',0,'Enable','on');

            set(findobj('Tag','rotate_phi'),'String','0','Enable','on');
            set(findobj('Tag','rotate_psi'),'String','0','Enable','on');
            set(findobj('Tag','rotate_theta'),'String','0','Enable','on');

        end

        if storage.rotation.phi ~= 0 | storage.rotation.psi ~= 0 | storage.rotation.theta ~= 0 & ~(storage.rotation.phi == 270 & storage.rotation.psi == 90 & storage.rotation.theta == 0)
            if isfield(storage,'Volume_Orig') == 0 | isempty(storage.Volume_Orig.Value)
                storage.Volume_Orig = storage.Volume;
            end
            storage.Volume.Value = single(tom_rotate(single(storage.Volume_Orig.Value), [storage.rotation.phi, storage.rotation.psi, storage.rotation.theta], 'linear'));
        end

        redraw_all_slices();

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Display slicenumbers checkbox                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_slicenumbers_Callback(hObject, eventdata, handles)

global storage;

if get(hObject, 'Value') == 1
   
    storage.display.slicenumbers = 1;
    
    axes(findobj('Tag','xy_axes'));
    delete(findobj('Tag','slicenumber_xy'));
    t = text(storage.actualaxis.xy(2)-(storage.actualaxis.xy(2)/100)*10,storage.actualaxis.xy(4)-(storage.actualaxis.xy(4)/100)*8,num2str(storage.position.xy));
    set(t,'Tag','slicenumber_xy','Color',[1 0 0],'FontWeight','bold','FontUnits','Normalized','FontSize',0.11,'HorizontalAlignment','right');
    set(t,'FontUnits','Pixel');
    storage.display.slicefontsize = get(t,'FontSize');
    
    axes(findobj('Tag','xz_axes'));
    delete(findobj('Tag','slicenumber_xz'));
    t = text(storage.actualaxis.xz(2)-(storage.actualaxis.xz(2)/100)*10,storage.actualaxis.xz(4)-(storage.actualaxis.xz(4)/100)*8,num2str(storage.position.xz));
    set(t,'Tag','slicenumber_xz','Color',[1 0 0],'FontWeight','bold','FontUnits','Pixel','FontSize',storage.display.slicefontsize,'HorizontalAlignment','right');
        
    axes(findobj('Tag','yz_axes'));
    delete(findobj('Tag','slicenumber_yz'));
    t = text(storage.actualaxis.yz(2)-(storage.actualaxis.yz(2)/100)*10,storage.actualaxis.yz(4)-(storage.actualaxis.yz(4)/100)*8,num2str(storage.position.yz));
   	set(t,'Tag','slicenumber_yz','Color',[1 0 0],'FontWeight','bold','FontUnits','Pixel','FontSize',storage.display.slicefontsize,'HorizontalAlignment','right');	
    
else

    storage.display.slicenumbers = 0;
    
    delete(findobj('Tag','slicenumber_xy'));
    delete(findobj('Tag','slicenumber_xz'));
    delete(findobj('Tag','slicenumber_yz'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Display Pixel info checkbox                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_pixelinfo_Callback(hObject, eventdata, handles)

global storage;

if get(hObject,'Value') == 1 
    storage.display.pixelinfo = 1;
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
else
    storage.display.pixelinfo = 0;
    delete(storage.display.pixelinfohandle);
end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Display 3Dbox checkbox                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_3dbox_Callback(hObject, eventdata, handles)

global storage;

if get(hObject,'Value') == 1 
    storage.display.box = 1;
else
    storage.display.box = 0;
end

render_3dbox();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Display Axis ticks checkbox                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_ticks_Callback(hObject, eventdata, handles)

global storage;

if get(hObject,'Value') == 1
    storage.display.ticks = 1;
else
    storage.display.ticks = 0;
end

redraw_all_slices();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Display cross hair checkbox                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_crosshair_Callback(hObject, eventdata, handles)

global storage;

if get(hObject,'Value') == 1
    storage.display.crosshair = 1;
else
    storage.display.crosshair = 0;
end

redraw_all_slices();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Display scalebar checkbox                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_scalebar_Callback(hObject, eventdata, handles)

global storage;

if get(hObject,'Value') == 1
    if storage.objpixelsize == 0
        set(hObject,'Value',0);
        msgbox('Define object pixel size first!','Error','error');
    else
        storage.display.scalebar = 1;
    end
else
    storage.display.scalebar = 0;
end

render_scalebar();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  scalefactor input field                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scalebar_factor_Callback(hObject, eventdata, handles)

global storage;

storage.objpixelsize = str2num(get(hObject,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rotateslider phi                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotateslider_phi_Callback(hObject, eventdata, handles)

global storage;

angle=round(get(hObject,'Value'));
set(findobj('Tag','rotate_phi'),'String',num2str(angle));
storage.rotate.phi = angle;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rotateslider psi                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotateslider_psi_Callback(hObject, eventdata, handles)

global storage;

angle=round(get(hObject,'Value'));
set(findobj('Tag','rotate_psi'),'String',num2str(angle));
storage.rotate.psi = angle;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rotateslider theta                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotateslider_theta_Callback(hObject, eventdata, handles)

global storage;

angle=round(get(hObject,'Value'));
set(findobj('Tag','rotate_theta'),'String',num2str(angle));
storage.rotate.theta = angle;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rotate phi input field                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotate_phi_Callback(hObject, eventdata, handles)

global storage;

angle=round(str2num(get(hObject,'String')));

set(hObject,'String',num2str(angle));
set(findobj('Tag','rotateslider_phi'),'Value',angle);
storage.rotate.phi = angle;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rotate psi input field                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotate_psi_Callback(hObject, eventdata, handles)

global storage;

angle=round(str2num(get(hObject,'String')));

set(hObject,'String',num2str(angle));
set(findobj('Tag','rotateslider_psi'),'Value',angle);
storage.rotate.psi = angle;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rotate theta input field                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotate_theta_Callback(hObject, eventdata, handles)

global storage;

angle=round(str2num(get(hObject,'String')));

set(hObject,'String',num2str(angle));
set(findobj('Tag','rotateslider_theta'),'Value',angle);
storage.rotate.theta = angle;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rotate x-Axis checkbox                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotate_xaxis_Callback(hObject, eventdata, handles)

global storage;

if get(hObject,'Value') == 1
    set(findobj('Tag','rotate_yaxis'),'Value',0);
    set(findobj('Tag','rotate_zaxis'),'Value',0);

    set(findobj('Tag','rotateslider_phi'),'Value',0,'Enable','off');
    set(findobj('Tag','rotateslider_psi'),'Value',0,'Enable','off');
    set(findobj('Tag','rotateslider_theta'),'Value',0,'Enable','on');

    set(findobj('Tag','rotate_phi'),'String','0','Enable','off');
    set(findobj('Tag','rotate_psi'),'String','0','Enable','off');
    set(findobj('Tag','rotate_theta'),'String','0','Enable','on');
    
    storage.rotate.phi = 0;
    storage.rotate.psi = 0;
else
    set(findobj('Tag','rotateslider_phi'),'Enable','on');
    set(findobj('Tag','rotateslider_psi'),'Enable','on');
    set(findobj('Tag','rotate_phi'),'Enable','on');
    set(findobj('Tag','rotate_psi'),'Enable','on');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rotate y-Axis checkbox                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotate_yaxis_Callback(hObject, eventdata, handles)

global storage;

if get(hObject,'Value') == 1
    set(findobj('Tag','rotate_xaxis'),'Value',0);
    set(findobj('Tag','rotate_zaxis'),'Value',0);

    set(findobj('Tag','rotateslider_phi'),'Value',270,'Enable','off');
    set(findobj('Tag','rotateslider_psi'),'Value',90,'Enable','off');
    set(findobj('Tag','rotateslider_theta'),'Enable','on');
    
    set(findobj('Tag','rotate_phi'),'String','270','Enable','off');
    set(findobj('Tag','rotate_psi'),'String','90','Enable','off');
    set(findobj('Tag','rotate_theta'),'Enable','on');
    
    storage.rotate.phi = 270;
    storage.rotate.psi = 90;
else
    set(findobj('Tag','rotateslider_phi'),'Enable','on');
    set(findobj('Tag','rotateslider_psi'),'Enable','on');
    set(findobj('Tag','rotate_phi'),'Enable','on');
    set(findobj('Tag','rotate_psi'),'Enable','on');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rotate z-Axis checkbox                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotate_zaxis_Callback(hObject, eventdata, handles)

global storage;

if get(hObject,'Value') == 1
    set(findobj('Tag','rotate_xaxis'),'Value',0);
    set(findobj('Tag','rotate_yaxis'),'Value',0);

    set(findobj('Tag','rotateslider_phi'),'Enable','on');
    set(findobj('Tag','rotateslider_psi'),'Value',0,'Enable','off');
    set(findobj('Tag','rotateslider_theta'),'Value',0,'Enable','off');

    set(findobj('Tag','rotate_phi'),'Enable','on');
    set(findobj('Tag','rotate_psi'),'String','0','Enable','off');
    set(findobj('Tag','rotate_theta'),'String','0','Enable','off');

    storage.rotate.psi = 0;
    storage.rotate.theta = 0;
else
    set(findobj('Tag','rotateslider_psi'),'Enable','on');
    set(findobj('Tag','rotateslider_theta'),'Enable','on');
    set(findobj('Tag','rotate_psi'),'Enable','on');
    set(findobj('Tag','rotate_theta'),'Enable','on');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  button rotate                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_rotate_Callback(hObject, eventdata, handles)

global storage;

if isfield(storage,'Volume_Orig') == 0 | isempty(storage.Volume_Orig.Value)
    storage.Volume_Orig = storage.Volume;
end
% added by SN
%storage.Volume.Value = single(tom_rotate(single(storage.Volume_Orig.Value),[storage.rotate.phi, storage.rotate.psi, storage.rotate.theta],'linear','taper'));
storage.Volume.Value = single(tom_rotate(single(storage.Volume_Orig.Value),[storage.rotate.phi, storage.rotate.psi, storage.rotate.theta],'linear'));
redraw_all_slices();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  button apply rotation                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_applyrotation_Callback(hObject, eventdata, handles)

global storage;

storage.Volume_Orig.Value = storage.Volume.Value;

set(findobj('Tag','rotateslider_phi'),'Value',0);
set(findobj('Tag','rotateslider_psi'),'Value',0);
set(findobj('Tag','rotateslider_theta'),'Value',0);

set(findobj('Tag','rotate_phi'),'String','0');
set(findobj('Tag','rotate_psi'),'String','0');
set(findobj('Tag','rotate_theta'),'String','0');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  button shift                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_shift_Callback(hObject, eventdata, handles)

global storage;

shiftx = -(storage.position.yz-(floor(storage.dimensions.x./2)+1));
shifty = -(storage.position.xz-(floor(storage.dimensions.y./2)+1));
shiftz = -(storage.position.xy-(floor(storage.dimensions.z./2)+1));

storage.Volume.Value = single(tom_shift(storage.Volume.Value, [shiftx, shifty, shiftz]));

center_position();
redraw_all_slices();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  button pixel profile                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_profile_Callback(hObject, eventdata, handles)

global storage;

thickness = str2num(get(findobj('Tag','profile_thickness'),'String'));

%clear previous profile drawings
delete(findobj('Tag','tom_volxyz_profile'));
%delete(findobj('Tag','profileline'));

axes(findobj('Tag',storage.profile.plane));

[xa,ya,a,state] = getimage;
prof = getline(gcf); % Get profile from user

% Parametric distance along segments
s = [0;cumsum(sqrt(sum((diff(prof).^2)')'))];

% Remove duplicate points if necessary.
killIdx = find(diff(s) == 0);
if (~isempty(killIdx))
    s(killIdx+1) = [];
    prof(killIdx+1,:) = [];
end

%Size of image
ma = size(a,1);
na = size(a,2);

xmin = min(xa(:)); ymin = min(ya(:));
xmax = max(xa(:)); ymax = max(ya(:));

%Size of axes
dx = max( (xmax-xmin)/(na-1), eps );  
xxa = xmin:dx:xmax;
dy = max( (ymax-ymin)/(ma-1), eps );
yya = ymin:dy:ymax;

%Determine number of points from profile
d = abs(diff(prof./(ones(size(prof,1),1)*[dx dy])));
n = max(sum(max(ceil(d)')),2); % In pixel coordinates

vec_orth_normed = zeros(size(prof,1)-1,2);
%calculate normalizes orthogonal vector to each segment
for i = 1:size(prof,1)-1
    vec_seg = [prof(i+1,1)-prof(i,1);prof(i+1,2)-prof(i,2)];
    vec_orth = [vec_seg(2);-vec_seg(1)];
    abs_vec_orth = sqrt((vec_orth(1)^2+vec_orth(2)^2));
    vec_orth_normed(i,:) = vec_orth ./ abs_vec_orth;
end
%Integrate over interpolated lines
factor = zeros(thickness,1);
lauf = 1;
for i = 2:thickness
    if mod(i,2) == 1
        factor(i) = -lauf;
        lauf = lauf + 1;
    else
        factor(i) = lauf;
    end
end

for j = 1:thickness
    % Interpolation points along segments
    if ~isempty(prof)
        %move all points by shift vectors
        prof2 = zeros(size(prof));
        prof2(1,:) = prof(1,:) + (vec_orth_normed(1,:) .* factor(j));
        for i = 2:size(prof,1)-1
            prof2(i,:) = prof(i,:) + (vec_orth_normed(i-1,:) .* factor(j)) + (vec_orth_normed(i,:) .* factor(j));
        end
        prof2(size(prof,1),:) = prof(size(prof,1),:) + (vec_orth_normed(size(prof,1)-1,:) .* factor(j));
        
        profi = interp1(s,prof2,0:(max(s)/(n-1)):max(s));
        xg = profi(:,1);
        yg = profi(:,2);
    else
        xg = []; yg = [];
    end

    if ~isempty(a) && ~isempty(xg)
        % Image values along interpolation points - the g stands for Grayscale
        zg = interp2(xxa,yya,a,xg,yg,'cubic'); %use nearest , linear, spline , cubic 
        if j == 1
            zg_integ = zg;
        else
            zg_integ = [zg_integ, zg];
        end
    end
    
end
zg = mean(zg_integ, 2);

%Draw profile line in slice
for i = 1:size(prof,1)-1
    line([prof(i,1) prof(i+1,1)],[prof(i,2) prof(i+1,2)],'Color',[0 1 0],'Tag','profileline','LineWidth',thickness);
end

% plot the profile
if ~isempty(zg) 
    fig = figure('Name','Profile','Tag','tom_volxyz_profile');
    set(fig,'CloseRequestFcn',@(h,e)button_profileclear_Callback(hObject, eventdata, handles));
    if length(prof)>2
        plot3(xg,yg,zg,'b');
        set(gca,'ydir','reverse');
        xlabel X, ylabel Y;grid on;
    else
        plot(sqrt((xg-xg(1)).^2+(yg-yg(1)).^2),zg,'b');
        xlabel('Distance along profile');
    end
end

storage.profile.valcell{end+1} = zg;
posx = zeros(size(xg,1),1)+storage.position.xy;
storage.profile.coordcell{end+1} = [xg yg posx];
valcell = storage.profile.valcell;
coordcell = storage.profile.coordcell;

save('./linescans4b.mat','valcell','coordcell');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  button pixel profile clear                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_profileclear_Callback(hObject, eventdata, handles)

delete(findobj('Tag','tom_volxyz_profile'));
delete(findobj('Tag','profileline'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  radio button profile xy plane                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function profile_xy_Callback(hObject, eventdata, handles)

global storage;

if get(hObject,'Value') == 1
   storage.profile.plane = 'xy_axes';
   set(findobj('Tag','profile_xz'),'Value',0);
   set(findobj('Tag','profile_yz'),'Value',0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  radio button profile xz plane                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function profile_xz_Callback(hObject, eventdata, handles)

global storage;

if get(hObject,'Value') == 1
   storage.profile.plane = 'xz_axes';
   set(findobj('Tag','profile_xy'),'Value',0);
   set(findobj('Tag','profile_yz'),'Value',0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  radio button profile yz plane                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function profile_yz_Callback(hObject, eventdata, handles)

global storage;

if get(hObject,'Value') == 1
   storage.profile.plane = 'yz_axes'; 
   set(findobj('Tag','profile_xy'),'Value',0);
   set(findobj('Tag','profile_xz'),'Value',0);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  radio button profile yz plane                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_profilehelp_Callback(hObject, eventdata, handles)

helpdlg('You specify the line or path using the mouse, by clicking points in the image. Press Backspace or Delete to remove the previously selected point. A shift-click, right-click, or double-click adds a final point and ends the selection; pressing Return finishes the selection without adding a point.','Help on Profile Tool');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  profile thickness                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function profile_thickness_Callback(hObject, eventdata, handles)

thickness = str2num(get(hObject,'String'));
set(hObject, 'String',num2str(round(thickness)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  measurement point 1                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_measurement_point1_Callback(hObject, eventdata, handles)

global storage;

figure(findobj('Tag','tom_volxyz_imagewindow'));
[x y] = ginput(1);

if x > 0 & y > 0 & strcmp(get(gca,'Tag'),'box') == 0 & ~isempty(get(gca,'Tag'))
    if strcmp(get(gca,'Tag'),'xy_axes') == 1
        storage.measurement.point1 = [x y storage.position.xy];
    elseif strcmp(get(gca,'Tag'),'xz_axes') == 1
        storage.measurement.point1 = [x storage.position.xz y];
    elseif strcmp(get(gca,'Tag'),'yz_axes') == 1
        storage.measurement.point1 = [storage.position.yz y x];
    end

    %update the measurement text field
    set(findobj('Tag','measuretable_x1'),'String',sprintf('%0.1f',storage.measurement.point1(1)));
    set(findobj('Tag','measuretable_y1'),'String',sprintf('%0.1f',storage.measurement.point1(2)));
    set(findobj('Tag','measuretable_z1'),'String',sprintf('%0.1f',storage.measurement.point1(3)));

    calculate_distance();
    calculate_angle();
    render_measuremarks();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  measurement point 2                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_measurement_point2_Callback(hObject, eventdata, handles)

global storage;

figure(findobj('Tag','tom_volxyz_imagewindow'));
[x y] = ginput(1);

if x > 0 & y > 0 & strcmp(get(gca,'Tag'),'box') == 0 & ~isempty(get(gca,'Tag'))
    if strcmp(get(gca,'Tag'),'xy_axes') == 1
        storage.measurement.point2 = [x y storage.position.xy];
    elseif strcmp(get(gca,'Tag'),'xz_axes') == 1
        storage.measurement.point2 = [x storage.position.xz y];
    elseif strcmp(get(gca,'Tag'),'yz_axes') == 1
        storage.measurement.point2 = [storage.position.yz y x];
    end

    %update the measurement text field
    set(findobj('Tag','measuretable_x2'),'String',sprintf('%0.1f',storage.measurement.point2(1)));
    set(findobj('Tag','measuretable_y2'),'String',sprintf('%0.1f',storage.measurement.point2(2)));
    set(findobj('Tag','measuretable_z2'),'String',sprintf('%0.1f',storage.measurement.point2(3)));

    calculate_distance();
    calculate_angle();
    render_measuremarks();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  measurement point 3                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_measurement_point3_Callback(hObject, eventdata, handles)

global storage;

figure(findobj('Tag','tom_volxyz_imagewindow'));
[x y] = ginput(1);

if x > 0 & y > 0 & strcmp(get(gca,'Tag'),'box') == 0 & ~isempty(get(gca,'Tag'))
    if strcmp(get(gca,'Tag'),'xy_axes') == 1
        storage.measurement.point3 = [x y storage.position.xy];
    elseif strcmp(get(gca,'Tag'),'xz_axes') == 1
        storage.measurement.point3 = [x storage.position.xz y];
    elseif strcmp(get(gca,'Tag'),'yz_axes') == 1
        storage.measurement.point3 = [storage.position.yz y x];
    end

    %update the measurement text field
    set(findobj('Tag','measuretable_x3'),'String',sprintf('%0.1f',storage.measurement.point3(1)));
    set(findobj('Tag','measuretable_y3'),'String',sprintf('%0.1f',storage.measurement.point3(2)));
    set(findobj('Tag','measuretable_z3'),'String',sprintf('%0.1f',storage.measurement.point3(3)));

    calculate_angle();
    render_measuremarks();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  measurement clear                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_measurement_clear_Callback(hObject, eventdata, handles)

global storage;

storage.measurement.point1 = [];
storage.measurement.point2 = [];
storage.measurement.point3 = [];

set(findobj('Tag','measuretable_x1'),'String','');
set(findobj('Tag','measuretable_y1'),'String','');
set(findobj('Tag','measuretable_z1'),'String','');

set(findobj('Tag','measuretable_x2'),'String','');
set(findobj('Tag','measuretable_y2'),'String','');
set(findobj('Tag','measuretable_z2'),'String','');

set(findobj('Tag','measuretable_x3'),'String','');
set(findobj('Tag','measuretable_y3'),'String','');
set(findobj('Tag','measuretable_z3'),'String','');

set(findobj('Tag','measurement_distance'),'String','');
set(findobj('Tag','measurement_angle'),'String','');

redraw_all_slices();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  gen_isosurface                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_2isosurface_Callback(hObject, eventdata, handles)

global storage;
    
%ice
if get(findobj('Tag','checkbox_negstain'),'Value') == 0
    tom_isosurface(-tom_norm(storage.Volume.Value,1));
%neg stain
else
    tom_isosurface(tom_norm(storage.Volume.Value,1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  isosurface 2 amira                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_2amira_Callback(hObject, eventdata, handles)

global storage;

thresh = str2double(get(findobj('Tag','edit_threshold'),'String'));

if isempty(thresh)
    errordlg('Enter threshold first!');
end

if isempty(storage.Volume.Header.Filename)
    [filename, pathname] = uiputfile({'*.em';'*.vol';'*.*'}, 'Save as EM-file');
    if isequal(filename,0) || isequal(pathname,0); disp('Data not saved.'); return; end;
    tom_emwrite([pathname '/' filename], storage.Volume);
    storage.Volume.Header.Filename = filename;
    storage.Volume.Header.Pathname = pathname;
end

tom_amira_loadfile([storage.Volume.Header.Pathname '/' storage.Volume.Header.Filename],storage.Volume.Header.Filename);
tom_amira_createisosurface(storage.Volume.Header.Filename, thresh, ['Isosurface_' storage.Volume.Header.Filename]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  copy point list 2 workspace                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_volxyz_points2workspace_Callback(hObject, eventdata, handles)

global storage;
assignin('base','points',storage.Align(1,2:end));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                     %%
%%                                                                     %%
%%  Helper functions                                                   %%
%%                                                                     %%
%%                                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate distance (measurement tool)                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calculate_distance()

global storage;

if ~isempty(storage.measurement.point1) & ~isempty(storage.measurement.point2)
    distance = sqrt((storage.measurement.point1(1) - storage.measurement.point2(1))^2 + (storage.measurement.point1(2) - storage.measurement.point2(2))^2 + (storage.measurement.point1(3) - storage.measurement.point2(3))^2);
    diststring = [sprintf('%0.1f',distance), ' Pixel'];
    if storage.objpixelsize ~= 0
        diststring = [diststring, ' (', num2str(round(distance .* storage.objpixelsize)), ' A)'];
    end
    set(findobj('Tag','measurement_distance'),'String',diststring);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate angle (measurement tool)                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calculate_angle()

global storage;

if ~isempty(storage.measurement.point1) & ~isempty(storage.measurement.point2) & ~isempty(storage.measurement.point3)
    vec1 = storage.measurement.point1 - storage.measurement.point2;
    vec2 = storage.measurement.point3 - storage.measurement.point2;
    angle = acosd((vec1(1) .* vec2(1) + vec1(2) .* vec2(2) + vec1(3) .* vec2(3)) ./ (sqrt((vec1(1)^2+vec1(2)^2+vec1(3)^2)) .* sqrt((vec2(1)^2+vec2(2)^2)+vec2(3)^2)));
    set(findobj('Tag','measurement_angle'),'String',[sprintf('%0.1f',angle) ' deg']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate histogram                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calculate_histogram()

global storage;

for lauf=1:100
    rand_pos=rand(3,1);
    rand_pos(1)=rand_pos(1) * storage.dimensions.x;
    rand_pos(2)=rand_pos(2) * storage.dimensions.y;
    rand_pos(3)=rand_pos(3) * storage.dimensions.z;
    
    rand_pos(1)=floor(rand_pos(1)+1);
    rand_pos(2)=floor(rand_pos(2)+1);
    rand_pos(3)=floor(rand_pos(3)+1);
       
    rand_pos(1)= rand_pos(1)-10;
    rand_pos(2)= rand_pos(2)-10;
    
    if rand_pos(1) < 1
    	rand_pos(1) = 1;
    end
    if rand_pos(2) < 1
    	rand_pos(2) = 1;
    end
      
    tmp = double(storage.Volume.Value(rand_pos(1):rand_pos(1)+8,rand_pos(2):rand_pos(2)+8,rand_pos(3)));
    storage.histogramdata(:,:,lauf)=tmp;
    
    [mean max min std] = tom_dev(storage.histogramdata,'noinfo');
    storage.DataScale = [mean-2.*std, mean+2.*std];
    storage.mean = mean;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  render histogram                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_histogram()

global storage;

[h,n] = tom_hist3d(storage.histogramdata);
h = 200 .* h ./ (100.*storage.dimensions.x .* storage.dimensions.y .* storage.dimensions.z);

axesobj = findobj('Tag','histogram');
axes(axesobj); 
bar(n,h); 
axis auto;
set(axesobj,'Tag','histogram');

set(findobj('Tag','histogram_low'),'String', num2str(storage.DataScale(1)));
set(findobj('Tag','histogram_high'),'String', num2str(storage.DataScale(2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  render scalebar                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_scalebar()

global storage;

if storage.display.scalebar == 1

    delete(findobj('Tag','scalebar'));
    delete(findobj('Tag','scalebar_text'));
    axes(findobj('Tag','xy_axes'));

    exactnumber = storage.objpixelsize .* ((storage.actualaxis.xy(4) - storage.actualaxis.xy(3)).*0.1);
    roundednumber = round(exactnumber ./ 5) .* 5;

    width = roundednumber ./ storage.objpixelsize;
    t = text((((storage.actualaxis.xy(4) - storage.actualaxis.xy(3)).*0.10).*2+width)./2,((storage.actualaxis.xy(2) - storage.actualaxis.xy(1)) .* .9) + (storage.actualaxis.xy(4) - storage.actualaxis.xy(3)).*0.05,[num2str(roundednumber) ' A'],'Color',[0 0 0],'Tag','scalebar_text','HorizontalAlignment','Center','FontUnits','normalized','FontSize',.05,'FontWeight','bold', 'BackgroundColor',[1 1 1],'Margin',(storage.actualaxis.xy(4) - storage.actualaxis.xy(3)).*0.05);
    l = line([(storage.actualaxis.xy(4) - storage.actualaxis.xy(3)).*0.10, ((storage.actualaxis.xy(4) - storage.actualaxis.xy(3)) .* .10)+width],[(storage.actualaxis.xy(2) - storage.actualaxis.xy(1)) .* .9, (storage.actualaxis.xy(2) - storage.actualaxis.xy(1)) .* .9],'Tag','scalebar','LineWidth',(storage.actualaxis.xy(2) - storage.actualaxis.xy(1)) .* 0.02,'Color',[0 0 0]);

else
    delete(findobj('Tag','scalebar'));
    delete(findobj('Tag','scalebar_text'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  render measurement marks/lines                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_measuremarks()

global storage;

%draw point marks if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

point1_vis = '';
point2_vis = '';
point3_vis = '';

delete(findobj('Tag','measureline'));
delete(findobj('Tag','measurepoint'));

%% point1
%%%%%%%%%%%%%%
if ~isempty(storage.measurement.point1)
    if storage.position.xy == round(storage.measurement.point1(3))
        axes(findobj('Tag','xy_axes'));
        drawmark(storage.measurement.point1(1),storage.measurement.point1(2),[0 0 1],'measurepoint','1');
        point1_vis = 'xy_axes';
        point1_coords = [storage.measurement.point1(1),storage.measurement.point1(2)];
    end
    if storage.position.xz == round(storage.measurement.point1(2))
        axes(findobj('Tag','xz_axes'));
        drawmark(storage.measurement.point1(1),storage.measurement.point1(3),[0 0 1],'measurepoint','1');
        point1_vis = 'xz_axes';
        point1_coords = [storage.measurement.point1(1),storage.measurement.point1(3)];
    end
    if storage.position.yz == round(storage.measurement.point1(1))
        axes(findobj('Tag','yz_axes'));
        drawmark(storage.measurement.point1(3),storage.measurement.point1(2),[0 0 1],'measurepoint','1');
        point1_vis = 'yz_axes';
        point1_coords = [storage.measurement.point1(3),storage.measurement.point1(2)];
    end
end

%% point2
%%%%%%%%%%%%%%
if ~isempty(storage.measurement.point2)
    if storage.position.xy == round(storage.measurement.point2(3))
        axes(findobj('Tag','xy_axes'));
        drawmark(storage.measurement.point2(1),storage.measurement.point2(2),[0 0 1],'measurepoint','2');
        point2_vis = 'xy_axes';
        point2_coords = [storage.measurement.point2(1),storage.measurement.point2(2)];
    end
    if storage.position.xz == round(storage.measurement.point2(2))
        axes(findobj('Tag','xz_axes'));
        drawmark(storage.measurement.point2(1),storage.measurement.point2(3),[0 0 1],'measurepoint','2');
        point2_vis = 'xz_axes';
        point2_coords = [storage.measurement.point2(1),storage.measurement.point2(3)];
    end
    if storage.position.yz == round(storage.measurement.point2(1))
        axes(findobj('Tag','yz_axes'));
        drawmark(storage.measurement.point2(3),storage.measurement.point2(2),[0 0 1],'measurepoint','2');
        point2_vis = 'yz_axes';
        point2_coords = [storage.measurement.point2(3),storage.measurement.point2(2)];
    end
end

%% point3
%%%%%%%%%%%%%%
if ~isempty(storage.measurement.point3)
    if storage.position.xy == round(storage.measurement.point3(3))
        axes(findobj('Tag','xy_axes'));
        drawmark(storage.measurement.point3(1),storage.measurement.point3(2),[0 0 1],'measurepoint','3');
        point3_vis = 'xy_axes';
        point3_coords = [storage.measurement.point3(1),storage.measurement.point3(2)];
    end
    if storage.position.xz == round(storage.measurement.point3(2))
        axes(findobj('Tag','xz_axes'));
        drawmark(storage.measurement.point3(1),storage.measurement.point3(3),[0 0 1],'measurepoint','3');
        point3_vis = 'xz_axes';
        point3_coords = [storage.measurement.point3(1),storage.measurement.point3(3)];
    end
    if storage.position.yz == round(storage.measurement.point3(1))
        axes(findobj('Tag','yz_axes'));
        drawmark(storage.measurement.point3(3),storage.measurement.point3(2),[0 0 1],'measurepoint','3');
        point3_vis = 'yz_axes';        
        point3_coords = [storage.measurement.point3(3),storage.measurement.point3(2)];
    end

end

%% connecting lines
%%%%%%%%%%%%%%%%%%%%%%%%
if (~isempty(storage.measurement.point1) & ~isempty(storage.measurement.point2)) | (~isempty(storage.measurement.point1) & ~isempty(storage.measurement.point3)) | (~isempty(storage.measurement.point2) & ~isempty(storage.measurement.point3))
    if strcmp(point1_vis,point2_vis) == 1 & ~isempty(point1_vis)
        axes(findobj('Tag',point1_vis));
        line('XData',[point1_coords(1),point2_coords(1)],'YData',[point1_coords(2),point2_coords(2)],'Color',[0 0 1],'Tag','measureline','LineWidth',2);
    end
    %if strcmp(point1_vis,point3_vis) == 1 & ~isempty(point1_vis)
    %    axes(findobj('Tag',point1_vis));
    %    line('XData',[point1_coords(1),point3_coords(1)],'YData',[point1_coords(2),point3_coords(2)],'Color',[0 0 1],'Tag','measureline','LineWidth',2);
    %end
    if strcmp(point2_vis,point3_vis) == 1 & ~isempty(point2_vis)
        axes(findobj('Tag',point1_vis));
        line('XData',[point2_coords(1),point3_coords(1)],'YData',[point2_coords(2),point3_coords(2)],'Color',[0 0 1],'Tag','measureline','LineWidth',2);
    end
    
%     if storage.display.box == 1
%         axes(findobj('Tag','box'));
%         if ~isempty(storage.measurement.point1) & ~isempty(storage.measurement.point2)
%             line('Xdata',[storage.measurement.point1(1),storage.measurement.point2(1)],'YData',[storage.measurement.point1(2),storage.measurement.point2(2)],'Zdata',[storage.measurement.point1(3),storage.measurement.point2(3)],'LineWidth',2,'Tag','measureline','Color',[0 0 1]);
%         end
%         if ~isempty(storage.measurement.point2) & ~isempty(storage.measurement.point3)
%             line('Xdata',[storage.measurement.point2(1),storage.measurement.point3(1)],'YData',[storage.measurement.point2(2),storage.measurement.point3(2)],'Zdata',[storage.measurement.point2(3),storage.measurement.point3(3)],'LineWidth',2,'Tag','measureline','Color',[0 0 1]);
%         end
%     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  redraw all slices                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function redraw_all_slices()

global storage;

render_slice('xy');
render_slice('yz');
render_slice('xz');

render_3dbox;

if storage.display.pixelinfo == 1
    storage.display.pixelinfohandle = impixelinfo(findobj('Tag','tom_volxyz_imagewindow'));
end

render_measuremarks();



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  redraw function for all slices                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_slice(orientation)

global storage;

% render xy slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(orientation,'xy') == 1
    
    sliceavg = storage.average.z;
    
    if storage.position.xy + sliceavg > storage.dimensions.z
    		sliceavg = storage.dimensions.z - storage.position.xy;
    elseif storage.position.xy - sliceavg < 1
		sliceavg=1; %FIXME: calculate the correct number of slices
    end
    
    if storage.position.xy > 0 & storage.position.xy <= storage.dimensions.z

        %cut out slice and make average
        slice = single(storage.Volume.Value(:,:,round(storage.position.xy - sliceavg ./ 2):round(storage.position.xy + sliceavg ./ 2) - 1));
		slice = single(mean(slice, 3));
        slice = squeeze(slice);
        storage.slices.xy = slice;        
        %bandpass
        if storage.filter.xy ~= 0
            slice = tom_bandpass(slice, storage.filter.low, storage.filter.high);
			DataScale = storage.DataScale ./ (storage.dimensions.x .* storage.dimensions.y);
            if storage.filter.low > 0
                DataScale = DataScale - storage.mean ./ (storage.dimensions.x .* storage.dimensions.y);
            end
        else
            DataScale = storage.DataScale;
        end
        
        %switchto xy axes
        axesobj = findobj('Tag','xy_axes');
        axes(axesobj); colormap(gray); axis ij;
        %create images
        imhandle = imagesc(slice',DataScale);
  		axis(storage.actualaxis.xy);
        set(imhandle, 'Tag', 'xy_image');
        set(axesobj,'Tag','xy_axes');
        if storage.display.slicenumbers == 1
            delete(findobj('Tag','slicenumber_xy'));
            t = text('Units','normalized','Position',[.9 .1],'String',num2str(storage.position.xy));
            set(t,'Tag','slicenumber_xy','Color',[1 0 0],'FontWeight','bold','FontUnits','Pixel','FontSize',storage.display.slicefontsize,'HorizontalAlignment','right');
        end

        if storage.display.ticks == 0;
            axis off;
        end
        
        render_scalebar();

    end

% render xz slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif strcmp(orientation,'xz') == 1
    
    sliceavg = storage.average.y;
    
    if storage.position.xz + sliceavg > storage.dimensions.y
    		sliceavg = storage.dimensions.y - storage.position.xz;
    elseif storage.position.xz - sliceavg < 1
		sliceavg=1; %FIXME: calculate the correct number of slices
    end
    
    if storage.position.xz > 0 & storage.position.xz <= storage.dimensions.y

        %cut out slice and make average
        slice = single(storage.Volume.Value(:,round(storage.position.xz-sliceavg./2):round(storage.position.xz+sliceavg./2)-1,:));
		slice = single(mean(slice, 2));
        slice = squeeze(slice);
        storage.slices.xz = slice;
        %bandpass
        if storage.filter.xz ~= 0
			slice = tom_bandpass(slice, storage.filter.low, storage.filter.high);
			DataScale = storage.DataScale ./ (storage.dimensions.x .* storage.dimensions.z);
            if storage.filter.low > 0
                DataScale = DataScale - storage.mean ./ (storage.dimensions.x .* storage.dimensions.y);
            end
        else
            DataScale = storage.DataScale;
        end
        
        %switchto xz axes
        axesobj = findobj('Tag','xz_axes');
        axes(axesobj); colormap(gray); axis ij;
        %create image
        imhandle = imagesc(slice',DataScale);
  		axis(storage.actualaxis.xz);
        set(imhandle, 'Tag', 'xz_image');
        set(axesobj,'Tag','xz_axes');
        if storage.display.slicenumbers == 1
            delete(findobj('Tag','slicenumber_xz'));
            t = text('Units','normalized','Position',[.9 .1],'String',num2str(storage.position.xz));
        	set(t,'Tag','slicenumber_xz','Color',[1 0 0],'FontWeight','bold','FontUnits','Pixel','FontSize',storage.display.slicefontsize,'HorizontalAlignment','right');
        end

        if storage.display.ticks == 0;
            axis off;
        end

        
    end
    
% render yz slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(orientation,'yz') == 1
    
    sliceavg = storage.average.x;
    
    if storage.position.yz + sliceavg > storage.dimensions.x
    		sliceavg = storage.dimensions.x - storage.position.yz;
    elseif storage.position.yz - sliceavg < 1
		sliceavg=1; %FIXME: calculate the correct number of slices
    end
    
    if storage.position.yz > 0 & storage.position.yz <= storage.dimensions.x

        %cut out slice and make average
        slice = single(storage.Volume.Value(round(storage.position.yz - sliceavg ./ 2):round(storage.position.yz + sliceavg ./ 2) - 1,:,:));
		slice = single(mean(slice, 1));
        slice = squeeze(slice);
        storage.slices.yz = slice;
        
        %bandpass
        if storage.filter.yz ~= 0
			slice = tom_bandpass(slice, storage.filter.low, storage.filter.high);
			DataScale = storage.DataScale ./ (storage.dimensions.y .* storage.dimensions.z);
            if storage.filter.low > 0
                DataScale = DataScale - storage.mean ./ (storage.dimensions.x .* storage.dimensions.y);
            end
        else
            DataScale = storage.DataScale;
        end
        
        %switchto yz axes
        axesobj = findobj('Tag','yz_axes');
        axes(axesobj); colormap(gray); axis ij;
        %create image
        imhandle = imagesc(slice,DataScale);
        axis(storage.actualaxis.yz);
        set(imhandle, 'Tag', 'yz_image');
        set(axesobj, 'Tag', 'yz_axes');
        set(axesobj,'YAxisLocation','right','Color','none');
        if storage.display.slicenumbers == 1
            delete(findobj('Tag','slicenumber_yz'));
            t = text('Units','normalized','Position',[.9 .1],'String',num2str(storage.position.yz));
        	set(t,'Tag','slicenumber_yz','Color',[1 0 0],'FontWeight','bold','FontUnits','Pixel','FontSize',storage.display.slicefontsize,'HorizontalAlignment','right');
        end
        
        if storage.display.ticks == 0;
            axis off;
        end
        
    end
else
    error('Unknown orientation parameter in render_slice()');
end


%Draw Position markers in all slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if storage.display.crosshair == 1

    %xy axes
    %%%%%%%%%%%%%%%
    axesobj = findobj('Tag','xy_axes');
    axes(axesobj);
    set(axesobj,'Units','pixel');
    delete(findobj('Tag','hline_xy'));
    delete(findobj('Tag','vline_xy'));
    line('Xdata', [0, storage.dimensions.x], 'Ydata', [storage.position.xz, storage.position.xz], 'Color', 'r', 'LineWidth',1, 'Tag', 'hline_xy');
    line('Xdata', [storage.position.yz, storage.position.yz], 'Ydata', [0, storage.dimensions.y], 'Color', 'r', 'LineWidth',1, 'Tag', 'vline_xy');
    set(axesobj,'Units','normalized');

    %xz axes
    %%%%%%%%%%%%%%%
    axesobj = findobj('Tag','xz_axes');
    axes(axesobj);
    set(axesobj,'Units','pixel');
    delete(findobj('Tag','hline_xz'));
    delete(findobj('Tag','vline_xz'));
    line('Xdata', [0, storage.dimensions.x], 'Ydata', [storage.position.xy, storage.position.xy], 'Color', 'r', 'LineWidth',1, 'Tag', 'hline_xz');
    line('Xdata', [storage.position.yz, storage.position.yz], 'Ydata', [0, storage.dimensions.z], 'Color', 'r', 'LineWidth',1, 'Tag', 'vline_xz');
    set(axesobj,'Units','normalized');

    %yz axes
    %%%%%%%%%%%%%%%
    axesobj = findobj('Tag','yz_axes');
    axes(axesobj);
    set(axesobj,'Units','pixel');
    delete(findobj('Tag','hline_yz'));
    delete(findobj('Tag','vline_yz'));
    line('Xdata', [storage.position.xy, storage.position.xy], 'Ydata', [0, storage.dimensions.y], 'Color', 'r', 'LineWidth',1, 'Tag', 'hline_yz');
    line('Xdata', [0 storage.dimensions.z], 'Ydata', [storage.position.xz, storage.position.xz], 'Color', 'r', 'LineWidth',1, 'Tag', 'vline_yz');
    set(axesobj,'Units','normalized');

end

%Draw subvolume box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if storage.subvol.data ~= 0
	x1 = storage.subvol.data(1);
	x2 = storage.subvol.data(2);
	y1 = storage.subvol.data(3);
	y2 = storage.subvol.data(4);
	z1 = storage.subvol.data(5);
	z2 = storage.subvol.data(6);
	axes(findobj('Tag','xy_axes'));
	rectangle('Position',[x1 y1 abs(x2 - x1) abs(y2 - y1)],'Edgecolor',[0 1 0], 'Linewidth',1,'Tag','subvol_rect');
	axes(findobj('Tag','xz_axes'));
	rectangle('Position',[x1 z1 abs(x2 - x1) abs(z2 - z1)],'Edgecolor',[0 1 0], 'Linewidth',1,'Tag','subvol_rect');
	axes(findobj('Tag','yz_axes'));
	rectangle('Position',[z1 y1 abs(z2 - z1) abs(y2 - y1)],'Edgecolor',[0 1 0], 'Linewidth',1,'Tag','subvol_rect');
end

draw_pickpoints(orientation);

set_callbacks();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  redraw function for 3D box                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_3dbox()

global storage;

boxobj = findobj('Tag','box');
axes(boxobj);

if storage.display.box == 1
    
    %newplot;
    cla;
    hold on;
    linewidth = [round(storage.dimensions.x ./ 64) round(storage.dimensions.x ./ 64)];

    %load xy slice
    slice_xy = storage.slices.xy';

    %load xz slice
    slice_xz = storage.slices.xz';

    %load yz slice
    slice_yz = storage.slices.yz;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %yz slice
    slice_yz = tom_mirror(slice_yz,'y');
    slice_yz = slice_yz';
    lauf = 0:5:storage.dimensions.y;
    y = [lauf;lauf];
    x = [zeros(1,size(y,2))+storage.position.yz;zeros(1,size(y,2))+storage.position.yz];
    z = [zeros(1,size(y,2))+1;zeros(1,size(y,2))+storage.dimensions.z];
    slice_yz(1:linewidth(1),:) = zeros()-1000;
    slice_yz(size(slice_yz,1)-linewidth(1):size(slice_yz,1),:)= zeros()-1000;
    slice_yz(:,size(slice_yz,2)-linewidth(1):size(slice_yz,2))= zeros()-1000;
    slice_yz(:,1:linewidth(1)) = zeros()-1000;
    warp(x,y,z,slice_yz);

    %xz slice
    lauf = 0:5:storage.dimensions.x;
    x = [lauf;lauf];
    y = [zeros(1,size(x,2))+storage.position.xz;zeros(1,size(x,2))+storage.position.xz];
    z = [zeros(1,size(x,2))+1;zeros(1,size(x,2))+storage.dimensions.z];

    pos_low = storage.position.yz-linewidth(2);
    if pos_low < 1
        pos_low = 1;
    end
    pos_high = storage.position.yz+linewidth(2);
    if pos_high > storage.dimensions.x
        pos_high = storage.dimensions.x
    end
    slice_xz = tom_mirror(slice_xz,'x');

    slice_xz(:,pos_low:pos_high) = zeros()-1000;
    slice_xz(1:linewidth(1),:) = zeros()-1000;
    slice_xz(size(slice_xz,1)-linewidth(1):size(slice_xz,1),:)= zeros()-1000;
    slice_xz(:,size(slice_xz,2)-linewidth(1):size(slice_xz,2))= zeros()-1000;
    slice_xz(:,1:linewidth(1)) = zeros()-1000;
    warp(x,y,z,slice_xz);

    %xy slice
    lauf = 0:5:storage.dimensions.x;
    x = [lauf;lauf];
    y = [zeros(1,size(x,2))+1;zeros(1,size(x,2))+storage.dimensions.y];
    z = [zeros(1,size(x,2))+storage.dimensions.z-storage.position.xy;zeros(1,size(x,2))+storage.dimensions.z-storage.position.xy];
    pos_low = storage.position.xz-linewidth(2);
    if pos_low < 1
        pos_low = 1;
    end
    pos_high = storage.position.xz+linewidth(2);
    if pos_high > storage.dimensions.y
        pos_high = storage.dimensions.y;
    end
    pos_low2 = storage.position.yz-linewidth(2);
    if pos_low2 < 1
        pos_low2 = 1;
    end
    pos_high2 = storage.position.yz+linewidth(2);
    if pos_high > storage.dimensions.x
        pos_high = storage.dimensions.x;
    end
    slice_xy(pos_low:pos_high,:) = zeros()-1000;
    slice_xy(:,pos_low2:pos_high2) = zeros()-1000;
    slice_xy(1:linewidth(1),:) = zeros()-1000;
    slice_xy(size(slice_xy,1)-linewidth(1):size(slice_xy,1),:)= zeros()-1000;
    slice_xy(:,size(slice_xy,2)-linewidth(1):size(slice_xy,2))= zeros()-1000;
    slice_xy(:,1:linewidth(1)) = zeros()-1000;
    warp(x,y,z,slice_xy);


    hold off;
    dx = storage.actualaxis.xy(1);
    dy = storage.actualaxis.xy(2);
    dz = storage.actualaxis.xz(2);
    axis([1,dx+5,1,dy+5,1,dz+5]);axis off;
    set(gca,'Projection','orthographic');
    caxis(storage.DataScale);
    axis tight;
    daspect([1,1,1]);
    set(boxobj,'Tag','box');
    
    
    %FIXME: This is not working correctly
    %if ~isempty(storage.measurement.point1) & ~isempty(storage.measurement.point2)
    %    line('Xdata',[storage.measurement.point1(1),storage.measurement.point2(1)],'YData',[storage.measurement.point1(2),storage.measurement.point2(2)],'Zdata',[storage.measurement.point1(3),storage.measurement.point2(3)],'LineWidth',2,'Tag','measureline','Color',[0 0 1]);
    %end
    %if ~isempty(storage.measurement.point2) & ~isempty(storage.measurement.point3)
    %    line('Xdata',[storage.measurement.point2(1),storage.measurement.point3(1)],'YData',[storage.measurement.point2(2),storage.measurement.point3(2)],'Zdata',[storage.measurement.point2(3),storage.measurement.point3(3)],'LineWidth',2,'Tag','measureline','Color',[0 0 1]);
    %end

else

    cla;
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  set callbacks for lines, images                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_callbacks()

set(findobj('Tag','xy_image'),'buttonDownFcn',@set_position);
set(findobj('Tag','xz_image'),'buttonDownFcn',@set_position);
set(findobj('Tag','yz_image'),'buttonDownFcn',@set_position);

set(findobj('Tag','vline_xy'),'buttonDownFcn',@set_position);
set(findobj('Tag','hline_xy'),'buttonDownFcn',@set_position);
set(findobj('Tag','vline_xz'),'buttonDownFcn',@set_position);
set(findobj('Tag','hline_xz'),'buttonDownFcn',@set_position);
set(findobj('Tag','vline_yz'),'buttonDownFcn',@set_position);
set(findobj('Tag','hline_yz'),'buttonDownFcn',@set_position);
set(findobj('Tag','pickpoint_xy'),'buttonDownFcn',@set_position);
set(findobj('Tag','pickpoint_xz'),'buttonDownFcn',@set_position);
set(findobj('Tag','pickpoint_yz'),'buttonDownFcn',@set_position);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  set callbacks for lines, images                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unset_callbacks(hObject, eventdata, handles)

set(findobj('Tag','xy_image'),'buttonDownFcn',@dummy);
set(findobj('Tag','xz_image'),'buttonDownFcn',@dummy);
set(findobj('Tag','yz_image'),'buttonDownFcn',@dummy);

set(findobj('Tag','vline_xy'),'buttonDownFcn',@dummy);
set(findobj('Tag','hline_xy'),'buttonDownFcn',@dummy);
set(findobj('Tag','vline_xz'),'buttonDownFcn',@dummy);
set(findobj('Tag','hline_xz'),'buttonDownFcn',@dummy);
set(findobj('Tag','vline_yz'),'buttonDownFcn',@dummy);
set(findobj('Tag','hline_yz'),'buttonDownFcn',@dummy);
set(findobj('Tag','pickpoint_xy'),'buttonDownFcn',@dummy);
set(findobj('Tag','pickpoint_xz'),'buttonDownFcn',@dummy);
set(findobj('Tag','pickpoint_yz'),'buttonDownFcn',@dummy);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  center the view                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function center_position()

global storage;

% adjust position sliders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpobj = findobj('Tag','navslider_x');
set(tmpobj,'Min',1);
set(tmpobj,'Max',storage.dimensions.x)
set(tmpobj,'SliderStep',[1 ./ storage.dimensions.x, 5 ./ storage.dimensions.x]);
set(tmpobj,'Value', floor(storage.dimensions.x ./ 2)+1);
set(findobj('Tag','nav_x'),'String',num2str(floor(storage.dimensions.x ./ 2)+1));

tmpobj = findobj('Tag','navslider_y');
set(tmpobj,'Min',1);
set(tmpobj,'Max',storage.dimensions.y)
set(tmpobj,'SliderStep',[1 ./ storage.dimensions.y, 5 ./ storage.dimensions.y]);
set(tmpobj,'Value',floor(storage.dimensions.y ./ 2)+1);
set(findobj('Tag','nav_y'),'String',num2str(floor(storage.dimensions.y ./ 2)+1));

tmpobj = findobj('Tag','navslider_z');
set(tmpobj,'Min',1);
set(tmpobj,'Max',storage.dimensions.z)
set(tmpobj,'SliderStep',[1 ./ storage.dimensions.z, 5 ./ storage.dimensions.z]);
set(tmpobj,'Value',floor(storage.dimensions.z ./ 2)+1);
set(findobj('Tag','nav_z'),'String',num2str(floor(storage.dimensions.z ./ 2)+1));

%set the position to the middle of the volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
storage.position.xy = floor(storage.dimensions.z ./ 2)+1;
storage.position.xz = floor(storage.dimensions.y ./ 2)+1;
storage.position.yz = floor(storage.dimensions.x ./ 2)+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  set_position callback function                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_position(a,b)

global storage;

point1=get(gca,'currentpoint');
orientation = get(gcbo,'Tag');
button = get(gcf,'selectiontype');

%button values:
%normal: left mouse button
%alt: right mouse button
%extend: middle mouse buttons

pt = point1(1,1:2);
x1 = round(pt(1));
y1 = round(pt(2));

%Handle Set Position event (left mouse button)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(button,'normal') == true

	if strcmp(orientation,'xy_image') == 1 || strcmp(orientation,'hline_xy') == 1 || strcmp(orientation,'vline_xy') == 1 || strcmp(orientation,'pickpoint_xy') == 1
		storage.position.xz = y1;
		storage.position.yz = x1;
	elseif strcmp(orientation,'xz_image') == 1 || strcmp(orientation,'hline_xz') == 1 || strcmp(orientation,'vline_xz') == 1 || strcmp(orientation,'pickpoint_xz') == 1
		storage.position.xy = y1;
		storage.position.yz = x1;
	elseif strcmp(orientation,'yz_image') == 1 || strcmp(orientation,'hline_yz') == 1 || strcmp(orientation,'vline_yz') == 1 || strcmp(orientation,'pickpoint_yz') == 1
		storage.position.xy = x1;
		storage.position.xz = y1;
	end

%Handle Set Marker Points event (right mouse button)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(button,'alt') == true
    
    pointnumber = size(storage.Align,2)+1;
    if strcmp(orientation,'xy_image') == 1 || strcmp(orientation,'hline_xy') == 1 || strcmp(orientation,'vline_xy') == 1
		xz = y1;
		yz = x1;
        xy = storage.position.xy;
	elseif strcmp(orientation,'xz_image') == 1 || strcmp(orientation,'hline_xz') == 1 || strcmp(orientation,'vline_xz') == 1
		xy = y1;
		yz = x1;
        xz = storage.position.xz;
	elseif strcmp(orientation,'yz_image') == 1 || strcmp(orientation,'hline_yz') == 1 || strcmp(orientation,'vline_yz') == 1
		xy = x1;
		xz = y1;
        yz = storage.position.yz;
    end
    storage.Align(1,pointnumber).Tomogram.Position.X = yz;
    storage.Align(1,pointnumber).Tomogram.Position.Y = xz;
    storage.Align(1,pointnumber).Tomogram.Position.Z = xy;
    
end

set(findobj('Tag','navslider_x'), 'Value', storage.position.yz);
set(findobj('Tag','nav_x'), 'String', num2str(storage.position.yz));
set(findobj('Tag','navslider_y'), 'Value', storage.position.xz);
set(findobj('Tag','nav_y'), 'String', num2str(storage.position.xz));
set(findobj('Tag','navslider_z'), 'Value', storage.position.xy);
set(findobj('Tag','nav_z'), 'String', num2str(storage.position.xy));

redraw_all_slices();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  drawmark paints a mark                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = drawmark(x,y,Color,Tag,Label)

global storage;

if nargin == 3 
    Tag = 'a';
    Label = '';
end
if nargin == 4 
    Label = '';
end

hold on;
Center= x + y*sqrt(-1);
%Radius = 5;
Radius = storage.dimensions.x .* 0.03;
%Gridpt = 100;
%[u,v]=circle(Center,Radius,Gridpt);
%line(u,v,'LineWidth',1,'Color',[1 0 0]);
uu = [x x x x-Radius x+Radius];
vv = [y-Radius y+Radius y y y];
h = line(uu,vv,'LineWidth',2,'color',Color,'Tag',Tag);
if ~isempty(Label)
    t = text(x+Radius.*1.3,y-Radius.*1.3,Label,'Color',Color,'Tag',Tag,'BackgroundColor',[1 1 1]);
end
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  renders all particle picking points                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_pickpoints(slice)

global storage;

%delete all existing points
tmpobj = findobj('Tag',['pickpoint_' slice]);

if tmpobj ~= 0
    delete(tmpobj);
end

if ~isfield(storage,'Align')
    return;
end

if strcmp(slice,'xy') == 1

    %find all points which have to be displayed in the xy image
    axesobj = findobj('Tag','xy_axes');
    axes(axesobj);
    for pointnumber=1:size(storage.Align,2)
        %exact match
        if storage.Align(1,pointnumber).Tomogram.Position.Z == storage.position.xy
            drawmark(storage.Align(1,pointnumber).Tomogram.Position.X,storage.Align(1,pointnumber).Tomogram.Position.Y,[0 1 0],'pickpoint_xy');
            %match near current position
        elseif sum(storage.Align(1,pointnumber).Tomogram.Position.Z == [storage.position.xy - storage.particleradius:storage.position.xy + storage.particleradius]) == 1
            drawmark(storage.Align(1,pointnumber).Tomogram.Position.X,storage.Align(1,pointnumber).Tomogram.Position.Y,[1 0 0],'pickpoint_xy');
        end
    end

elseif strcmp(slice,'xz')

    %find all points which have to be displayed in the xz image
    axesobj = findobj('Tag','xz_axes');
    axes(axesobj);
    for pointnumber=1:size(storage.Align,2)
        %exact match
        if storage.Align(1,pointnumber).Tomogram.Position.Y == storage.position.xz
            drawmark(storage.Align(1,pointnumber).Tomogram.Position.X,storage.Align(1,pointnumber).Tomogram.Position.Z,[0 1 0],'pickpoint_xz');
            %match near current position
        elseif sum(storage.Align(1,pointnumber).Tomogram.Position.Y == [storage.position.xz - storage.particleradius:storage.position.xz + storage.particleradius]) == 1
            drawmark(storage.Align(1,pointnumber).Tomogram.Position.X,storage.Align(1,pointnumber).Tomogram.Position.Z,[1 0 0],'pickpoint_xz');
        end
    end

elseif strcmp(slice,'yz')

    %find all points which have to be displayed in the yz image
    axesobj = findobj('Tag','yz_axes');
    axes(axesobj);
    for pointnumber=1:size(storage.Align,2)
        %exact match
        if storage.Align(1,pointnumber).Tomogram.Position.X == storage.position.yz
            drawmark(storage.Align(1,pointnumber).Tomogram.Position.Z,storage.Align(1,pointnumber).Tomogram.Position.Y,[0 1 0],'pickpoint_yz');
            %match near current position
        elseif sum(storage.Align(1,pointnumber).Tomogram.Position.X == [storage.position.yz - storage.particleradius:storage.position.yz + storage.particleradius]) == 1
            drawmark(storage.Align(1,pointnumber).Tomogram.Position.Z,storage.Align(1,pointnumber).Tomogram.Position.Y,[1 0 0],'pickpoint_yz');
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Recalculate image positions on window resize                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resize_window(hObject,eventdata,handles)

global storage;

figure(findobj('Tag','tom_volxyz_imagewindow'));
figuredimensions = get(findobj('Tag','tom_volxyz_imagewindow'),'Position');

%figuredimensions(3) = figuredimensions(4).*storage.figureaspect;

figurewidth = figuredimensions(3) - storage.display.outer_padding * 2 - storage.display.inner_padding;
figureheight = figuredimensions(4) - storage.display.outer_padding * 2 - storage.display.inner_padding;
yfrac = storage.dimensions.y ./ (storage.dimensions.z + storage.dimensions.y);
xfrac = storage.dimensions.x ./ (storage.dimensions.z + storage.dimensions.x);
zfrac = storage.dimensions.z ./ (storage.dimensions.z + storage.dimensions.y);

zsize = figureheight * zfrac;
xsize = figurewidth * xfrac;
ysize = figureheight * yfrac;

%Position of xy Slice
axesobj = findobj('Tag','xy_axes');
axes(axesobj);
left = storage.display.outer_padding;
bottom = storage.display.outer_padding + zsize + storage.display.inner_padding;
width = xsize;
height = ysize;
set(axesobj,'Units','pixel');
set(axesobj,'Tag','xy_axes');
set(axesobj, 'Position', [left bottom width height],'visible','on');
set(axesobj,'Units','normalized');
if storage.display.ticks == 0;
    axis off;
end


%Position of xz slice
axesobj = findobj('Tag','xz_axes');
axes(axesobj);
left = storage.display.outer_padding;
bottom = storage.display.outer_padding;
width = xsize;
height = zsize;
set(axesobj,'Units','pixel');
set(axesobj, 'Tag', 'xz_axes');
set(axesobj, 'Position', [left bottom width height],'visible','on');
set(axesobj,'Units','normalized');
if storage.display.ticks == 0;
    axis off;
end


%Position of yz slice
axesobj = findobj('Tag','yz_axes');
axes(axesobj);

left = storage.display.outer_padding + storage.display.inner_padding + xsize;
bottom = storage.display.outer_padding + zsize + storage.display.inner_padding;
width = zsize;
height = ysize;
set(axesobj,'Units','pixel');
set(axesobj,'Tag','yz_axes');
set(axesobj, 'Position', [left bottom width height],'visible','on');
set(axesobj,'Units','normalized');
set(axesobj,'YAxisLocation','right','Color','none');
set(axesobj,'Units','normalized');
if storage.display.ticks == 0;
    axis off;
end


%Position of 3D box
axesobj = findobj('Tag','box');
axes(axesobj);
left = storage.display.outer_padding + storage.display.inner_padding + xsize;
bottom = storage.display.outer_padding + 1;
width = zsize;
height = zsize;
set(axesobj,'Units','pixel');
set(axesobj,'Tag','box');
set(axesobj,'Position',[left bottom width height],'visible','on');
set(axesobj,'Units','normalized');axis off;
%set(findobj('Tag','tom_volxyz_imagewindow'),'Position',figuredimensions);

if get(findobj('Tag','display_slicenumbers'), 'Value') == 1
   
    axes(findobj('Tag','xy_axes'));
    delete(findobj('Tag','slicenumber_xy'));
    t = text(storage.actualaxis.xy(2)-(storage.actualaxis.xy(2)/100)*10,storage.actualaxis.xy(4)-(storage.actualaxis.xy(4)/100)*8,num2str(storage.position.xy));
    set(t,'Tag','slicenumber_xy','Color',[1 0 0],'FontWeight','bold','FontUnits','Normalized','FontSize',0.15,'HorizontalAlignment','right','Units','Normalized','Position',[.9 .1]);
    set(t,'FontUnits','Pixel');
    storage.display.slicefontsize = get(t,'FontSize');
        
    axes(findobj('Tag','xz_axes'));
    delete(findobj('Tag','slicenumber_xz'));
    t = text(storage.actualaxis.xz(2)-(storage.actualaxis.xz(2)/100)*10,storage.actualaxis.xz(4)-(storage.actualaxis.xz(4)/100)*8,num2str(storage.position.xz));
    set(t,'Tag','slicenumber_xz','Color',[1 0 0],'FontWeight','bold','FontUnits','Pixel','FontSize',storage.display.slicefontsize,'HorizontalAlignment','right','Units','Normalized','Position',[.9 .1]);
            
    axes(findobj('Tag','yz_axes'));
    delete(findobj('Tag','slicenumber_yz'));
    t = text(storage.actualaxis.yz(2)-(storage.actualaxis.yz(2)/100)*10,storage.actualaxis.yz(4)-(storage.actualaxis.yz(4)/100)*8,num2str(storage.position.yz));
   	set(t,'Tag','slicenumber_yz','Color',[1 0 0],'FontWeight','bold','FontUnits','Pixel','FontSize',storage.display.slicefontsize,'HorizontalAlignment','right','Units','Normalized','Position',[.9 .1]);	
    
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function histogram_low_Callback(hObject, eventdata, handles)
function histogram_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function histogram_high_Callback(hObject, eventdata, handles)
function histogram_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function bandpass_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function bandpass_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function navslider_x_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function navslider_y_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function navslider_z_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function nav_x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function nav_y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function nav_z_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function averageslider_x_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function averageslider_y_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function average_x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function averageslider_z_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function average_y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function average_z_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function scalebar_factor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rotate_phi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rotateslider_psi_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function rotate_psi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rotate_theta_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rotateslider_theta_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function rotateslider_phi_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function profile_thickness_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dummy(a,b)
function edit_threshold_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_threshold_Callback(hObject, eventdata, handles)
function checkbox_negstain_Callback(hObject, eventdata, handles)




