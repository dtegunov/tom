function varargout = viewer3d(varargin)
% VIEWER3D a Matlab 3D volume renderer using the fast shearwarp algorithm.
%
% VIEWER3D(V, RENDERTYPE, SCALES);
%
% inputs,
% V : 3D Input image volume, of type double, single, uint8, 
%            uint16, uint32, int8, int16 or int32 
%            (the render process uses only double calculations)
% RENDERTYPE: 'MIP' Maximum Intensity Render (default)
%             'VR' Volume Rendering
%             'VRC' Volume Rendering Color
%             'VRS' Volume Rendering with Shading
%             'SLICE', Slice viewing
% SCALES: The sizes(height, width, depth) of one voxel. (default [1 1 1])
%
% Volume Data, 
%  Range of V must be [0 1] in case of double or single. Volume Data of 
%  type double has shorter render times than data of uint8 or uint16.
%
% example,
%   % Load data
%   load TestData;
%   viewer3d(V);
%
% See also: render
%
% Function is written by D.Kroon University of Twente (March 2009)

% Edit the above text to modify the response to help viewer3d

% Last Modified by GUIDE v2.5 07-Jun-2010 11:23:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @viewer3d_OpeningFcn, ...
                   'gui_OutputFcn',  @viewer3d_OutputFcn, ...
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


% --- Executes just before viewer3d is made visible.
function viewer3d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to viewer3d (see VARARGIN)

% Choose default command line output for viewer3d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% addpath mexcode and help
functionname='viewer3d.m';
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
try addpath(functiondir); catch me, disp(me.message); end
try addpath([functiondir '/help']); catch me, disp(me.message); end


% Initialized data storage structure
data.mouse_pressed=false;
data.mouse_button='arrow';
data.autocontrast=false;

% Save the default config
filename_config=[functiondir '/default_config.mat'];
if(exist(filename_config,'file'))
    load(filename_config,'config')
    data.config=config;
else
    data.config.VolumeScaling=100;
    data.config.PreviewVolumeSize=32;
    data.config.ImageSizeRender=400;
    data.config.ShearInterpolation= 'bilinear';
    data.config.WarpInterpolation='bilinear';
    data.config.PreRender= 0;
    data.config.StoreXYZ=0;
    data.config.ColorSlice=false;
end

% Get input voxel volume and convert to double
if (isempty(varargin)), 
    data.volume_original=zeros(3,3,3);
    data.volume_preview=zeros(3,3,3);
else
    if(ndims(varargin{1})==3)
        data.volume_original=varargin{1};
        data=checkvolumetype(data);
    else
        error('viewer3d:inputs', 'Input image not 3 dimensional');
    end
end

% Make Preview, Render volume
data=makePreviewVolume(data);
data=makeRenderVolume(data);
% Some speedup data (if set by config)
data=makeVolumeXY(data);
data=computeNormals(data);

% Viewer vector and Light Vector
data.ViewerVector = [0 0 1];
data.LightVector = [0.5 -0.5 -0.67];

% Get input render type
if(length(varargin)>1)
    switch lower(varargin{2})
    case 'mip'
        data.render_type='mip';
    case 'vr'
        data.render_type='vr';
    case 'vrc'
        data.render_type='vrc';
    case 'vrs'
        data.render_type='vrs';
    case 'slice'
        data.render_type='slicez';
    otherwise
        error('viewer3d:inputs', 'Render type unknown');
    end
else
    data.render_type='mip';
end

% Get input voxelvolume scaling
if(length(varargin)>2)
    Scales=varargin{3}; 
    data.Scales=Scales;
else
    data.Scales=[1 1 1];
end

data.Zoom=(sqrt(3)./sqrt(sum(data.Scales.^2)));
data=set_initial_view_matrix(data);

% 2D initialization
data.SliceSelected=round(size(data.volume_original)/2);
data.handles=handles;
data.handle_viewer3d=gcf;
data.handle_histogram=[];
data.handle_console=[];
data.handle_voxelsize=[];
data.handle_lightvector=[];
data.handle_contrast=[];
data.handle_qualityspeed=[];
data.contrast=0; 
data.brightness=0; 
data.measure_distance=false;
data.measure_roi=false;
data.measure_landmark=false;
data.MeasureList=[];
data.tVolumemm=0;
data.VoxelLocation=[0 0 0];
data.histogram_positions = [0.2 0.4 0.6 0.9]; 
data.histogram_alpha = [0 0.5 0.35 1]; 
data.histogram_colors= [0 0 0; 1 0 0; 1 1 0; 1 1 1];
data.first_render=true;
data.histogram_pointselected=[];
data.mouse_position_pressed=[0 0];
data.mouse_position=[0 0];
data.mouse_position_last=[0 0];
data.shading_material='shiny';
data=loadmousepointershapes(data);
data.handles=handles;
set_menu_checks(data);
setMyData(data);
createAlphaColorTable();
% Show the data
show3d(false,true)

% UIWAIT makes viewer3d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function data=makePreviewVolume(data)
if(data.config.PreviewVolumeSize==100)
    data.volume_preview=data.volume_original;
else
    t=data.config.PreviewVolumeSize;
    data.volume_preview=imresize3d(data.volume_original,[],[t t t],'linear');
end

function data=makeRenderVolume(data)
if(data.config.VolumeScaling==100)
    data.volume=data.volume_original;
else
    data.volume=imresize3d(data.volume_original,data.config.VolumeScaling/100,[],'linear');
end

function createAlphaColorTable()
% This function creates a Matlab colormap and alphamap from the markers
data=getMyData(); if(isempty(data)), return, end
    data.colortable=zeros(1000,3); 
    data.alphatable=zeros(1000,1);
    % Loop through all 256 color/alpha indexes
    for j=0:999
        i=j/999;
        if (i<data.histogram_positions(1)), alpha=0; color=data.histogram_colors(1,:);
        elseif(i>data.histogram_positions(end)), alpha=0; color=data.histogram_colors(end,:);
        elseif(i==data.histogram_positions(1)), alpha=data.histogram_alpha(1); color=data.histogram_colors(1,:);
        elseif(i==data.histogram_positions(end)), alpha=data.histogram_alpha(end); color=data.histogram_colors(end,:);
        else
            % Linear interpolate the color and alpha between markers
            index_down=find(data.histogram_positions<=i); index_down=index_down(end);
            index_up=find(data.histogram_positions>i); index_up=index_up(1);
            perc=(i-data.histogram_positions(index_down)) / (data.histogram_positions(index_up) - data.histogram_positions(index_down));
            color=(1-perc)*data.histogram_colors(index_down,:)+perc*data.histogram_colors(index_up,:);
            alpha=(1-perc)*data.histogram_alpha(index_down)+perc*data.histogram_alpha(index_up);
        end
        data.colortable(j+1,:)=color;
        data.alphatable(j+1)=alpha;
    end
setMyData(data);

function data=loadmousepointershapes(data)
I=1-(imread('icon_mouse_rotate1.png')>0); I(I==0)=NaN;
data.icon_mouse_rotate1=I;
I=1-(imread('icon_mouse_rotate2.png')>0); I(I==0)=NaN;
data.icon_mouse_rotate2=I;
I=1-(imread('icon_mouse_zoom.png')>0); I(I==0)=NaN;
data.icon_mouse_zoom=I;
I=1-(imread('icon_mouse_pan.png')>0); I(I==0)=NaN;
data.icon_mouse_pan=I;

function show3d(preview,render_new_image)
data=getMyData(); if(isempty(data)), return, end
tic;
if(render_new_image)
    datarender=struct();
    datarender.ImageSize=[data.config.ImageSizeRender data.config.ImageSizeRender];
    switch data.render_type
    case 'mip'
        datarender.RenderType='mip';
        datarender.ShearInterp=data.config.ShearInterpolation;
        datarender.WarpInterp=data.config.WarpInterpolation;
    case 'vr'
        datarender.RenderType='bw';
        datarender.AlphaTable=data.alphatable;
        datarender.ShearInterp=data.config.ShearInterpolation;
        datarender.WarpInterp=data.config.WarpInterpolation;
    case 'vrc'
        datarender.RenderType='color';
        datarender.AlphaTable=data.alphatable;
        datarender.ColorTable=data.colortable;
        datarender.ShearInterp=data.config.ShearInterpolation;
        datarender.WarpInterp=data.config.WarpInterpolation;
    case 'vrs'
        datarender.RenderType='shaded';
        datarender.AlphaTable=data.alphatable; datarender.ColorTable=data.colortable;
        datarender.LightVector=data.LightVector; datarender.ViewerVector=data.ViewerVector;
        datarender.ShadingMaterial=data.shading_material;
        datarender.ShearInterp=data.config.ShearInterpolation;
        datarender.WarpInterp=data.config.WarpInterpolation;
    case 'slicex'
        datarender.RenderType='slicex';
        datarender.ColorTable=data.colortable;
        datarender.SliceSelected=data.SliceSelected(1);
        datarender.WarpInterp='bilinear';  
        datarender.ColorSlice=data.config.ColorSlice;
    case 'slicey'
        datarender.RenderType='slicey';
        datarender.ColorTable=data.colortable;
        datarender.SliceSelected=data.SliceSelected(2);
        datarender.WarpInterp='bilinear';  
        datarender.ColorSlice=data.config.ColorSlice;
    case 'slicez'
        datarender.RenderType='slicez';
        datarender.ColorTable=data.colortable;
        datarender.SliceSelected=data.SliceSelected(3);
        datarender.WarpInterp='bilinear';  
        datarender.ColorSlice=data.config.ColorSlice;
    end

    if(preview)
        switch data.render_type
            case {'slicex','slicey','slicez'}
                datarender.WarpInterp='nearest';
                datarender.Mview=data.viewer_matrix;
                data.render_image = render(data.volume_original, datarender);
            otherwise
                datarender.Mview=data.viewer_matrix*ResizeMatrix(size(data.volume_preview)./size(data.volume_original));
                data.render_image = render(data.volume_preview,datarender);
        end
    else
        mouse_button_old=data.mouse_button;
        set_mouse_shape('watch',data); drawnow('expose');
        switch data.render_type
            case {'slicex','slicey','slicez'}
                datarender.Mview=data.viewer_matrix;
                data.render_image = render(data.volume_original, datarender);
            otherwise
                datarender.Mview=data.viewer_matrix*ResizeMatrix(size(data.volume)./size(data.volume_original));
                datarender.VolumeX=data.volumex;
                datarender.VolumeY=data.volumey;
                datarender.Normals=data.normals;
                data.render_image = render(data.volume, datarender);
        end
        set_mouse_shape(mouse_button_old,data); drawnow('expose'); 
    end
end

% Add position information etc. to the rendered image
data.total_image=data.render_image;
if(data.autocontrast)
    data.total_image=data.total_image-min(data.total_image(:));
    data.total_image=data.total_image./(max(data.total_image(:))+eps);
end
if(data.contrast~=0||data.brightness~=0)
    data.total_image=data.total_image+data.brightness;
    data.total_image=data.total_image*10^data.contrast;
end

data=InfoOnScreen(data);
data=showMeasureList(data);

% To range
data.total_image(data.total_image<0)=0;
data.total_image(data.total_image>1)=1;

if(data.first_render)
    data.imshow_handle=imshow(data.total_image); drawnow('expose')
    data.first_render=false;
else
    set(data.imshow_handle,'Cdata',data.total_image);
end
data.axes_size=get(data.handles.axes3d,'PlotBoxAspectRatio');
set(get(data.handles.axes3d,'Children'),'ButtonDownFcn','viewer3d(''axes3d_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
data=console_addline(data,['Render Time : ' num2str(toc)]);

setMyData(data);

function data=set_initial_view_matrix(data)
switch data.render_type
    case 'slicex'
        data.viewer_matrix=[data.Scales(1)*data.Zoom 0 0 0; 0 data.Scales(2)*data.Zoom 0 0; 0 0 data.Scales(3)*data.Zoom 0; 0 0 0 1];
        data.viewer_matrix=[0 0 1 0;0 1 0 0; -1 0 0 0;0 0 0 1]*data.viewer_matrix;
    case 'slicey'
        data.viewer_matrix=[data.Scales(1)*data.Zoom 0 0 0; 0 data.Scales(2)*data.Zoom 0 0; 0 0 data.Scales(3)*data.Zoom 0; 0 0 0 1];
        data.viewer_matrix=[1 0 0 0;0 0 -1 0; 0 1 0 0;0 0 0 1]*data.viewer_matrix;
    case 'slicez'
        data.viewer_matrix=[data.Scales(1)*data.Zoom 0 0 0; 0 data.Scales(2)*data.Zoom 0 0; 0 0 data.Scales(3)*data.Zoom 0; 0 0 0 1];
        data.viewer_matrix=data.viewer_matrix*[1 0 0 0;0 1 0 0; 0 0 1 0;0 0 0 1];
    otherwise
        data.viewer_matrix=[data.Scales(1)*data.Zoom 0 0 0; 0 data.Scales(2)*data.Zoom 0 0; 0 0 data.Scales(3)*data.Zoom 0; 0 0 0 1];
end


% --- Outputs from this function are returned to the command line.
function varargout = viewer3d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function set_menu_checks(data)
set(data.handles.menu_render_mip,'Checked','off');
set(data.handles.menu_render_vr,'Checked','off');
set(data.handles.menu_render_vrc,'Checked','off');
set(data.handles.menu_render_vrs,'Checked','off');
set(data.handles.menu_render_xslice,'Checked','off');
set(data.handles.menu_render_yslice,'Checked','off');
set(data.handles.menu_render_zslice,'Checked','off');
set(data.handles.menu_measure,'Enable','off');
set(data.handles.menu_shiny,'Checked','off');
set(data.handles.menu_dull,'Checked','off');
set(data.handles.menu_metal,'Checked','off');

set(data.handles.menu_shiny,'Enable','off');
set(data.handles.menu_dull,'Enable','off');
set(data.handles.menu_metal,'Enable','off');
set(data.handles.menu_lightvector,'Enable','off');
set(data.handles.menu_config_slicescolor,'Enable','off');
set(data.handles.menu_config_slicescolor,'Checked','off');

switch data.render_type
case 'mip'
    set(data.handles.menu_render_mip,'Checked','on');
case 'vr'
    set(data.handles.menu_render_vr,'Checked','on');
case 'vrc'
    set(data.handles.menu_render_vrc,'Checked','on');
case 'vrs'
    set(data.handles.menu_shiny,'Enable','on');
    set(data.handles.menu_dull,'Enable','on');
    set(data.handles.menu_metal,'Enable','on');
    set(data.handles.menu_lightvector,'Enable','on');

    set(data.handles.menu_render_vrs,'Checked','on');
    switch data.shading_material
        case 'shiny'
            set(data.handles.menu_shiny,'Checked','on');
        case 'dull'
            set(data.handles.menu_dull,'Checked','on');
        case 'metal'
            set(data.handles.menu_metal,'Checked','on');
    end
case 'slicex'
    set(data.handles.menu_render_xslice,'Checked','on');
    set(data.handles.menu_measure,'Enable','on');
    set(data.handles.menu_config_slicescolor,'Enable','on');
    if(data.config.ColorSlice), set(data.handles.menu_config_slicescolor,'Checked','on'); end
case 'slicey'
    set(data.handles.menu_render_yslice,'Checked','on');
    set(data.handles.menu_measure,'Enable','on');
    set(data.handles.menu_config_slicescolor,'Enable','on');
    if(data.config.ColorSlice), set(data.handles.menu_config_slicescolor,'Checked','on'); end
case 'slicez'
    set(data.handles.menu_render_zslice,'Checked','on');
    set(data.handles.menu_measure,'Enable','on');
    set(data.handles.menu_config_slicescolor,'Enable','on');
    if(data.config.ColorSlice), set(data.handles.menu_config_slicescolor,'Checked','on'); end
end

% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cursor_position_in_axes(hObject,handles);
data=getMyData(); if(isempty(data)), return, end
if(isempty(data)), return, end;
if((length(data.render_type)>5)&&strcmp(data.render_type(1:5),'slice'))
    data=mouseposition_to_voxelposition(data);
    setMyData(data);
end

if(data.mouse_pressed)
    switch(data.mouse_button)
    case 'rotate1'
        r1=-360*(data.mouse_position_last(1)-data.mouse_position(1));
        r2=360*(data.mouse_position_last(2)-data.mouse_position(2));
        R=RotationMatrix([r1 r2 0]);
        data.viewer_matrix=R*data.viewer_matrix;
        setMyData(data);
        show3d(true,true)
    case 'rotate2'
        r1=100*(data.mouse_position_last(1)-data.mouse_position(1));
        r2=100*(data.mouse_position_last(2)-data.mouse_position(2));
        if(data.mouse_position(2)>0.5), r1=-r1; end
        if(data.mouse_position(1)<0.5), r2=-r2; end
        r3=r1+r2;
        R=RotationMatrix([0 0 r3]);
        data.viewer_matrix=R*data.viewer_matrix;
        setMyData(data);
        show3d(true,true)
    case 'pan'
        t2=200*(data.mouse_position_last(1)-data.mouse_position(1));
        t1=200*(data.mouse_position_last(2)-data.mouse_position(2));
        M=TranslateMatrix([t1 t2 0]);
        data.viewer_matrix=M*data.viewer_matrix;
        setMyData(data);
        show3d(true,true)      
    case 'zoom'
        z1=1+2*(data.mouse_position_last(1)-data.mouse_position(1));
        z2=1+2*(data.mouse_position_last(2)-data.mouse_position(2));
        z=0.5*(z1+z2); %sqrt(z1.^2+z2.^2);
        R=ResizeMatrix([z z z]); 
        data.Zoom=data.Zoom*(1/z);
        data.viewer_matrix=R*data.viewer_matrix;
        setMyData(data);
        show3d(true,true)        
    otherwise
    end
end

function R=RotationMatrix(r)
% Determine the rotation matrix (View matrix) for rotation angles xyz ...
    Rx=[1 0 0 0; 0 cosd(r(1)) -sind(r(1)) 0; 0 sind(r(1)) cosd(r(1)) 0; 0 0 0 1];
    Ry=[cosd(r(2)) 0 sind(r(2)) 0; 0 1 0 0; -sind(r(2)) 0 cosd(r(2)) 0; 0 0 0 1];
    Rz=[cosd(r(3)) -sind(r(3)) 0 0; sind(r(3)) cosd(r(3)) 0 0; 0 0 1 0; 0 0 0 1];
    R=Rx*Ry*Rz;
    
function M=ResizeMatrix(s)
	M=[1/s(1) 0 0 0;
	   0 1/s(2) 0 0;
	   0 0 1/s(3) 0;
	   0 0 0 1];

function M=TranslateMatrix(t)
	M=[1 0 0 -t(1);
	   0 1 0 -t(2);
	   0 0 1 -t(3);
	   0 0 0 1];
 
function cursor_position_in_axes(hObject,handles)
data=getMyData(); if(isempty(data)), return, end;
data.mouse_position_last=data.mouse_position;
% Get position of the mouse in the large axes
% p = get(0, 'PointerLocation');
% pf = get(hObject, 'pos');
% p(1:2) = p(1:2)-pf(1:2);
% set(gcf, 'CurrentPoint', p(1:2));
p = get(handles.axes3d, 'CurrentPoint');
data.mouse_position=[p(1, 1) p(1, 2)]./data.axes_size(1:2);
setMyData(data);

function setMyData(data)
% Store data struct in figure
setappdata(gcf,'data3d',data);

function data=getMyData()
% Get data struct stored in figure
data=getappdata(gcf,'data3d');

function VolumeROI=roi2binaryvolume(data)
VolumeROI=zeros(size(data.volume_original),'uint8');
if(~isfield(data,'MeasureList')), return; end

% if(data.render_type(6)==data.MeasureList(i).RenderSelected&&SliceSelected==data.MeasureList(i).SliceSelected)
for i=1:length(data.MeasureList)
    x=data.MeasureList(i).x;
    y=data.MeasureList(i).y;
    z=data.MeasureList(i).z;
    if(data.MeasureList(i).type=='r')
        switch (data.MeasureList(i).RenderSelected)
        case {'x'}
            J=squeeze(VolumeROI(x(1),:,:,:));
            J=bitmapplot(y,z,J,struct('FillColor',[1 1 1 1],'Color',[1 1 1 1]))>0;
            VolumeROI(x(1),:,:,:)=J;
        case {'y'}
            J=squeeze(VolumeROI(:,y(1),:,:));
            J=bitmapplot(x,z,J,struct('FillColor',[1 1 1 1],'Color',[1 1 1 1]))>0;
            VolumeROI(:,y(1),:,:)=J;
        case {'z'}
            J=squeeze(VolumeROI(:,:,z(1),:));
            J=bitmapplot(x,y,J,struct('FillColor',[1 1 1 1],'Color',[1 1 1 1]))>0;
            VolumeROI(:,:,z(1),:)=J;
        end
    end
end
    

          
        

% --- Executes on mouse press over axes background.
function axes3d_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.mouse_pressed=true;
data.mouse_button=get(handles.figure1,'SelectionType');
data.mouse_position_pressed=data.mouse_position;
if(strcmp(data.mouse_button,'normal'))
    if(data.measure_distance&&getnumberofpoints(data)==0)
        % Do measurement
        data=addMeasureList('p',data.VoxelLocation(1),data.VoxelLocation(2),data.VoxelLocation(3),0,data);
        data.mouse_button='select_distance';
        data.mouse_pressed=false;
        setMyData(data);
        show3d(false,false);
        return
     elseif(data.measure_distance&&getnumberofpoints(data)>0)
        VoxelLocation1=[data.MeasureList(end).x data.MeasureList(end).y data.MeasureList(end).z];
        VoxelLocation2=data.VoxelLocation;
        % First remove the point (will be replaced by distance)
        data=rmvMeasureList(data.MeasureList(end).id,data);
        % Do measurement
        x=[VoxelLocation1(1) VoxelLocation2(1)];
        y=[VoxelLocation1(2) VoxelLocation2(2)];
        z=[VoxelLocation1(3) VoxelLocation2(3)];
        dx=data.Scales(1)*(VoxelLocation1(1)-VoxelLocation2(1));
        dy=data.Scales(2)*(VoxelLocation1(2)-VoxelLocation2(2));
        dz=data.Scales(3)*(VoxelLocation1(3)-VoxelLocation2(3));
        distance=sqrt(dx.^2+dy.^2+dz.^2);
        data=addMeasureList('d',x,y,z,distance,data);
        data.measure_distance=false;
        data.mouse_pressed=false;
        setMyData(data);
        show3d(false,false);
        return
    end
    
    if(data.measure_landmark)
        data=addMeasureList('l',data.VoxelLocation(1),data.VoxelLocation(2),data.VoxelLocation(3),0,data);
        data.mouse_button='select_landmark';
        data.measure_landmark=false;
        data.mouse_pressed=false;
        setMyData(data);
        show3d(false,false);
        return
    end
    
    if(data.measure_roi)
        % Do measurement
        data=addMeasureList('p',data.VoxelLocation(1),data.VoxelLocation(2),data.VoxelLocation(3),0,data);
        data.mouse_button='select_roi';
        data.mouse_pressed=false;
        setMyData(data);
        show3d(false,false);
        return
    end
    
    distance_center=sum((data.mouse_position-[0.5 0.5]).^2);
    if((distance_center<0.15)&&data.render_type(1)~='s')
        data.mouse_button='rotate1';
        set_mouse_shape('rotate1',data)
    else
        data.mouse_button='rotate2';
        set_mouse_shape('rotate2',data)
    end
end
if(strcmp(data.mouse_button,'open'))
    if(data.measure_roi)
        data=addMeasureList('p',data.VoxelLocation(1),data.VoxelLocation(2),data.VoxelLocation(3),0,data);
        np=getnumberofpoints(data);
        x=zeros(1,np); y=zeros(1,np); z=zeros(1,np);
        
        for i=1:np;
            x(i)=data.MeasureList(end).x;
            y(i)=data.MeasureList(end).y;
            z(i)=data.MeasureList(end).z;
            data=rmvMeasureList(data.MeasureList(end).id,data);
        end
        % Close the loop
        x=[x x(1)]; y=[y y(1)]; z=[z z(1)];
      
        % Generate spline interpolated contour points
          [points] = fnplt(cscvn([x;y;z]));
          x = round(squeeze(points(1,1:end-1))); 
          y = round(squeeze(points(2,1:end-1)));
          z = round(squeeze(points(3,1:end-1)));
   
                  
        switch (data.render_type)
        case {'slicex'}
            x_2d=y; y_2d=z; 
            sizeI=[size(data.volume_original,2) size(data.volume_original,3)];
        case {'slicey'}
            x_2d=x; y_2d=z; 
            sizeI=[size(data.volume_original,1) size(data.volume_original,3)];
        case {'slicez'}
            x_2d=x; y_2d=y; 
            sizeI=[size(data.volume_original,1) size(data.volume_original,2)];
        end

        I=bitmapplot(x_2d,y_2d,zeros(sizeI),struct('FillColor',[1 1 1 1],'Color',[1 1 1 1]))>0;
        volume=sum(I(:))*prod(data.Scales);
        data=addMeasureList('r',x,y,z,volume,data);
        data.measure_roi=false;
        data.mouse_pressed=false;
        setMyData(data);
        show3d(false,false);
        return            
    end
end

if(strcmp(data.mouse_button,'extend'))
    data.mouse_button='pan';
    set_mouse_shape('pan',data)
end
if(strcmp(data.mouse_button,'alt'))
    if(data.render_type(1)=='s')
        % Get the mouse position
        x_2d=data.mouse_position(2);
        y_2d=data.mouse_position(1); 

        % To rendered image position
        x_2d=round(x_2d*data.config.ImageSizeRender); 
        y_2d=round(y_2d*data.config.ImageSizeRender);
        x_2d_start=x_2d-3; x_2d_start(x_2d_start<1)=1;
        x_2d_end=x_2d+3; x_2d_end(x_2d_end>size(data.hitmap,1))=size(data.hitmap,1);
        y_2d_start=y_2d-3; y_2d_start(y_2d_start<1)=1;
        y_2d_end=y_2d+3; y_2d_end(y_2d_end>size(data.hitmap,2))=size(data.hitmap,2);
        hitmap_part=data.hitmap(x_2d_start:x_2d_end,y_2d_start:y_2d_end);
        id_detect=max(hitmap_part(:));
        if(id_detect>0)
            data=rmvMeasureList(id_detect,data);
            setMyData(data);
            show3d(false,false);
            return;
        end
    end
        data.mouse_button='zoom';
        set_mouse_shape('zoom',data);
end
setMyData(data);

function p=getnumberofpoints(data)
p=0;
for i=length(data.MeasureList):-1:1,
    if(data.MeasureList(i).type=='p'), p=p+1; else return; end
end

function data=showMeasureList(data)
    data.hitmap=zeros([size(data.total_image,1) size(data.total_image,2)]);
    if(~isfield(data,'MeasureList')), return; end
    if(length(data.render_type)<6), return; end
    SliceSelected=data.SliceSelected(uint8(data.render_type(6))-119);
    for i=1:length(data.MeasureList)
        if(data.render_type(6)==data.MeasureList(i).RenderSelected&&SliceSelected==data.MeasureList(i).SliceSelected)
            id=data.MeasureList(i).id;
            x=data.MeasureList(i).x;
            y=data.MeasureList(i).y;
            z=data.MeasureList(i).z;

            [x,y]=voxelposition_to_imageposition(x,y,z,data);
            
            if(data.MeasureList(i).type=='d')
                distancemm=data.MeasureList(i).varmm;
                data=plotDistance(x,y,distancemm,id,data);
            elseif(data.MeasureList(i).type=='r')
                volumemm=data.MeasureList(i).varmm;
                data=plotRoi(x,y,volumemm,id,data);
            elseif(data.MeasureList(i).type=='p')
                data=plotPoint(x,y,id,data);
            elseif(data.MeasureList(i).type=='l')
                data=plotPointBlue(x,y,id,data);
            end
        end
    end
    
function data=plotPoint(x,y,id,data)
I=data.total_image;
I=bitmapplot(x,y,I,struct('Marker','*','MarkerColor',[1 0 0 1],'Color',[0 0 1 1]));
data.hitmap=bitmapplot(x,y,data.hitmap,struct('Marker','*','MarkerColor',[id id id 1],'Color',[id id id 1]));
data.total_image=I;

function data=plotPointBlue(x,y,id,data)
I=data.total_image;
I=bitmapplot(x,y,I,struct('Marker','*','MarkerColor',[0 0 1 1],'Color',[1 0 0 1]));
data.hitmap=bitmapplot(x,y,data.hitmap,struct('Marker','*','MarkerColor',[id id id 1],'Color',[id id id 1]));
data.total_image=I;

function data=plotDistance(x,y,distancemm,id,data)
I=data.total_image;
I=bitmapplot(x,y,I,struct('Marker','*','MarkerColor',[1 0 0 1],'Color',[0 0 1 1]));
data.hitmap=bitmapplot(x,y,data.hitmap,struct('Color',[id id id 1]));
info=[num3str(distancemm,0,2) ' mm']; 
I=bitmaptext(info,I,[mean(x)-5 mean(y)-5],struct('Color',[0 1 0 1]));
data.total_image=I;

function data=plotRoi(x,y,volumemm,id,data)
I=data.total_image;
I=bitmapplot(x,y,I,struct('FillColor',[0 0 1 0.1],'Color',[1 0 0 1]));
data.hitmap=bitmapplot(x,y,data.hitmap,struct('Color',[id id id 1]));
info=[num3str(volumemm,0,2) ' mm^3']; 
I=bitmaptext(info,I,[mean(x)-5 mean(y)-5],struct('Color',[0 1 0 1]));
data.total_image=I;

function numstr=num3str(num,bef,aft)
    numstr=num2str(num,['%.' num2str(aft) 'f']);
    if(aft>0), maxlen = bef + aft +1; else maxlen = bef; end
    while (length(numstr)<maxlen), numstr=['0' numstr]; end
    

function data=addMeasureList(type,x,y,z,varmm,data)
    if(isempty(data.MeasureList)), p=1; else p=length(data.MeasureList)+1; end
    data.MeasureList(p).id=(rand+sum(clock)-floor(sum(clock)))/2;
    data.MeasureList(p).type=type;
    data.MeasureList(p).RenderSelected=data.render_type(6);
    SliceSelected=data.SliceSelected(uint8(data.render_type(6))-119);
    data.MeasureList(p).SliceSelected=SliceSelected;
    data.MeasureList(p).x=x;
    data.MeasureList(p).y=y;
    data.MeasureList(p).z=z;
    data.MeasureList(p).varmm=varmm;
    data=calcTotalVolume(data);

function data=rmvMeasureList(id,data)
    index=-1;
    for i=1:length(data.MeasureList)
        if(data.MeasureList(i).id==id), index=i;end
    end
    if(index>-1)
        data.MeasureList(index)=[];
    end
   data=calcTotalVolume(data);

function data=calcTotalVolume(data)
    data.tVolumemm=0;
    for i=1:length(data.MeasureList)
        if(data.MeasureList(i).type=='r')
            data.tVolumemm=data.tVolumemm+data.MeasureList(i).varmm;
        end
    end
   
    
function set_mouse_shape(type,data)

switch(type)
case 'rotate1'
    set(gcf,'Pointer','custom','PointerShapeCData',data.icon_mouse_rotate1,'PointerShapeHotSpot',round(size(data.icon_mouse_rotate1)/2))
    set(data.handles.figure1,'Pointer','custom');
case 'rotate2'
    set(gcf,'Pointer','custom','PointerShapeCData',data.icon_mouse_rotate2,'PointerShapeHotSpot',round(size(data.icon_mouse_rotate2)/2))
    set(data.handles.figure1,'Pointer','custom');
case 'select_distance'
    set(data.handles.figure1,'Pointer','crosshair')
case 'select_landmark'
    set(data.handles.figure1,'Pointer','crosshair')
case 'select_roi'
    set(data.handles.figure1,'Pointer','crosshair')
case 'normal'
    set(data.handles.figure1,'Pointer','arrow')
case 'alt'
    set(data.handles.figure1,'Pointer','arrow')
case 'open'
    set(data.handles.figure1,'Pointer','arrow')
case 'zoom'
    set(gcf,'Pointer','custom','PointerShapeCData',data.icon_mouse_zoom,'PointerShapeHotSpot',round(size(data.icon_mouse_zoom)/2))
    set(data.handles.figure1,'Pointer','custom');
case 'pan'
    set(gcf,'Pointer','custom','PointerShapeCData',data.icon_mouse_pan,'PointerShapeHotSpot',round(size(data.icon_mouse_pan)/2))
    set(data.handles.figure1,'Pointer','custom');
    
otherwise
    set(data.handles.figure1,'Pointer',type);
end
  
% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
if(data.mouse_pressed)
    data.mouse_pressed=false;
    setMyData(data);
    show3d(false,true)
end
if(data.measure_distance)
    set_mouse_shape('select_distance',data)
elseif(data.measure_roi)
    set_mouse_shape('select_roi',data)
else
    set_mouse_shape('arrow',data)
end


function A=imresize3d(V,scale,tsize,ntype,npad)
% This function resizes a 3D image volume to new dimensions
% Vnew = imresize3d(V,scale,nsize,ntype,npad);
%
% inputs,
%   V: The input image volume
%   scale: scaling factor, when used set tsize to [];
%   nsize: new dimensions, when used set scale to [];
%   ntype: Type of interpolation ('nearest', 'linear', or 'cubic')
%   npad: Boundary condition ('replicate', 'symmetric', 'circular', 'fill', or 'bound')  
%
% outputs,
%   Vnew: The resized image volume
%
% example,
%   load('mri','D'); D=squeeze(D);
%   Dnew = imresize3d(D,[],[80 80 40],'nearest','bound');
%
% This function is written by D.Kroon University of Twente (July 2008)
% Check the inputs
if(exist('ntype', 'var') == 0), ntype='nearest'; end
if(exist('npad', 'var') == 0), npad='bound'; end
if(exist('scale', 'var')&&~isempty(scale)), tsize=round(size(V)*scale); end
if(exist('tsize', 'var')&&~isempty(tsize)),  scale=(tsize./size(V)); end

% Make transformation structure   
T = makehgtform('scale',scale);
tform = maketform('affine', T);

% Specify resampler
R = makeresampler(ntype, npad);

% Resize the image volueme
A = tformarray(V, tform, R, [1 2 3], [1 2 3], tsize, [], 0);

% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_config_Callback(hObject, eventdata, handles)
% hObject    handle to menu_config (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_change_alpha_colors_Callback(hObject, eventdata, handles)
% hObject    handle to menu_change_alpha_colors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.handle_histogram=viewer3d_histogram(data.handle_viewer3d);
handles_histogram=guidata(data.handle_histogram);
data.handle_histogram_axes=handles_histogram.axes_histogram;
setMyData(data);
createHistogram();
drawHistogramPoints();

% --------------------------------------------------------------------
function menu_load_view_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataold=getMyData();
dataold.volume=[];
if(ishandle(dataold.handle_histogram)), close(dataold.handle_histogram); end
uiload();
if(exist('data','var'))
    data.first_render=true;
    
    % The handles must be from the current figure
    data.handle_viewer3d=dataold.handle_viewer3d;
    data.handles.axes3d=dataold.handles.axes3d;
    data.handles.figure1=dataold.handles.figure1;
    data.handle_console=dataold.handle_console;
    data.handle_contrast=dataold.handle_contrast;
    data.handle_voxelsize=dataold.handle_voxelsize;
    data.handle_lightvector=dataold.handle_lightvector;
    data.handles=handles;
    
    % Add fields which are missing in the stored data, because of an older
    % version or something like that.
    tags = fieldnames(dataold);
    for i=1:length(tags)
         if(~isfield(data,tags{i})),  
             data.(tags{i})=dataold.(tags{i}); 
         else
            if(isstruct(data.(tags{i}))&&isstruct(dataold.(tags{i})))
                    tags2 = fieldnames(dataold.(tags{i}));
                    for j=1:length(tags2)
                        if(~isfield(data.(tags{i}),tags2{j})),  
                            data.(tags{i}).(tags2{j})=dataold.(tags{i}).(tags2{j}); 
                        end
                    end
            end
         end
    end

    data=makeVolumeXY(data);
    data=computeNormals(data);
    data=makePreviewVolume(data);
    data=makeRenderVolume(data);
    
    setMyData(data);
    set_menu_checks(data);
    createAlphaColorTable();
    show3d(false,true);
else
    viewer3d_error({'Matlab File does not contain','data from "Save View"'})
end

function load_variable_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
if(ishandle(data.handle_histogram)), close(data.handle_histogram); end
data.volume_original=eventdata;
data.SliceSelected=round(size(data.volume_original)/2);
data.Scales=[1 1 1];
data.first_render=true;
data.handles=handles;
data.Zoom=(sqrt(3)./sqrt(sum(data.Scales.^2)));
data=set_initial_view_matrix(data);

data=makeVolumeXY(data);
data=computeNormals(data);
data=makePreviewVolume(data);
data=makeRenderVolume(data);

setMyData(data);
set_menu_checks(data);
show3d(false,true);

% --------------------------------------------------------------------
function menu_save_view_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.volume_preview=[];
data.volume=[];
data.volumex=[];
data.volumey=[];
data.normals=[];
uisave('data');

function menu_load_histogram_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
alpha=0;
uiload;
if(exist('positions','var'))
    data.histogram_positions=positions;
    data.histogram_colors=colors;
    data.histogram_alpha=alpha;
    setMyData(data);
    drawHistogramPoints();
    createAlphaColorTable();
    show3d(false,true);
else
    viewer3d_error({'Matlab File does not contain','data from "Save AlphaColors"'})
end
% --------------------------------------------------------------------
function menu_save_histogram_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
positions=data.histogram_positions;
colors=data.histogram_colors;
alpha=data.histogram_alpha;
uisave({'positions','colors','alpha'});

% --------------------------------------------------------------------
function menu_render_Callback(hObject, eventdata, handles)
% hObject    handle to menu_render (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_render_mip_Callback(hObject, eventdata, handles)
% hObject    handle to menu_render_mip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.render_type='mip';
set_menu_checks(data);
data.first_render=true;
setMyData(data);
show3d(false,true);

% --------------------------------------------------------------------
function menu_render_vr_Callback(hObject, eventdata, handles)
% hObject    handle to menu_render_vr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.render_type='vr';
set_menu_checks(data);
data.first_render=true;
setMyData(data);
show3d(false,true);

% --------------------------------------------------------------------
function menu_render_vrc_Callback(hObject, eventdata, handles)
% hObject    handle to menu_render_vrc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.render_type='vrc';
set_menu_checks(data);
data.first_render=true;
setMyData(data);
show3d(false,true);

% --------------------------------------------------------------------
function menu_render_vrs_Callback(hObject, eventdata, handles)
% hObject    handle to menu_render_vrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.render_type='vrs';
set_menu_checks(data);
data.first_render=true;
setMyData(data);
show3d(false,true);

% --------------------------------------------------------------------
function menu_info_Callback(hObject, eventdata, handles)
% hObject    handle to menu_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_save_picture_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_picture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile({'*.png';'*.jpg'}, 'Save Rendered Image as');
data=getMyData(); if(isempty(data)), return, end
imwrite(data.total_image,[pathname filename]);

% --------------------------------------------------------------------
function menu_about_Callback(hObject, eventdata, handles)
% hObject    handle to menu_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
viewer3d_about
   
function createHistogram()
% This function creates and show the (log) histogram of the data    
data=getMyData(); if(isempty(data)), return, end
    % Get histogram
    switch(class(data.volume_original))
    case 'int8'
        volumepart=uint8(single(data.volume(1:8:end))+2^7);
    case 'int16'
        volumepart=uint16(single(data.volume(1:8:end))+2^15);
    case 'int32';
        volumepart=uint32(single(data.volume(1:8:end))+2^31);
    otherwise
        volumepart=data.volume(:);
    end
     
    [data.histogram_countsy, data.histogram_countsx]=imhist(volumepart);
    % Log the histogram data
    data.histogram_countsy=log(data.histogram_countsy+100); data.histogram_countsy=data.histogram_countsy-min(data.histogram_countsy);
    data.histogram_countsx=data.histogram_countsx./max(data.histogram_countsx(:));
    data.histogram_countsy=data.histogram_countsy./max(data.histogram_countsy(:));
    % Focus on histogram axes
    figure(data.handle_histogram)    
    % Display the histogram
    stem(data.handle_histogram_axes,data.histogram_countsx,data.histogram_countsy,'Marker', 'none'); 
    hold(data.handle_histogram_axes,'on'); 
    % Set the axis of the histogram axes
    data.histogram_maxy=max(data.histogram_countsy(:));
    data.histogram_maxx=max(data.histogram_countsx(:));
    
    set(data.handle_histogram_axes,'yLim', [0 1]);
    set(data.handle_histogram_axes,'xLim', [0 1]);
setMyData(data);

% --- Executes on selection change in popupmenu_colors.
function popupmenu_colors_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_colors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_colors contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_colors
data=getMyData(); if(isempty(data)), return, end
    % Generate the new color markers
    c_choice=get(handles.popupmenu_colors,'Value');
    ncolors=length(data.histogram_positions);
    switch c_choice,
        case 1,new_colormap=jet(1000); 
        case 2, new_colormap=hsv(1000);
        case 3, new_colormap=hot(1000);
        case 4, new_colormap=cool(1000);
        case 5, new_colormap=spring(1000);
        case 6, new_colormap=summer(1000);
        case 7, new_colormap=autumn(1000);
        case 8, new_colormap=winter(1000);
        case 9, new_colormap=gray(1000);
        case 10, new_colormap=bone(1000);
        case 11, new_colormap=copper(1000);
        case 12, new_colormap=pink(1000);
        otherwise, new_colormap=hot(1000);
    end
    new_colormap=new_colormap(round(1:(end-1)/(ncolors-1):end),:);
    data.histogram_colors=new_colormap;
    
    % Draw the new color markers and make the color and alpha map
    setMyData(data);
    drawHistogramPoints();
    createAlphaColorTable();
    show3d(false,true);

function drawHistogramPoints()
data=getMyData(); if(isempty(data)), return, end
    % Delete old points and line
    try
        delete(data.histogram_linehandle), 
        for i=1:length(data.histogram_pointhandle), 
           delete(data.histogram_pointhandle(i)), 
        end, 
    catch
    end
    stem(data.handle_histogram_axes,data.histogram_countsx,data.histogram_countsy,'Marker', 'none'); 
    hold(data.handle_histogram_axes,'on');
    
    % Display the markers and line through the markers.
    data.histogram_linehandle=plot(data.handle_histogram_axes,data.histogram_positions,data.histogram_alpha*data.histogram_maxy,'m');
    set(data.histogram_linehandle,'ButtonDownFcn','viewer3d(''lineHistogramButtonDownFcn'',gcbo,[],guidata(gcbo))');
    for i=1:length(data.histogram_positions)
        data.histogram_pointhandle(i)=plot(data.handle_histogram_axes,data.histogram_positions(i),data.histogram_alpha(i)*data.histogram_maxy,'bo','MarkerFaceColor',data.histogram_colors(i,:));
        set(data.histogram_pointhandle(i),'ButtonDownFcn','viewer3d(''pointHistogramButtonDownFcn'',gcbo,[],guidata(gcbo))');
    end
    
    % For detection of mouse up, down and motion in histogram figure.
    set(data.handle_histogram, 'WindowButtonDownFcn','viewer3d(''HistogramButtonDownFcn'',gcbo,[],guidata(gcbo))');
    set(data.handle_histogram, 'WindowButtonMotionFcn','viewer3d(''HistogramButtonMotionFcn'',gcbo,[],guidata(gcbo))');
    set(data.handle_histogram, 'WindowButtonUpFcn','viewer3d(''HistogramButtonUpFcn'',gcbo,[],guidata(gcbo))');
setMyData(data);    

function pointHistogramButtonDownFcn(hObject, eventdata, handles)
data=getMyData(); if(isempty(data)), return, end
data.mouse_button=get(data.handle_histogram,'SelectionType');
if(strcmp(data.mouse_button,'normal'))
    data.histogram_pointselected=find(data.histogram_pointhandle==gcbo);
    data.histogram_pointselectedhandle=gcbo;
    set(data.histogram_pointselectedhandle, 'MarkerSize',8);
    setMyData(data);
elseif(strcmp(data.mouse_button,'extend'))
    data.histogram_pointselected=find(data.histogram_pointhandle==gcbo);
    data.histogram_colors(data.histogram_pointselected,:)=rand(1,3);
    data.histogram_pointselected=[];
    setMyData(data);
    drawHistogramPoints();
    createAlphaColorTable();    
    % Show the data
    histogram_handles=guidata(data.handle_histogram);
    if(get(histogram_handles.checkbox_auto_update,'value'))
        show3d(false,true);
    else
        show3d(true,true);    
    end

elseif(strcmp(data.mouse_button,'alt'))
    data.histogram_pointselected=find(data.histogram_pointhandle==gcbo);

    data.histogram_positions(data.histogram_pointselected)=[];
    data.histogram_colors(data.histogram_pointselected,:)=[];
    data.histogram_alpha(data.histogram_pointselected)=[];

    data.histogram_pointselected=[];
    setMyData(data);
    drawHistogramPoints();
    createAlphaColorTable();
    % Show the data
    histogram_handles=guidata(data.handle_histogram);
    if(get(histogram_handles.checkbox_auto_update,'value'))
        show3d(false,true);
    else
        show3d(true,true);    
    end
end

function HistogramButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function HistogramButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
if(~isempty(data.histogram_pointselected))
    set(data.histogram_pointselectedhandle, 'MarkerSize',6);
    data.histogram_pointselected=[];
    setMyData(data);
    createAlphaColorTable();
    % Show the data
    histogram_handles=guidata(data.handle_histogram);
    if(get(histogram_handles.checkbox_auto_update,'value'))
        show3d(false,true);
    else
        show3d(true,true);
    end
 end

function Histogram_pushbutton_update_view_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_update_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
show3d(false,true)

function HistogramButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cursor_position_in_histogram_axes(hObject,handles);
data=getMyData(); if(isempty(data)), return, end
if(~isempty(data.histogram_pointselected))
 % Set point to location mouse
        data.histogram_positions(data.histogram_pointselected)=data.histogram_mouse_position(1,1); 
        data.histogram_alpha(data.histogram_pointselected)=data.histogram_mouse_position(1,2);
        
        % Correct new location
        if(data.histogram_alpha(data.histogram_pointselected)<0), data.histogram_alpha(data.histogram_pointselected)=0; end
        if(data.histogram_alpha(data.histogram_pointselected)>1), data.histogram_alpha(data.histogram_pointselected)=1; end
        if(data.histogram_positions(data.histogram_pointselected)<0), data.histogram_positions(data.histogram_pointselected)=0; end
        if(data.histogram_positions(data.histogram_pointselected)>1), data.histogram_positions(data.histogram_pointselected)=1; end
        if((data.histogram_pointselected>1)&&(data.histogram_positions(data.histogram_pointselected-1)>data.histogram_positions(data.histogram_pointselected)))
            data.histogram_positions(data.histogram_pointselected)=data.histogram_positions(data.histogram_pointselected-1);
        end
        
        if((data.histogram_pointselected<length(data.histogram_positions))&&(data.histogram_positions(data.histogram_pointselected+1)<data.histogram_positions(data.histogram_pointselected)))
            data.histogram_positions(data.histogram_pointselected)=data.histogram_positions(data.histogram_pointselected+1);
        end

        % Move point
        set(data.histogram_pointselectedhandle, 'xdata', data.histogram_positions(data.histogram_pointselected));
        set(data.histogram_pointselectedhandle, 'ydata', data.histogram_alpha(data.histogram_pointselected));
        
        % Move line
        set(data.histogram_linehandle, 'xdata',data.histogram_positions);
        set(data.histogram_linehandle, 'ydata',data.histogram_alpha);
end
setMyData(data);

function lineHistogramButtonDownFcn(hObject, eventdata, handles)
data=getMyData(); if(isempty(data)), return, end
        % New point on mouse location
        newposition=data.histogram_mouse_position(1,1);
        
        % List for the new markers
        newpositions=zeros(1,length(data.histogram_positions)+1);
        newalphas=zeros(1,length(data.histogram_alpha)+1);
        newcolors=zeros(size(data.histogram_colors,1)+1,3);

        % Check if the new point is between old points
        index_down=find(data.histogram_positions<=newposition); 
        if(isempty(index_down)) 
        else
            index_down=index_down(end);
            index_up=find(data.histogram_positions>newposition); 
            if(isempty(index_up)) 
            else
                index_up=index_up(1);
                
                % Copy the (first) old markers to the new lists
                newpositions(1:index_down)=data.histogram_positions(1:index_down);
                newalphas(1:index_down)=data.histogram_alpha(1:index_down);
                newcolors(1:index_down,:)=data.histogram_colors(1:index_down,:);
                
                % Add the new interpolated marker
                perc=(newposition-data.histogram_positions(index_down)) / (data.histogram_positions(index_up) - data.histogram_positions(index_down));
                color=(1-perc)*data.histogram_colors(index_down,:)+perc*data.histogram_colors(index_up,:);
                alpha=(1-perc)*data.histogram_alpha(index_down)+perc*data.histogram_alpha(index_up);
                
                newpositions(index_up)=newposition; 
                newalphas(index_up)=alpha; 
                newcolors(index_up,:)=color;
              
                % Copy the (last) old markers to the new lists
                newpositions(index_up+1:end)=data.histogram_positions(index_up:end);
                newalphas(index_up+1:end)=data.histogram_alpha(index_up:end);
                newcolors(index_up+1:end,:)=data.histogram_colors(index_up:end,:);
        
                % Make the new lists the used marker lists
                data.histogram_positions=newpositions; 
                data.histogram_alpha=newalphas; 
                data.histogram_colors=newcolors;
            end
        end
        
        % Update the histogram window
        cla(data.handle_histogram_axes);
setMyData(data);
drawHistogramPoints();
createAlphaColorTable();
% Show the data
histogram_handles=guidata(data.handle_histogram);
if(get(histogram_handles.checkbox_auto_update,'value'))
    show3d(false,true);
else
    show3d(true,true);    
end
       
function cursor_position_in_histogram_axes(hObject,handles)
data=getMyData(); if(isempty(data)), return, end
%     % Get position of the mouse in the large axes
%     p = get(0, 'PointerLocation');
%     pf = get(hObject, 'pos');
%     p(1:2) = p(1:2)-pf(1:2);
%     set(data.handle_histogram, 'CurrentPoint', p(1:2));
    p = get(data.handle_histogram_axes, 'CurrentPoint');
    data.histogram_mouse_position=[p(1, 1) p(1, 2)];
setMyData(data);

% --------------------------------------------------------------------
function menu_help_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('info.html');

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
data=getMyData(); if(isempty(data)), delete(hObject); return, end
try 
if(ishandle(data.handle_histogram)), delete(data.handle_histogram); end
if(ishandle(data.handle_qualityspeed)), delete(data.handle_qualityspeed); end
if(ishandle(data.handle_console)), delete(data.handle_console);  end
if(ishandle(data.handle_contrast)), delete(data.handle_contrast); end
if(ishandle(data.handle_voxelsize)), delete(data.handle_voxelsize); end
if(ishandle(data.handle_lightvector)),  delete(data.handle_lightvector); end
catch me
    disp(me.message);
end

% Remove the data of this figure
try rmappdata(gcf,'data3d'); catch end
delete(hObject);


% --------------------------------------------------------------------
function menu_shiny_Callback(hObject, eventdata, handles)
% hObject    handle to menu_shiny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.shading_material='shiny';
set_menu_checks(data);
setMyData(data);
show3d(false,true);

% --------------------------------------------------------------------
function menu_dull_Callback(hObject, eventdata, handles)
% hObject    handle to menu_dull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.shading_material='dull';
set_menu_checks(data);
setMyData(data);
show3d(false,true);

% --------------------------------------------------------------------
function menu_metal_Callback(hObject, eventdata, handles)
% hObject    handle to menu_metal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.shading_material='metal';
set_menu_checks(data);
setMyData(data);
show3d(false,true);

function button_voxelsize_apply_Callback(hObject, eventdata, handles)
data=getMyData(); if(isempty(data)), return, end
handles_voxelsize=guidata(data.handle_voxelsize);

Scales_old=data.Scales;
Zoom_old=data.Zoom;

data.Scales(1)=str2double(get(handles_voxelsize.edit_scax,'String'));
data.Scales(2)=str2double(get(handles_voxelsize.edit_scay,'String'));
data.Scales(3)=str2double(get(handles_voxelsize.edit_scaz,'String'));
data.first_render=true;
data.Zoom=(sqrt(3)./sqrt(sum(data.Scales.^2)));

data.viewer_matrix=data.viewer_matrix*ResizeMatrix((Scales_old.*Zoom_old)./(data.Scales.*data.Zoom));
%data=set_initial_view_matrix(data);

setMyData(data);
show3d(false,true);



% --------------------------------------------------------------------
function menu_voxelsize_Callback(hObject, eventdata, handles)
% hObject    handle to menu_voxelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.handle_voxelsize=viewer3d_voxelsize;
setMyData(data);
handles_voxelsize=guidata(data.handle_voxelsize);
set(handles_voxelsize.edit_volx,'String',num2str(size(data.volume,1)));
set(handles_voxelsize.edit_voly,'String',num2str(size(data.volume,2)));
set(handles_voxelsize.edit_volz,'String',num2str(size(data.volume,3)));
set(handles_voxelsize.edit_scax,'String',num2str(data.Scales(1)));
set(handles_voxelsize.edit_scay,'String',num2str(data.Scales(2)));
set(handles_voxelsize.edit_scaz,'String',num2str(data.Scales(3)));


function button_lightvector_apply_Callback(hObject, eventdata, handles)
data=getMyData(); if(isempty(data)), return, end
handles_lightvector=guidata(data.handle_lightvector);
data.LightVector(1)=str2double(get(handles_lightvector.edit_lightx,'String'));
data.LightVector(2)=str2double(get(handles_lightvector.edit_lighty,'String'));
data.LightVector(3)=str2double(get(handles_lightvector.edit_lightz,'String'));
setMyData(data);
show3d(false,true);

% --------------------------------------------------------------------
function menu_lightvector_Callback(hObject, eventdata, handles)
% hObject    handle to menu_voxelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.handle_lightvector=viewer3d_lightvector;
setMyData(data);
handles_lightvector=guidata(data.handle_lightvector);
set(handles_lightvector.edit_lightx,'String',num2str(data.LightVector(1)));
set(handles_lightvector.edit_lighty,'String',num2str(data.LightVector(2)));
set(handles_lightvector.edit_lightz,'String',num2str(data.LightVector(3)));


% --------------------------------------------------------------------
function menu_load_worksp_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_worksp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.handle_workspacevars=viewer3d_workspacevars;
setMyData(data)
% Get variables in the workspace
vars = evalin('base','who');
% Select only variables with 3 dimensions
vars3d=[];
for i=1:length(vars),
    if(evalin('base',['ndims(' vars{i} ')'])==3), vars3d{length(vars3d)+1}=vars{i}; end
end
% Show the 3D variables in the workspace
handles_workspacevars=guidata(data.handle_workspacevars);
set(handles_workspacevars.listbox_vars,'String',vars3d);

% --- Executes on button press in pushbutton1.
function workspacevars_button_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
handles_workspacevars=guidata(data.handle_workspacevars);
list_entries = get(handles_workspacevars.listbox_vars,'String');
index_selected = get(handles_workspacevars.listbox_vars,'Value');
if length(index_selected) ~= 1
	errordlg('You must select one variable')
    return;
else
    var1 = list_entries{index_selected(1)};
    evalin('base',['viewer3d(''load_variable_Callback'',gcf,' var1 ',guidata(gcf))']);
end 
if(ishandle(data.handle_workspacevars)), close(data.handle_workspacevars); end


function console_button_clear_Callback(hObject, eventdata, handles)
data=getMyData(); if(isempty(data)), return, end
data.consoletext=[];
data.consolelines=0;
set(data.handle_console_edit,'String','');
setMyData(data);

function data=console_addline(data,newline)
if(ishandle(data.handle_console)),
    data.consolelines=data.consolelines+1;
    data.consoletext{data.consolelines}=newline;
    if(data.consolelines>14), 
        data.consolelines=14; 
        data.consoletext={data.consoletext{2:end}}; 
    end
    set(data.handle_console_edit,'String',data.consoletext);
end

% --------------------------------------------------------------------
function menu_console_Callback(hObject, eventdata, handles)
% hObject    handle to menu_console (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.handle_console=viewer3d_console;
handles_console=guidata(data.handle_console);
data.handle_console_edit=handles_console.edit_console;
data.consoletext=[];
data.consolelines=0;
setMyData(data)
set(data.handle_console_edit,'String','');


% --------------------------------------------------------------------
function menu_quality_speed_Callback(hObject, eventdata, handles)
% hObject    handle to menu_quality_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.handle_qualityspeed=viewer3d_qualityspeed;
setMyData(data);
handles_qualityspeed=guidata(data.handle_qualityspeed);
switch(data.config.VolumeScaling)
    case{25}
        set(handles_qualityspeed.uipanel_VolumeScaling,'SelectedObject',handles_qualityspeed.radiobutton_scaling25);
    case{50}
        set(handles_qualityspeed.uipanel_VolumeScaling,'SelectedObject',handles_qualityspeed.radiobutton_scaling50);
    case{100}
        set(handles_qualityspeed.uipanel_VolumeScaling,'SelectedObject',handles_qualityspeed.radiobutton_scaling100);
    case{200}
        set(handles_qualityspeed.uipanel_VolumeScaling,'SelectedObject',handles_qualityspeed.radiobutton_scaling200);
end
switch(data.config.PreviewVolumeSize)
    case{32}
        set(handles_qualityspeed.uipanel_PreviewVolumeSize,'SelectedObject',handles_qualityspeed.radiobutton_preview_32);
    case{64}
        set(handles_qualityspeed.uipanel_PreviewVolumeSize,'SelectedObject',handles_qualityspeed.radiobutton_preview_64);
    case{100}
        set(handles_qualityspeed.uipanel_PreviewVolumeSize,'SelectedObject',handles_qualityspeed.radiobutton_preview_100);
end
switch(data.config.ImageSizeRender)
    case{150}
        set(handles_qualityspeed.uipanel_ImageSizeRender,'SelectedObject',handles_qualityspeed.radiobutton_rendersize150);
    case{250}
        set(handles_qualityspeed.uipanel_ImageSizeRender,'SelectedObject',handles_qualityspeed.radiobutton_rendersize250);
    case{400}
        set(handles_qualityspeed.uipanel_ImageSizeRender,'SelectedObject',handles_qualityspeed.radiobutton_rendersize400);
    case{600}
        set(handles_qualityspeed.uipanel_ImageSizeRender,'SelectedObject',handles_qualityspeed.radiobutton_rendersize600);
    case{800}
        set(handles_qualityspeed.uipanel_ImageSizeRender,'SelectedObject',handles_qualityspeed.radiobutton_rendersize800);
    case{1400}
        set(handles_qualityspeed.uipanel_ImageSizeRender,'SelectedObject',handles_qualityspeed.radiobutton_rendersize1400);
    case{2500}
        set(handles_qualityspeed.uipanel_ImageSizeRender,'SelectedObject',handles_qualityspeed.radiobutton_rendersize2500);
end
switch(data.config.ShearInterpolation)
    case{'bilinear'}
        set(handles_qualityspeed.uipanel_ShearInterpolation,'SelectedObject',handles_qualityspeed.radiobutton_shear_int_bilinear);
    case{'nearest'}
        set(handles_qualityspeed.uipanel_ShearInterpolation,'SelectedObject',handles_qualityspeed.radiobutton_shear_int_nearest);
end
switch(data.config.WarpInterpolation)
    case{'bilinear'}
        set(handles_qualityspeed.uipanel_WarpInterpolation,'SelectedObject',handles_qualityspeed.radiobutton_warp_int_bilinear);
    case{'nearest'}
        set(handles_qualityspeed.uipanel_WarpInterpolation,'SelectedObject',handles_qualityspeed.radiobutton_warp_int_nearest);
end
set(handles_qualityspeed.checkbox_prerender,'Value',data.config.PreRender);
set(handles_qualityspeed.checkbox_storexyz,'Value',data.config.StoreXYZ);


% --- Executes on button press in pushbutton_applyconfig.
function qualityspeed_pushbutton_applyconfig_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_applyconfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
handles_qualityspeed=guidata(data.handle_qualityspeed);

VolumeScaling=get(get(handles_qualityspeed.uipanel_VolumeScaling,'SelectedObject'),'Tag');
PreviewVolumeSize=get(get(handles_qualityspeed.uipanel_PreviewVolumeSize,'SelectedObject'),'Tag');
ShearInterpolation=get(get(handles_qualityspeed.uipanel_ShearInterpolation,'SelectedObject'),'Tag');
WarpInterpolation=get(get(handles_qualityspeed.uipanel_WarpInterpolation,'SelectedObject'),'Tag');
ImageSizeRender=get(get(handles_qualityspeed.uipanel_ImageSizeRender,'SelectedObject'),'Tag');

VolumeScaling=str2double(VolumeScaling(20:end));
ImageSizeRender=str2double(ImageSizeRender(23:end));

PreviewVolumeSize=str2double(PreviewVolumeSize(21:end));
data.config.ShearInterpolation=ShearInterpolation(23:end);
data.config.WarpInterpolation=WarpInterpolation(22:end);
data.config.PreRender=get(handles_qualityspeed.checkbox_prerender,'Value');
data.config.StoreXYZ=get(handles_qualityspeed.checkbox_storexyz,'Value');

if(ImageSizeRender~=data.config.ImageSizeRender)
    s=data.config.ImageSizeRender/ImageSizeRender;
    data.viewer_matrix=data.viewer_matrix*ResizeMatrix([s s s]);
    data.config.ImageSizeRender=ImageSizeRender;
end

scale_change=data.config.VolumeScaling~=VolumeScaling;
if(scale_change)
      data.config.VolumeScaling=VolumeScaling;
      data=makeRenderVolume(data);
end

if(data.config.PreviewVolumeSize~=PreviewVolumeSize)
    data.config.PreviewVolumeSize=PreviewVolumeSize;
    data=makePreviewVolume(data);
end
if((isempty(data.volumey)||scale_change)&&data.config.StoreXYZ)
    data=makeVolumeXY(data);
end
if(~data.config.StoreXYZ)
    data.volumex=[]; data.volumey=[];
end

if((isempty(data.normals)||scale_change)&&data.config.PreRender)
    % Make normals
    data=computeNormals(data);
end
if(~data.config.PreRender)
    data.normals=[];
end

data.first_render=true;
setMyData(data);
show3d(false,true);

function data=makeVolumeXY(data)
if(data.config.StoreXYZ)
    data.volumex=shiftdim(data.volume,1);
    data.volumey=shiftdim(data.volume,2);
else
    data.volumex=[];
    data.volumey=[];
end

function data=computeNormals(data)
if(data.config.PreRender)
    % Pre computer Normals for faster shading rendering.
    [fy,fx,fz]=gradient(imgaussian(double(data.volume),1/2));
    flength=sqrt(fx.^2+fy.^2+fz.^2)+1e-6;
    data.normals=zeros([size(fx) 3]);
    data.normals(:,:,:,1)=fx./flength;
    data.normals(:,:,:,2)=fy./flength;
    data.normals(:,:,:,3)=fz./flength;
else
    data.normals=[];
end

function I=imgaussian(I,sigma,siz)
% IMGAUSSIAN filters an 1D, 2D or 3D image with an gaussian filter.
% This function uses IMFILTER, for the filtering but instead of using
% a multidimensional gaussian kernel, it uses the fact that a gaussian
% filter can be separated in 1D gaussian kernels.
%
% J=IMGAUSSIAN(I,SIGMA,SIZE)
%
% inputs,
%   I: The 1D, 2D, or 3D input image
%   SIGMA: The sigma used for the gaussian
%   SIZE: Kernel size (single value) (default: sigma*6)
% 
% outputs,
%   J: The gaussian filterd image
%
% example,
%   I = im2double(rgb2gray(imread('peppers.png')));
%   figure, imshow(imgaussian(I,3));
% 
% Function is written by D.Kroon University of Twente (October 2008)

if(~exist('siz','var')), siz=sigma*6; end

% Make 1D gaussian kernel
x=-(siz/2)+0.5:siz/2;
H = exp(-(x.^2/(2*sigma^2))); 
H = H/sum(H(:));

% Filter each dimension with the 1D gaussian kernels
if(ndims(I)==1)
    I=imfilter(I,H);
elseif(ndims(I)==2)
    Hx=reshape(H,[length(H) 1]); 
    Hy=reshape(H,[1 length(H)]); 
    I=imfilter(imfilter(I,Hx),Hy);
elseif(ndims(I)==3)
    Hx=reshape(H,[length(H) 1 1]); 
    Hy=reshape(H,[1 length(H) 1]); 
    Hz=reshape(H,[1 1 length(H)]);
    I=imfilter(imfilter(imfilter(I,Hx),Hy),Hz);
else
    error('imgaussian:input','unsupported input dimension');
end

             
% --- Executes on button press in pushbutton_saveconfig.
function qualityspeed_pushbutton_saveconfig_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_saveconfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
handles_qualityspeed=guidata(data.handle_qualityspeed);
VolumeScaling=get(get(handles_qualityspeed.uipanel_VolumeScaling,'SelectedObject'),'Tag');
PreviewVolumeSize=get(get(handles_qualityspeed.uipanel_PreviewVolumeSize,'SelectedObject'),'Tag');
ShearInterpolation=get(get(handles_qualityspeed.uipanel_ShearInterpolation,'SelectedObject'),'Tag');
WarpInterpolation=get(get(handles_qualityspeed.uipanel_WarpInterpolation,'SelectedObject'),'Tag');
ImageSizeRender=get(get(handles_qualityspeed.uipanel_ImageSizeRender,'SelectedObject'),'Tag');
VolumeScaling=str2double(VolumeScaling(20:end));
PreviewVolumeSize=str2double(PreviewVolumeSize(21:end));
data.config.ImageSizeRender=str2double(ImageSizeRender(23:end));
data.config.ShearInterpolation=ShearInterpolation(23:end);
data.config.WarpInterpolation=WarpInterpolation(22:end);
data.config.PreRender=get(handles_qualityspeed.checkbox_prerender,'Value');
data.config.StoreXYZ=get(handles_qualityspeed.checkbox_storexyz,'Value');
data.config.VolumeScaling=VolumeScaling;
data.config.PreviewVolumeSize=PreviewVolumeSize;

% Save the default config
config=data.config;
functiondir=which('viewer3d.m'); functiondir=functiondir(1:end-length('viewer3d.m'));
save([functiondir '/default_config.mat'],'config')


% --------------------------------------------------------------------
function menu_measure_Callback(hObject, eventdata, handles)
% hObject    handle to menu_measure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_config_contrast_Callback(hObject, eventdata, handles)
% hObject    handle to menu_config_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.handle_contrast=viewer3d_contrast(data.handle_viewer3d);
handles_contrast=guidata(data.handle_contrast);
set(handles_contrast.slider_contrast,'value',data.contrast);
set(handles_contrast.slider_brightness,'value',data.brightness);
set(handles_contrast.checkbox_autocontrast,'value',data.autocontrast);
setMyData(data);


% --------------------------------------------------------------------
function menu_measure_distance_Callback(hObject, eventdata, handles)
% hObject    handle to menu_measure_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.mouse_button='select_distance';
data.measure_distance=true;
setMyData(data);
set_mouse_shape('select_distance',data)
    

% --------------------------------------------------------------------
function menu_measure_roi_Callback(hObject, eventdata, handles)
% hObject    handle to menu_measure_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.mouse_button='select_roi';
data.measure_roi=true;
setMyData(data);
set_mouse_shape('select_roi',data)

% --------------------------------------------------------------------
function menu_render_xslice_Callback(hObject, eventdata, handles)
% hObject    handle to menu_render_xslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.render_type='slicex';
data=set_initial_view_matrix(data);
set_menu_checks(data);
data.first_render=true;
setMyData(data);
show3d(false,true);

% --------------------------------------------------------------------
function menu_render_yslice_Callback(hObject, eventdata, handles)
% hObject    handle to menu_render_yslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.render_type='slicey';
data=set_initial_view_matrix(data);
set_menu_checks(data);
data.first_render=true;
setMyData(data);
show3d(false,true);

% --------------------------------------------------------------------
function menu_render_zslice_Callback(hObject, eventdata, handles)
% hObject    handle to menu_render_zslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.render_type='slicez';
data=set_initial_view_matrix(data);
set_menu_checks(data);
data.first_render=true;
setMyData(data);
show3d(false,true);

% --------------------------------------------------------------------
function menu_load_dicom_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_dicom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
if(ishandle(data.handle_histogram)), close(data.handle_histogram); end
[volume,scales] = viewer3d_dicom;
if(exist('volume','var'))
    data.volume_original=volume;
    data=checkvolumetype(data);
    data.SliceSelected=round(size(data.volume_original)/2);
    data.Scales=scales;
    data.first_render=true;
    data.handles=handles;
    data.Zoom=(sqrt(3)./sqrt(sum(data.Scales.^2)));
    data.viewer_matrix=[data.Scales(1)*data.Zoom 0 0 0; 0 data.Scales(2)*data.Zoom 0 0; 0 0 data.Scales(3)*data.Zoom 0; 0 0 0 1];
    
    data=makeVolumeXY(data);
    data=computeNormals(data);
    data=makePreviewVolume(data);
    data=makeRenderVolume(data);
    
    setMyData(data);
    set_menu_checks(data);
    show3d(false,true);
else
    viewer3d_error({'Matlab Dicom Data Load Error'})
end

function data=checkvolumetype(data)
     switch(class(data.volume_original))
     case 'uint8'
         if(max(data.volume_original(:))<2^7),  data.autocontrast=true; end
     case 'uint16'
         if(max(data.volume_original(:))<2^15), data.autocontrast=true; end
     case 'uint32'
         if(max(data.volume_original(:))<2^31), data.autocontrast=true; end
     case 'int8'
         if(max(data.volume_original(:))<2^6),  data.autocontrast=true; end
     case 'int16',
         if(max(data.volume_original(:))<2^14), data.autocontrast=true; end
     case 'int32',
         if(max(data.volume_original(:))<2^30), data.autocontrast=true; end
     case 'single'
         if((min(data.volume_original(:))<0)||(max(data.volume_original(:))>1)), 
               warning('viewer3d:inputs', 'Input not using valid data range [0..1]');    
         end
         data.volume_original(data.volume_original<0)=0; 
         data.volume_original(data.volume_original>1)=1;
     case 'double'
         if((min(data.volume_original(:))<0)||(max(data.volume_original(:))>1)), 
               warning('viewer3d:inputs', 'Input not using valid data range [0..1]');    
         end
        data.volume_original(data.volume_original<0)=0; 
        data.volume_original(data.volume_original>1)=1;
     otherwise
        warning('viewer3d:inputs', 'Unsupported input datatype converted to double');
        data.volume_original=double(data.volume_original);
        % To range [0 1]
        data.volume_original=data.volume_original-min(data.volume_original(:));
        data.volume_original=data.volume_original./max(data.volume_original(:));
     end

     

% --------------------------------------------------------------------
function menu_save_rois_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
VolumeROI=roi2binaryvolume(data);
uisave('VolumeROI');

% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
switch data.render_type
    case {'slicex','slicey','slicez'}
        handles=guidata(hObject);
        data=changeslice(eventdata.VerticalScrollCount,handles,data);
        setMyData(data);
        show3d(false,true);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
if(strcmp(eventdata.Key,'uparrow')), eventdata.Character='+'; end
if(strcmp(eventdata.Key,'downarrow')), eventdata.Character='-'; end
    
switch data.render_type
    case {'slicex','slicey','slicez'}
        handles=guidata(hObject);
        switch(eventdata.Character)
            case '+'
                data=changeslice(1,handles,data);
                setMyData(data); show3d(true,true);
            case '-'
                data=changeslice(-1,handles,data);
                setMyData(data); show3d(true,true);
            case 'r'
                menu_measure_roi_Callback(hObject, eventdata, handles);
            case 'd'
                menu_measure_distance_Callback(hObject, eventdata, handles);
            case 'l'
                menu_measure_landmark_Callback(hObject, eventdata, handles);
            otherwise
        end
     otherwise        
end

function data=changeslice(updown,handles,data)
switch data.render_type
case 'slicex'
    data.SliceSelected(1)=data.SliceSelected(1)+updown;
    if(data.SliceSelected(1)>size(data.volume_original,1)),  data.SliceSelected(1)=size(data.volume_original,1); end
case 'slicey'
    data.SliceSelected(2)=data.SliceSelected(2)+updown;
    if(data.SliceSelected(2)>size(data.volume_original,2)),  data.SliceSelected(2)=size(data.volume_original,2); end
case 'slicez'
    data.SliceSelected(3)=data.SliceSelected(3)+updown;
    if(data.SliceSelected(3)>size(data.volume_original,3)),  data.SliceSelected(3)=size(data.volume_original,3); end
end
% Boundary limit
data.SliceSelected(data.SliceSelected<1)=1;
% Stop measurement
data.measure_distance=false;
data.measure_roi=false;

% --- Executes on key release with focus on figure1 and none of its controls.
function figure1_KeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
switch data.render_type
    case {'slicex','slicey','slicez'}
        show3d(false,true);
end

function data=InfoOnScreen(data)
if (size(data.total_image,3)==3)
    I=data.total_image;
else
    % Greyscale to color
    I(:,:,1)=data.total_image; 
    I(:,:,2)=data.total_image; 
    I(:,:,3)=data.total_image;
end

if(data.render_type(1)=='s')
    info=[];
    switch data.render_type
    case 'slicex'
        info{1}=['Slice X : ' num2str(data.SliceSelected(1))];
    case 'slicey'
        info{1}=['Slice y : ' num2str(data.SliceSelected(2))];
    case 'slicez'
        info{1}=['Slice Z : ' num2str(data.SliceSelected(3))];
    end
    info{2}=['ROIs mm^3: ' num2str(data.tVolumemm)];
    info{3}=['x,y,z px: ' num2str(data.VoxelLocation(1)) ' - ' num2str(data.VoxelLocation(2)) ' - ' num2str(data.VoxelLocation(3))]; 
    info{4}=['x,y,z mm: ' num2str(data.VoxelLocation(1)*data.Scales(1)) ' - ' num2str(data.VoxelLocation(2)*data.Scales(2)) ' - ' num2str(data.VoxelLocation(3)*data.Scales(3))]; 
    
    I=bitmaptext(info,I,[1 1],struct('Color',[0 1 0 1]));
end

data.total_image=I;



function I=bitmapplot(x,y,Ibackground,options)
% BITMAPPLOT, Linear plot in bitmap.
%
% Iplot=bitmapplot(x,y,Ibackground,options)
%
% inputs,
%   x : a vector with x values
%   y : a vector with y values, with same length as x 
%   Ibackground: the bitmap used as background when a m x n x 3 matrix
%       color plots are made, when m x n a greyscale plot.
%   options: struct with options such as color
% 
% outputs,
%   Iplot: The bitmap containing the plotted lines
%
% note,
%   Colors are always [r(ed) g(reen) b(lue) a(pha)], with range 0..1.
%   when Ibackground is grayscale, the mean of r,g,b is used as grey value.
%
% options,
%   options.Color: The color of the line.
%   options.FillColor: If this color is set, the region between 
%          the x and y coordnates will be filled with this color.
%   options.LineWidth: Thickness of the line in pixels 1,2,3..n
%   options.Marker: The marker type: 'o', '+' or '*'.
%   options.MarkerColor: The color of the markers used.
%   options.MarkerSize: The size of the markers used
%
% example,
%   % Make empty bitmap
%   I = zeros([320 256 3]);
%   
%   % Add a line
%   x=rand(1,10)*50+50; y=linspace(1,512,10);
%   I=bitmapplot(x,y,I);
%
%   % Add a thick red line
%   x=rand(1,10)*50+100; y=linspace(1,256,10);
%   I=bitmapplot(x,y,I,struct('LineWidth',5,'Color',[1 0 0 1]));
%
%   % Add a line with markers
%   x=rand(1,10)*50+150; y=linspace(1,256,10);
%   I=bitmapplot(x,y,I,struct('Marker','*','MarkerColor',[1 0 1 1],'Color',[1 1 0 1]));
%
%   % Add a filled polygon
%   x=[1 100 30 100]+200; y=[30 1 250 200];
%   I=bitmapplot(x,y,I,struct('FillColor',[0 1 0 0.5],'Color',[1 1 0 1]));
%
%   % Add a filled polygon on top
%   x=[30 80 70 120]+200; y=[30 1 250 200];
%   I=bitmapplot(x,y,I,struct('FillColor',[1 0 0 0.5],'Color',[1 0 0 1]));
%
%   lines={'Plot Test,','BitmapPlot version 1.0'};
%   % Plot text into background image
%   I=bitmaptext(lines,I,[1 1],struct('Color',[1 1 1 1]));
%   
%   % Show the bitmap
%   figure, imshow(I);
%
% Function is written by D.Kroon University of Twente (March 2009)


% Process inputs
defaultoptions=struct('Color',[0 0 1 1],'FillColor',[],'LineWidth',1,'Grid',[],'MarkerColor',[1 0 0 1],'Marker',[],'MarkerSize',6);
if(~exist('options','var')), 
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(options))), 
        warning('register_images:unknownoption','unknown options found');
    end
end

% The function works with double values (store class for ouput)
Classb=class(Ibackground);
Ibackground=im2double(Ibackground);

% x and y to row vectors
x=x(:)'; y=y(:)';
x=round(x); y=round(y);

% Make line, marker an fill bitmap
I_line=zeros([size(Ibackground,1) size(Ibackground,2)]);
I_marker=zeros([size(Ibackground,1) size(Ibackground,2)]);
I_fill = zeros([size(Ibackground,1)+2 size(Ibackground,2)+2]);

% Close the line if, fill color is set
if(~isempty(options.FillColor)), x=[x x(1)]; y=[y y(1)]; end

% Loop through all line coordinates
for i=1:(length(x)-1)
   % Calculate the pixels needed to construct a line of 1 pixel thickness
   % between two coordinates.
   xp=[x(i) x(i+1)];  yp=[y(i) y(i+1)]; 
   dx=abs(xp(2)-xp(1)); dy=abs(yp(2)-yp(1));
   if(dx==dy)
     if(xp(2)>xp(1)), xline=xp(1):xp(2); else xline=xp(1):-1:xp(2); end
     if(yp(2)>yp(1)), yline=yp(1):yp(2); else yline=yp(1):-1:yp(2); end
   elseif(dx>dy)
     if(xp(2)>xp(1)), xline=xp(1):xp(2); else xline=xp(1):-1:xp(2); end
     yline=linspace(yp(1),yp(2),length(xline));
   else
     if(yp(2)>yp(1)), yline=yp(1):yp(2); else yline=yp(1):-1:yp(2); end
     xline=linspace(xp(1),xp(2),length(yline));   
   end
   
   % Make closed line structure for fill if FillColor specified.
   if(~isempty(options.FillColor))
        xline_fill=xline; yline_fill=yline;
        % Limit to boundaries
        xline_fill(xline_fill<1)=1; yline_fill(yline_fill<1)=1;
        xline_fill(xline_fill>size(I_line,1))=size(I_line,1); 
        yline_fill(yline_fill>size(I_line,2))=size(I_line,2);
        % I_fill is one pixel larger than I_line to allow background fill
        xline_fill=xline_fill+1; yline_fill=yline_fill+1;
        % Insert all pixels in the fill image
        I_fill(round(xline_fill)+(round(yline_fill)-1)*size(I_fill,1))=1;
   end

   if(options.LineWidth==1)
       % Remove pixels outside image
       xline1=xline; yline1=yline;
       check=(xline1<1)|(yline1<1)|(xline1>size(I_line,1))|(yline1>size(I_line,2));
       xline1(check)=[]; yline1(check)=[];
       % Insert all pixels in the line image
       I_line(round(xline1)+round(yline1-1)*size(I_line,1))=1;
   elseif(options.LineWidth>1) % Add more pixel is line-width is larger than 1...    
       % Calculate normal on line
       ang=[yline(end)-yline(1) xline(end)-xline(1)]; ang=ang./(0.00001+sqrt(sum(ang.^2)));
       for j=-((options.LineWidth-1)/2):((options.LineWidth-1)/2);
           % Make lines close to the other lines
           xline1=xline+(ang(1)*j); yline1=yline-(ang(2)*j);
           % Remove pixels outside image
           check=(xline1<1)|(yline1<1)|(xline1>size(I_line,1))|(yline1>size(I_line,2));
           xline1(check)=[]; yline1(check)=[];
           % Insert all pixels in the line image
           I_line(ceil(xline1)+floor(yline1-1)*size(I_line,1))=1;
           I_line(floor(xline1)+floor(yline1-1)*size(I_line,1))=1;
           I_line(ceil(xline1)+ceil(yline1-1)*size(I_line,1))=1;
           I_line(floor(xline1)+ceil(yline1-1)*size(I_line,1))=1;
       end
   end
end

% Fill the line image I_fill
if(~isempty(options.FillColor))
    I_fill=bwfill(I_fill,1,1); I_fill=1-I_fill(2:end-1,2:end-1);
end

% Make marker image
if(~isempty(options.Marker))
    % Make marker pixels (center 0,0)
    switch(options.Marker)
        case '+'
            markerx=[-(options.MarkerSize/2):(options.MarkerSize/2) zeros(1,options.MarkerSize+1)];
            markery=[zeros(1,options.MarkerSize+1) -(options.MarkerSize/2):(options.MarkerSize/2)]; 
        case '*'
            markerx=[-(options.MarkerSize/2):(options.MarkerSize/2) zeros(1,options.MarkerSize+1)];
            markery=[zeros(1,options.MarkerSize+1) -(options.MarkerSize/2):(options.MarkerSize/2)]; 
            markerx=[markerx -(options.MarkerSize/2):(options.MarkerSize/2) -(options.MarkerSize/2):(options.MarkerSize/2)];
            markery=[markery -(options.MarkerSize/2):(options.MarkerSize/2) (options.MarkerSize/2):-1:-(options.MarkerSize/2)];
        case 'o'
            step=360/(2*pi*options.MarkerSize);
            markerx=options.MarkerSize/2*sind(0:step:90);
            markery=options.MarkerSize/2*cosd(0:step:90);
            markerx=[markerx -markerx markerx -markerx];
            markery=[markery markery -markery -markery];
    end
    % Add all line markers to the marker image
    for i=1:length(x);
        % Move marker to line coordinate
        xp=round(markerx)+round(x(i));  yp=round(markery)+round(y(i)); 
        % Remove outside marker pixels
        check=(xp<1)|(yp<1)|(xp>size(I_line,1))|(yp>size(I_line,2));
        xp(check)=[]; yp(check)=[];
        I_marker(xp+(yp-1)*size(I_line,1))=1;
    end
end    


% Adjust the lines and markers  with alpha value
I_line=I_line*options.Color(4);
if(~isempty(options.FillColor)), I_fill=I_fill*options.FillColor(4); end
if(~isempty(options.Marker)), I_marker=I_marker*options.MarkerColor(4); end

% Add lines, markers and fill in the right colors in the image
I=Ibackground;
if(size(Ibackground,3)==3) 
    % Color image
    for i=1:3
        if(~isempty(options.FillColor)),
            I(:,:,i)=I(:,:,i).*(1-I_fill)+options.FillColor(i)*(I_fill);
        end
        I(:,:,i)=I(:,:,i).*(1-I_line)+options.Color(i)*(I_line);
        if(~isempty(options.Marker)),
            I(:,:,i)=I(:,:,i).*(1-I_marker)+options.MarkerColor(i)*(I_marker);
        end
    end
else
    % Grey scale
    if(~isempty(options.FillColor)),
        I=I.*(1-I_fill)+mean(options.FillColor(1:3))*(I_fill);
    end
    I=I.*(1-I_line)+mean(options.Color(1:3))*(I_line);
    if(~isempty(options.Marker)),
        I=I.*(1-I_marker)+mean(options.MarkerColor(1:3))*(I_marker);
    end
end

% Set to range 0..1
I(I>1)=1; I(I<0)=0;

% Back to class background
switch (Classb)
    case 'single', I=im2single(I);
    case 'int16', I=im2int16(I);
    case 'uint8', I=im2uint8(I);
    case 'uint16', I=im2uint16(I);
end

function I=bitmaptext(lines,I,pos,options)
% The function BITMAPTEXT will insert textline(s) on the specified position
% in the image.
%
% I=bitmaptext(Text,Ibackground,Position,options)
%
% inputs,
%   Text : Cell array with text lines
%   Ibackground: the bitmap used as background when a m x n x 3 matrix
%       color plots are made, when m x n a greyscale plot. If empty []
%       autosize to fit text.
%   Position: x,y position of the text
%   options: struct with options such as color
%
% outputs,
%   Iplot: The bitmap containing the plotted text
%
% note,
%   Colors are always [r(ed) g(reen) b(lue) a(pha)], with range 0..1.
%   when Ibackground is grayscale, the mean of r,g,b is used as grey value.
%
% options,
%   options.Color: The color of the text.
%   options.FontSize: The size of the font, 1,2 or 3 (small,medium,large).
%
% example,
%
%  % The text consisting of 2 lines
%  lines={'a_A_j_J?,','ImageText version 1.1'};
%  % Background image
%  I=ones([256 256 3]);
%  % Plot text into background image
%  I=bitmaptext(lines,I,[1 1],struct('FontSize',3));
%  % Show the result
%  figure, imshow(I),
%
% Function is written by D.Kroon University of Twente (March 2009)
global character_images;

% Process inputs
defaultoptions=struct('Color',[0 0 1 1],'FontSize',1);
if(~exist('options','var')), 
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(options))), 
        warning('register_images:unknownoption','unknown options found');
    end
end

% If single line make it a cell array
if(~iscell(lines)), lines={lines}; end

if(~exist('I','var')), I=[]; end
if(exist('pos','var')),
     if(length(pos)~=2)
         error('imagtext:inputs','position must have x,y coordinates'); 
     end
else
    pos=[1 1];
end
% Round the position
pos=round(pos);

% Set the size of the font
fsize=options.FontSize;

% The character bitmap and character set;
character_set='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890!@#$%^&*()_+-=[]\;'''',./{}|:"<>?';
if(isempty(character_images)), character_images=load_font(); end

% Font parameters
Font_sizes_x=[8 10 11];
Font_sizes_y=[13 15 18];
Font_startx=[1 1 1];
Font_starty=[1 14 29];

% Get maximum sentence length
lengths=zeros(1,length(lines));
for i=1:length(lines), lengths(i)=length(lines{i}); end 
max_line_length=max(lengths);

% Make text image from the lines
lines_image=zeros([(Font_sizes_y(fsize)+4)*length(lines),max_line_length*Font_sizes_x(fsize)],'double');
for j=1:length(lines)
    line=lines{j};
    for i=1:length(line),
        [t,p]=find(character_set==line(i));
        if(~isempty(p))
            p=p(1)-1;
            character_bitmap=character_images(Font_starty(fsize):(Font_starty(fsize)+Font_sizes_y(fsize)-1),Font_startx(fsize)+(1+p*Font_sizes_x(fsize)):Font_startx(fsize)+((p+1)*Font_sizes_x(fsize)));
            posx=Font_sizes_x(fsize)*(i-1);
            posy=(Font_sizes_y(fsize)+4)*(j-1);
            lines_image((1:Font_sizes_y(fsize))+posy,(1:Font_sizes_x(fsize))+posx)=character_bitmap;
        end
    end
end

if(isempty(I)), I=zeros([size(lines_image) 3]); end

% Remove part of textimage which will be outside of the output image
if(pos(1)<1), lines_image=lines_image(2-pos(1):end,:); pos(1)=1; end
if(pos(2)<2), lines_image=lines_image(:,2-pos(2):end); pos(2)=1; end
if((pos(1)+size(lines_image,1))>size(I,1)), dif=size(I,1)-(pos(1)+size(lines_image,1)); lines_image=lines_image(1:end+dif,:); end
if((pos(2)+size(lines_image,2))>size(I,2)), dif=size(I,2)-(pos(2)+size(lines_image,2)); lines_image=lines_image(:,1:end+dif); end
% Make text image the same size as background image
I_line=zeros([size(I,1) size(I,2)]);
I_line(pos(1):(pos(1)+size(lines_image,1)-1),pos(2):(pos(2)+size(lines_image,2)-1))=lines_image;
I_line=I_line*options.Color(4);
% Insert the text image into the output image
if(~isempty(lines_image))
    if(size(I,3)==3)
        for i=1:3
            I(:,:,i)=I(:,:,i).*(1-I_line)+options.Color(i)*(I_line);
        end
    else
        I=I.*(1-I_line)+mean(options.Color(1:3))*(I_line);
    end
end

function character_images=load_font()
character_images=uint8([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 9 9 6 0 0 0 0 1 4 2 3 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 4 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 6 0 0 0 0 0 0 0 0 0 0 0 0 4 9 9 9 2 0 0 0 0 0 0 0 0 3 9 2 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 1 7 1 0 0 0 5 9 1 0 0 0 0 0 3 9 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 9 9 6 0 0 0 3 9 9 9 9 9 1 0 0 0 4 9 9 6 3 2 1 8 9 9 9 8 1 0 3 9 9 9 9 9 6 0 0 5 9 9 9 9 9 5 0 0 5 9 9 6 4 1 2 9 9 2 2 9 9 2 0 4 9 9 9 9 5 0 0 0 1 8 9 9 9 8 3 9 9 6 1 8 9 6 2 9 9 9 8 1 0 0 8 9 1 0 0 2 9 6 4 9 3 0 3 9 9 6 0 0 4 9 9 6 0 0 0 5 9 9 9 9 3 0 0 0 4 9 9 6 0 0 2 9 9 9 9 6 0 0 0 0 5 9 9 6 5 0 3 9 9 9 9 9 9 1 4 9 9 3 3 9 9 9 8 9 8 1 1 8 9 9 9 9 9 1 2 9 9 5 5 9 6 0 2 9 9 2 3 9 9 1 0 5 9 5 0 3 9 9 9 9 3 0 0 0 4 9 3 0 0 0 0 0 5 9 9 3 0 0 0 0 3 9 9 6 0 0 0 0 0 0 8 5 0 0 0 1 8 9 9 8 1 0 0 0 0 2 9 9 8 1 0 5 9 9 9 9 3 0 0 0 4 9 9 5 0 0 0 0 4 9 9 5 0 0 0 0 3 9 9 5 0 0 0 0 0 4 5 0 0 0 0 1 8 1 0 4 2 0 0 0 1 4 3 2 0 0 0 0 4 9 9 9 2 0 0 2 9 9 5 0 0 0 0 0 0 8 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 3 2 0 0 0 0 0 0 0 0 5 1 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 9 1 0 0 0 8 9 3 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 6 0 0 0 0 0 0 8 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 1 8 3 0 0 0 0 4 8 1 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 5 8 1 5 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 5 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 3 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 1 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 4 1 0 0 0 3 2 0 0 2 6 0 0 4 5 0 0 3 9 2 0 2 5 0 0 2 8 1 0 3 3 0 0 1 5 0 0 0 5 1 0 0 1 4 0 4 3 0 0 3 9 1 0 3 2 0 0 2 5 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 2 0 0 3 2 0 0 4 5 0 0 0 4 1 0 0 0 0 1 8 5 0 0 4 5 1 0 4 5 1 0 1 5 0 0 4 5 0 0 3 6 0 0 0 7 1 0 0 5 1 0 4 5 0 0 3 8 1 0 3 3 0 0 4 3 0 0 3 3 0 0 4 6 0 3 3 0 4 2 0 5 1 0 4 1 0 0 0 5 0 0 5 0 0 0 0 5 1 2 3 0 0 0 0 5 1 0 5 3 0 0 4 2 0 0 2 3 0 0 2 5 0 0 3 3 0 0 4 2 0 0 3 5 2 3 0 0 0 0 4 3 0 0 8 2 0 0 4 6 0 0 3 3 0 0 0 0 5 3 4 0 0 0 1 4 0 0 0 0 0 0 0 3 8 1 0 0 0 0 5 1 0 0 3 3 0 0 3 5 0 0 4 2 0 0 2 6 0 0 4 3 0 0 2 6 0 0 4 2 0 0 0 0 3 5 0 0 0 0 3 3 0 0 2 3 0 0 0 2 5 3 2 0 0 0 2 5 0 0 4 2 0 0 4 1 0 5 1 0 0 0 0 5 2 2 6 0 0 0 0 1 8 9 6 0 0 0 1 1 3 2 1 3 0 0 0 0 0 2 5 0 0 0 0 4 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 3 3 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 6 0 0 0 0 0 0 5 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 3 3 0 0 0 0 0 0 4 2 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 2 9 1 2 8 1 0 0 0 0 0 0 3 8 1 0 8 3 0 0 0 0 0 0 2 6 0 0 4 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 2 9 9 9 3 0 0 0 5 2 8 9 9 1 0 0 0 4 9 9 6 5 1 0 1 8 9 9 2 5 0 0 1 8 9 9 6 0 0 0 4 9 9 9 9 6 0 0 0 8 9 9 3 9 3 0 3 3 8 9 8 1 0 0 3 9 9 3 0 0 0 0 5 9 9 9 3 0 0 0 0 5 1 3 9 9 8 1 0 0 3 3 0 0 1 8 7 8 8 2 9 8 1 3 9 4 9 9 8 1 0 0 0 5 9 9 8 1 0 5 9 3 9 9 9 1 0 0 0 5 9 9 2 8 6 0 5 9 3 3 9 8 1 0 0 5 9 9 7 5 0 1 8 9 9 9 9 2 0 3 9 2 0 4 9 3 0 4 9 9 3 2 9 9 6 5 9 5 0 0 4 9 6 1 8 9 2 2 9 9 1 1 8 9 1 0 3 9 6 0 3 9 9 9 9 6 0 0 0 3 2 2 3 0 0 0 3 2 0 0 2 5 0 1 4 0 0 0 0 0 0 0 2 5 0 0 0 3 3 0 3 3 0 5 0 0 0 0 0 5 1 3 3 0 0 1 5 0 0 0 0 0 0 0 3 2 0 0 2 5 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 2 0 0 3 2 0 5 5 0 0 0 0 4 1 0 0 0 0 1 4 5 1 1 4 5 1 0 4 3 5 0 1 5 0 1 5 0 0 0 0 4 2 0 0 7 1 0 0 3 2 1 4 0 0 0 0 3 3 0 3 3 0 0 1 4 0 0 3 3 0 0 0 0 0 0 0 0 4 2 0 0 0 0 4 1 0 0 0 5 0 0 4 2 0 0 2 3 0 1 4 0 5 3 0 5 1 0 0 8 2 3 5 0 0 0 0 5 2 0 7 1 0 0 0 0 0 3 5 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 3 2 0 0 0 0 0 0 3 3 0 0 0 4 3 1 4 0 0 0 1 4 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 5 1 0 0 3 2 0 0 3 3 0 0 3 2 0 0 0 5 0 0 4 1 0 0 1 4 0 0 0 0 3 3 0 0 0 0 4 1 0 8 9 3 0 0 4 9 9 9 9 6 0 0 3 2 0 0 0 0 0 0 2 9 9 5 0 0 0 0 4 3 0 0 3 5 0 0 0 4 2 0 0 0 0 0 2 9 9 9 9 1 0 0 0 0 0 5 2 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 0 3 9 9 9 9 9 9 3 0 0 0 3 3 0 0 0 0 0 0 3 3 0 0 0 0 0 3 3 0 0 0 0 0 0 0 8 9 1 0 0 0 0 0 5 5 0 0 0 0 0 0 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 3 3 0 0 0 0 0 0 3 2 0 0 0 0 0 0 3 3 0 0 0 0 0 0 8 8 1 0 0 0 5 3 0 5 2 0 0 0 0 0 1 8 6 0 0 0 0 5 9 1 0 0 0 0 0 0 0 0 1 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 1 0 0 5 2 0 0 5 9 1 0 1 7 1 0 4 5 0 0 3 9 1 0 7 1 0 0 8 6 0 0 5 1 0 0 3 6 0 0 0 1 5 0 0 0 0 0 5 2 0 1 8 3 0 0 3 9 1 0 3 3 0 0 0 0 3 3 0 0 0 0 0 0 0 2 3 0 0 0 0 5 1 3 8 1 0 0 0 0 3 3 0 0 0 1 8 2 2 8 1 3 2 0 3 8 1 0 3 3 0 0 5 3 0 0 2 6 0 0 5 8 1 0 1 7 1 0 5 3 0 0 8 6 0 0 0 2 9 6 0 2 2 0 3 3 0 0 4 6 0 0 0 5 1 0 0 0 0 0 3 2 0 0 2 3 0 0 3 3 0 0 2 5 0 1 4 0 0 0 0 4 1 0 1 8 1 1 8 2 0 0 2 3 0 0 0 5 1 0 3 2 0 1 8 1 0 0 0 5 1 0 5 0 0 0 3 9 9 9 8 1 0 2 3 0 0 0 0 0 0 0 2 5 0 0 0 2 5 0 3 9 9 6 0 0 0 0 0 5 9 9 3 0 0 2 3 0 0 0 0 0 0 0 3 9 9 9 9 3 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 2 0 0 3 2 8 3 0 0 0 0 0 4 1 0 0 0 0 1 4 2 3 4 2 5 1 0 4 1 4 1 1 5 0 3 3 0 0 0 0 2 3 0 0 7 1 0 0 5 1 2 3 0 0 0 0 2 3 0 3 3 0 0 5 2 0 0 0 5 9 9 6 0 0 0 0 0 4 2 0 0 0 0 4 1 0 0 0 5 0 0 1 4 0 0 4 1 0 1 5 0 5 5 0 5 0 0 0 1 8 6 0 0 0 0 0 0 7 6 2 0 0 0 0 0 1 5 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 5 0 0 0 0 0 8 9 6 0 0 0 3 5 0 1 4 0 0 0 1 8 9 9 8 1 0 0 2 5 8 9 9 1 0 0 0 0 0 1 4 0 0 0 0 8 9 9 6 0 0 0 2 6 0 0 3 8 1 0 4 1 0 0 1 5 0 0 0 0 3 3 0 0 0 0 4 1 4 2 2 3 0 0 0 2 3 4 1 0 0 0 0 8 9 1 0 0 0 0 0 0 1 8 9 8 1 0 0 0 0 0 0 0 0 0 0 2 3 0 0 0 0 0 0 3 6 5 1 0 0 0 0 0 0 7 1 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 3 3 0 0 0 0 0 0 7 1 0 0 0 0 0 0 8 9 1 0 0 0 0 0 4 3 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 3 3 0 0 0 0 0 0 3 2 0 0 0 0 0 0 3 3 0 0 0 0 0 0 8 8 1 0 0 0 0 0 0 0 0 0 0 0 0 5 9 1 0 0 0 0 0 0 1 8 6 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 8 9 9 9 2 0 0 5 2 0 0 0 3 2 0 5 0 0 0 0 0 0 2 3 0 0 0 1 5 0 1 8 9 9 9 9 8 1 0 0 1 5 0 0 0 0 1 5 0 0 0 3 3 0 0 3 2 0 0 2 3 0 0 0 0 3 3 0 0 0 0 0 0 0 2 3 0 0 0 0 5 9 6 0 0 0 0 0 0 3 3 0 0 0 1 5 0 2 5 0 3 2 0 3 2 0 0 2 3 0 1 4 0 0 0 0 4 1 0 5 2 0 0 0 3 2 1 5 0 0 0 2 6 0 0 0 2 3 0 0 0 0 0 0 8 9 9 8 1 0 0 0 5 1 0 0 0 0 0 3 2 0 0 2 3 0 0 0 7 1 0 5 1 0 0 5 1 5 6 0 5 0 0 0 1 8 9 1 0 0 0 0 5 1 0 2 3 0 0 0 0 1 8 1 0 0 0 2 5 0 0 4 2 0 0 3 2 0 0 2 9 1 2 3 0 0 0 0 0 0 0 2 5 0 0 0 2 5 0 3 3 0 5 0 0 0 0 0 5 1 3 3 0 0 2 3 0 0 8 9 9 5 0 3 2 0 0 2 3 0 0 0 0 3 3 0 0 0 0 7 1 0 0 3 2 0 0 3 9 2 5 3 0 0 0 0 4 1 0 0 0 0 1 4 0 8 6 0 5 1 0 4 1 1 5 1 5 0 3 3 0 0 0 0 2 3 0 0 8 9 9 9 3 0 3 3 0 0 0 0 2 3 0 3 9 9 9 3 0 0 0 0 0 0 0 3 6 0 0 0 0 4 2 0 0 0 0 4 1 0 0 0 5 0 0 0 5 1 1 5 0 0 0 5 2 3 5 2 5 0 0 0 1 8 8 1 0 0 0 0 0 2 5 0 0 0 0 0 0 7 1 0 0 0 0 0 0 2 3 0 0 0 0 0 0 1 8 2 0 0 0 0 0 0 0 4 3 0 0 4 9 9 9 9 3 0 0 0 0 0 0 2 5 0 0 3 9 3 0 1 5 0 0 0 0 0 3 2 0 0 0 2 2 0 0 2 1 0 0 0 3 9 9 6 5 1 0 4 1 0 0 1 5 0 0 0 0 3 3 0 0 0 0 4 1 5 1 2 3 0 0 0 3 3 4 1 0 0 0 0 0 1 8 9 2 0 0 5 9 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 1 0 0 0 0 0 5 1 1 5 0 0 0 0 0 1 7 1 0 0 0 0 0 5 2 0 0 0 0 0 0 0 0 0 0 0 1 8 9 9 9 9 9 2 0 8 9 9 9 9 8 1 3 9 9 9 9 9 9 3 0 0 0 3 3 0 0 0 0 0 0 3 3 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 4 6 0 0 0 0 0 0 0 0 8 3 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 3 0 0 0 0 0 0 0 0 0 0 3 9 3 0 0 0 2 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 7 1 0 0 3 2 0 0 5 1 0 0 0 3 2 0 5 0 0 0 0 0 0 2 3 0 0 0 1 5 0 0 5 0 0 0 0 0 0 0 0 1 5 0 0 0 0 1 5 0 0 0 2 3 0 0 3 2 0 0 2 3 0 0 0 0 3 3 0 0 0 0 0 0 0 2 3 0 0 0 0 5 4 8 1 0 0 0 0 0 3 3 0 0 0 1 5 0 2 5 0 3 2 0 3 2 0 0 2 3 0 1 4 0 0 0 0 4 1 0 5 1 0 0 0 3 2 1 5 0 0 0 1 5 0 0 0 2 3 0 0 0 0 0 0 0 0 0 2 6 0 0 0 5 1 0 0 0 0 0 3 2 0 0 2 3 0 0 0 3 2 2 5 0 0 0 4 2 5 4 2 4 0 0 0 2 9 9 2 0 0 0 0 2 3 0 5 1 0 0 0 1 8 1 0 0 0 0 4 9 9 9 9 5 0 0 3 2 0 0 0 3 2 1 5 0 0 0 0 0 0 0 2 5 0 0 0 3 3 0 3 3 0 0 0 0 0 0 0 5 1 0 0 0 0 1 4 0 0 0 0 4 1 0 3 2 0 0 2 3 0 0 0 0 3 3 0 0 0 0 7 1 0 0 4 2 0 0 3 2 0 0 5 1 0 0 0 4 1 0 0 3 2 1 4 0 1 1 0 5 1 0 4 1 0 3 3 5 0 1 5 0 0 0 0 4 2 0 0 7 1 0 0 0 0 1 4 0 0 0 0 4 2 0 3 3 0 2 8 1 0 0 0 0 0 0 0 7 1 0 0 0 4 2 0 0 0 0 4 1 0 0 1 5 0 0 0 2 3 3 2 0 0 0 5 5 2 3 5 4 0 0 0 7 1 3 6 0 0 0 0 0 2 5 0 0 0 0 0 5 2 0 0 0 0 0 0 0 2 3 0 0 0 0 0 3 9 1 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 1 4 0 0 0 0 0 0 0 0 5 0 0 3 6 0 0 0 5 1 0 0 0 0 5 1 0 0 0 4 1 0 0 2 3 0 0 0 0 0 0 0 5 0 0 4 2 0 0 1 4 0 0 0 0 0 0 0 0 0 0 4 1 1 8 9 6 0 0 8 9 9 9 9 3 0 0 4 1 0 0 2 3 0 0 0 0 3 9 9 2 0 0 0 0 0 0 0 0 0 0 3 3 1 5 4 5 0 0 0 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 5 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 3 2 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 3 0 0 0 0 0 0 0 0 3 6 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 1 5 0 0 1 8 2 0 0 5 6 0 0 1 7 1 0 4 3 0 0 0 8 2 0 7 1 0 0 4 6 0 0 4 3 0 0 0 0 0 0 0 1 5 0 0 0 0 0 5 1 0 0 5 3 0 0 3 2 0 0 2 3 0 0 0 0 3 3 0 0 0 0 0 0 0 2 3 0 0 0 0 5 1 2 8 1 0 0 0 0 3 3 0 0 0 1 5 0 2 5 0 3 2 0 3 2 0 0 2 3 0 0 5 1 0 0 0 5 0 0 5 6 0 0 0 5 1 0 5 1 0 0 4 6 0 0 0 2 3 0 0 0 0 0 5 1 0 0 2 5 0 0 0 5 1 0 2 8 1 0 3 3 0 1 8 3 0 0 0 0 5 5 1 0 0 0 2 4 3 2 6 2 0 0 3 8 1 0 8 3 0 0 0 0 5 2 3 0 0 0 1 8 1 0 0 5 0 0 5 0 0 0 0 5 1 0 3 2 0 0 0 7 1 0 4 3 0 0 0 8 2 0 2 5 0 0 0 7 1 0 3 3 0 0 0 5 1 0 0 5 1 0 0 0 0 0 5 1 0 0 0 4 1 0 3 2 0 0 2 3 0 0 0 0 3 3 0 0 0 0 5 3 0 1 8 1 0 0 3 2 0 0 2 5 0 0 0 4 1 0 0 3 2 1 4 0 0 0 0 5 1 0 4 1 0 0 8 6 0 0 4 5 0 0 3 6 0 0 0 7 1 0 0 0 0 0 4 3 0 0 3 6 0 0 3 3 0 0 3 5 0 0 7 1 0 0 1 4 0 0 0 0 4 2 0 0 0 0 3 5 0 0 3 5 0 0 0 0 8 8 1 0 0 0 5 5 1 1 8 3 0 0 5 2 0 0 4 5 0 0 0 0 2 5 0 0 0 0 3 3 0 0 0 5 0 0 0 0 2 3 0 0 0 0 4 6 0 0 0 0 0 0 5 2 0 0 4 5 0 0 0 0 0 1 4 0 0 0 5 3 0 0 3 5 0 0 0 5 2 0 0 5 0 0 0 0 1 4 0 0 0 0 3 2 0 0 2 3 0 0 0 0 0 0 4 2 0 0 1 3 0 0 4 2 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 3 2 5 1 0 0 0 4 9 9 9 8 1 0 0 0 0 5 0 0 5 0 0 0 0 0 0 0 0 0 0 3 5 0 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 3 3 0 0 0 0 2 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 9 5 0 0 0 0 0 0 8 8 1 0 0 0 0 7 1 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 3 2 0 0 0 0 0 0 3 3 0 0 0 0 0 0 8 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 1 0 0 0 0 1 8 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 3 9 9 9 4 9 2 4 9 4 9 9 9 1 0 0 0 5 9 9 9 2 0 0 2 9 9 9 6 8 6 0 0 5 9 9 9 8 1 0 5 9 9 9 9 5 0 0 1 8 9 9 6 3 0 3 9 9 2 2 9 9 3 0 8 9 9 9 9 8 1 0 0 0 0 2 3 0 0 0 5 9 1 4 9 9 3 0 8 9 9 9 9 8 1 8 9 5 2 9 3 3 9 6 9 9 2 2 9 9 3 0 1 8 9 9 9 2 0 0 5 4 9 9 9 3 0 0 1 8 9 9 6 5 0 0 8 9 9 9 9 2 0 0 5 9 9 9 8 1 0 0 0 2 9 9 8 1 0 0 1 8 9 9 3 9 3 0 0 0 3 5 0 0 0 0 1 8 1 0 7 1 0 3 9 9 2 2 9 9 3 0 0 0 2 8 1 0 0 0 4 9 9 9 9 6 1 8 9 9 1 0 8 9 9 4 9 9 9 9 9 3 0 0 0 5 9 9 9 2 0 1 8 9 9 9 9 2 0 3 9 9 9 9 9 8 1 0 5 9 9 9 1 0 0 0 1 8 9 9 9 5 0 3 9 9 2 2 9 9 3 0 4 9 9 9 9 5 0 0 0 5 9 9 2 0 0 3 9 9 6 0 0 8 6 2 9 9 9 9 9 9 3 8 9 9 1 1 8 9 8 5 9 9 3 0 2 6 0 0 0 4 9 9 6 0 0 0 5 9 9 9 1 0 0 0 0 5 9 9 6 0 0 2 9 9 6 0 0 5 6 0 8 9 9 9 9 1 0 0 3 9 9 9 8 1 0 0 0 4 9 9 6 0 0 0 0 0 3 3 0 0 0 0 4 5 0 0 5 3 0 5 9 8 1 1 8 9 5 0 1 8 9 9 9 2 0 0 4 9 9 9 9 6 0 0 4 9 9 9 9 6 0 1 8 9 9 9 9 2 0 0 0 8 9 9 6 0 0 0 0 0 4 9 9 3 0 0 0 5 9 9 6 0 0 0 0 0 8 9 9 2 0 0 0 0 3 2 0 0 0 0 0 8 9 9 8 1 0 0 2 9 9 9 5 0 0 0 0 5 9 9 5 0 0 0 0 0 5 6 0 0 0 0 1 7 1 0 5 2 0 0 0 4 1 5 1 0 0 0 0 0 3 2 0 0 0 0 0 0 3 9 9 2 0 0 0 0 0 0 0 0 0 0 0 5 9 9 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 5 2 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 7 1 0 0 0 4 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 1 0 0 0 0 0 0 8 8 1 0 0 0 3 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 4 2 0 0 0 0 0 0 3 3 0 0 0 0 0 0 8 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 6 0 0 8 9 1 0 0 0 0 0 0 0 8 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 1 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 9 3 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 9 5 0 0 0 0 4 2 5 0 0 0 0 0 0 3 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 5 0 0 0 0 4 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 3 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 3 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 1 8 3 0 0 0 0 4 8 1 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 9 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 9 6 0 0 0 0 0 0 0 0 4 9 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 9 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 1 0 8 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 1 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 9 1 0 0 0 8 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 9 9 9 9 9 9 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 3 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 5 2 5 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 4 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 9 9 9 1 0 0 0 0 0 0 0 0 0 0 3 9 6 0 0 0 0 0 0 0 0 0 0 0 5 3 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 3 9 6 0 0 0 0 0 0 0 0 8 9 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 3 0 0 0 0 0 0 1 8 9 9 3 0 0 0 0 0 0 8 9 9 8 1 0 0 0 0 0 0 0 5 8 1 0 0 0 0 4 9 9 9 9 5 0 0 0 0 0 0 1 8 9 9 5 0 0 2 9 9 9 9 9 8 1 0 0 0 0 8 9 9 6 0 0 0 0 0 0 5 9 9 6 0 0 0 0 0 0 3 9 9 6 0 0 0 0 0 0 0 4 6 0 0 0 0 0 0 0 4 9 9 8 1 0 0 0 0 0 3 3 2 5 0 0 0 0 0 0 5 9 9 9 6 0 0 0 0 4 9 9 3 0 0 0 0 0 0 0 0 8 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 2 6 0 0 0 0 4 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 9 5 0 0 0 0 3 9 9 3 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 9 1 0 0 0 0 0 0 1 8 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 6 0 2 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 9 9 8 1 0 0 0 2 9 9 9 9 9 9 2 0 0 0 0 0 8 9 9 8 1 7 1 0 8 9 9 9 9 9 3 0 0 2 9 9 9 9 9 9 9 3 0 0 4 9 9 9 9 9 9 9 2 0 0 0 8 9 9 9 2 5 0 0 8 9 9 1 1 8 9 8 1 0 2 9 9 9 9 9 9 2 0 0 0 0 2 9 9 9 9 9 9 3 9 9 9 5 3 9 9 8 1 2 9 9 9 9 9 2 0 0 0 5 9 8 1 0 0 0 5 9 6 3 9 9 2 0 1 8 9 9 5 0 0 0 5 9 9 9 2 0 0 0 3 9 9 9 9 9 6 0 0 0 0 0 5 9 9 9 2 0 0 1 8 9 9 9 9 9 1 0 0 0 0 0 8 9 9 8 3 3 0 1 8 9 9 9 9 9 9 6 0 4 9 9 9 2 2 9 9 9 6 8 9 9 8 1 0 8 9 9 6 5 9 9 6 0 1 8 9 9 3 4 9 9 3 0 0 8 9 9 1 2 9 9 6 0 0 4 9 9 3 0 0 8 9 9 9 9 9 1 0 0 2 9 9 3 3 0 0 0 0 0 1 8 2 0 0 8 3 0 0 0 2 9 2 0 0 3 8 1 0 0 0 0 0 3 7 7 1 0 0 0 0 4 2 0 0 0 0 0 0 0 0 0 3 9 1 0 0 0 0 0 2 5 0 0 0 0 7 1 0 0 0 5 2 0 0 3 5 0 0 0 0 4 5 0 0 3 6 0 0 0 0 3 6 0 0 3 5 0 0 0 0 0 0 4 6 0 0 0 0 0 0 4 6 0 0 3 6 0 0 0 0 0 3 3 2 3 0 0 0 0 0 5 3 0 0 2 6 0 0 0 2 5 0 0 5 1 0 0 0 0 0 0 5 5 4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 8 1 0 0 0 0 0 0 1 8 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 8 6 0 0 0 0 0 0 5 8 1 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 1 0 5 6 0 0 0 0 0 0 0 0 0 5 6 0 0 5 5 0 0 0 0 0 0 0 0 0 0 5 9 9 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 5 4 2 0 0 0 0 0 7 1 0 0 1 8 1 0 0 2 9 2 0 0 2 9 8 1 0 0 5 1 0 0 0 5 3 0 0 0 7 1 0 0 0 3 3 0 0 0 2 5 0 0 0 0 4 2 0 1 8 2 0 0 1 8 6 0 0 1 5 0 0 0 0 5 1 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 7 1 0 0 5 5 0 0 0 0 0 7 1 0 0 0 0 0 0 7 6 3 0 0 2 7 5 0 0 2 9 6 0 0 0 4 3 0 0 1 8 3 0 0 1 8 3 0 0 0 2 5 0 0 0 3 6 0 0 0 8 3 0 0 0 8 5 0 0 0 7 1 0 0 2 9 1 0 0 0 8 2 0 0 2 9 3 0 1 5 0 0 4 2 0 1 5 0 0 2 5 0 0 0 0 3 3 0 0 3 3 0 0 0 0 3 5 0 1 5 0 0 0 0 0 2 5 0 0 3 6 0 0 0 2 8 1 0 0 0 7 1 0 0 0 5 1 0 0 0 7 1 0 0 1 7 1 0 0 0 0 0 2 3 0 0 0 0 0 3 3 0 0 0 1 5 0 0 0 0 0 0 0 0 0 5 1 0 0 0 0 1 8 2 7 1 0 0 0 0 4 2 0 0 0 0 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 2 5 0 0 0 1 5 0 0 0 0 7 1 0 0 0 7 1 0 0 0 3 2 0 0 0 5 0 0 0 0 5 1 0 0 0 0 0 4 5 0 0 0 0 0 0 7 1 0 0 0 7 1 0 0 0 0 4 2 3 3 0 0 0 0 0 5 0 0 0 0 0 0 0 0 2 5 0 0 5 1 0 0 0 0 0 4 6 0 0 5 5 0 0 0 0 0 2 9 9 9 1 0 0 0 2 9 1 3 3 1 8 2 0 0 0 0 0 0 3 6 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 6 0 0 0 0 0 0 0 0 8 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 4 2 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 3 0 2 9 2 0 0 0 0 0 0 0 3 9 3 0 0 0 0 4 9 2 0 0 0 0 0 0 0 8 3 0 0 3 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 8 9 9 9 5 0 0 0 0 3 3 2 9 9 9 3 0 0 0 0 0 8 9 9 9 3 5 0 0 0 3 9 9 9 3 3 3 0 0 0 4 9 9 9 8 1 0 0 0 1 8 9 9 9 9 9 3 0 0 0 3 9 9 9 2 8 9 2 0 1 5 3 9 9 9 2 0 0 0 0 8 9 9 3 0 0 0 0 0 3 9 9 9 9 8 1 0 0 0 0 1 5 0 3 9 9 9 1 0 0 0 0 3 3 0 0 0 2 9 6 5 9 5 2 9 9 2 0 2 9 8 1 8 9 9 2 0 0 0 0 1 8 9 9 9 1 0 0 4 9 3 3 9 9 9 3 0 0 0 0 2 9 9 9 3 3 9 5 0 2 9 9 2 2 9 9 6 0 0 0 1 8 9 9 8 5 2 0 0 8 9 9 9 9 9 6 0 0 3 9 6 0 0 3 9 9 1 0 3 9 9 9 2 1 8 9 9 5 5 9 9 2 0 0 2 9 9 5 0 8 9 9 1 1 8 9 8 1 0 8 9 9 1 0 1 8 9 5 0 0 8 9 9 9 9 9 2 0 0 0 0 4 2 2 6 0 0 0 0 0 7 1 0 0 0 4 2 0 0 5 1 0 0 0 0 1 7 1 0 0 5 1 0 0 0 0 7 1 0 0 7 1 0 0 0 3 3 0 0 0 2 5 0 0 0 0 4 2 0 5 1 0 0 0 0 0 0 0 0 1 5 0 0 0 0 5 1 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 7 1 0 8 3 0 0 0 0 0 0 7 1 0 0 0 0 0 0 7 2 8 1 0 4 3 5 0 0 2 5 5 3 0 0 4 3 0 0 4 3 0 0 0 0 1 8 1 0 0 2 5 0 0 0 0 5 1 0 4 3 0 0 0 0 0 5 2 0 0 7 1 0 0 0 4 2 0 0 2 6 0 0 0 0 3 3 0 1 5 0 0 4 2 0 1 5 0 0 2 5 0 0 0 0 3 3 0 0 1 5 0 0 0 0 5 1 0 0 7 1 0 8 3 0 2 5 0 0 0 3 6 0 2 8 1 0 0 0 0 2 6 0 0 4 3 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 2 5 0 0 0 0 0 0 0 0 3 6 0 0 0 0 0 5 3 0 7 1 0 0 0 0 4 2 0 0 0 0 0 0 0 0 5 2 0 0 0 0 0 0 0 0 0 0 0 0 4 2 0 0 0 0 7 1 0 0 1 5 0 0 0 0 7 1 0 0 0 3 3 0 0 1 4 0 0 0 0 4 2 0 0 0 0 0 4 5 0 0 0 0 0 2 6 0 0 5 9 8 1 0 0 2 9 9 9 9 9 9 2 0 0 0 7 1 0 0 0 0 0 0 0 0 4 9 9 3 0 0 0 0 0 2 8 1 0 0 0 8 2 0 0 0 0 7 1 0 0 0 0 0 0 0 1 8 9 9 9 1 0 0 0 0 0 0 0 5 3 0 0 0 0 0 0 2 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 9 9 9 9 9 9 9 9 2 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 8 9 2 0 0 0 0 0 0 0 5 5 0 0 0 0 0 0 0 0 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 2 9 9 1 0 0 0 0 3 6 0 0 5 5 0 0 0 0 0 0 0 8 6 0 0 0 0 0 0 0 0 8 8 1 0 0 0 0 0 7 1 0 0 0 4 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 4 6 0 0 0 3 6 8 1 0 0 5 3 0 0 1 8 2 0 0 0 8 6 0 0 3 6 0 0 0 5 9 3 0 0 3 5 0 0 0 2 8 1 0 0 0 0 1 5 0 0 0 0 0 0 3 6 0 0 1 8 9 1 0 0 1 8 6 0 0 2 9 1 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 1 5 0 0 1 1 0 0 0 0 0 0 3 3 0 0 0 0 1 8 3 0 5 8 1 2 8 1 0 0 8 9 1 0 1 8 1 0 0 2 9 1 0 0 1 8 2 0 0 3 6 6 0 0 0 5 3 0 0 3 8 1 0 0 5 6 3 0 0 0 0 4 9 9 1 0 1 1 0 0 8 2 0 0 2 9 2 0 0 0 1 5 0 0 0 0 0 0 0 1 5 0 0 0 0 5 1 0 0 0 7 1 0 0 0 5 2 0 0 5 1 0 0 0 0 1 5 0 0 0 5 5 0 0 4 6 0 0 0 1 8 1 0 0 0 2 6 0 0 0 7 1 0 0 3 6 0 0 0 0 1 5 0 0 5 1 0 0 0 0 7 1 0 0 1 8 1 0 1 5 0 0 0 0 0 0 0 0 0 0 5 1 0 0 0 0 4 3 0 0 7 1 0 7 1 0 0 0 0 0 2 5 0 3 5 0 0 0 1 7 1 0 0 0 0 0 0 0 0 1 5 0 0 0 0 5 1 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 7 1 8 2 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 7 1 5 3 1 5 1 5 0 0 2 5 1 7 1 0 4 3 0 0 7 1 0 0 0 0 0 4 3 0 0 2 5 0 0 0 0 5 1 0 7 1 0 0 0 0 0 3 3 0 0 7 1 0 0 3 8 1 0 0 0 8 2 0 0 0 0 0 0 1 5 0 0 4 2 0 1 5 0 0 2 5 0 0 0 0 3 3 0 0 0 4 2 0 0 2 6 0 0 0 5 1 1 6 5 0 3 3 0 0 0 0 4 7 8 1 0 0 0 0 0 0 3 5 3 5 0 0 0 0 0 0 0 0 5 2 0 0 0 0 0 0 0 2 5 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 1 8 9 6 0 0 0 0 0 3 6 0 0 7 1 0 0 0 0 4 9 9 9 9 2 0 0 0 0 7 1 5 9 9 3 0 0 0 0 0 0 0 0 7 1 0 0 0 0 1 8 9 9 9 1 0 0 0 0 4 5 0 0 3 9 5 0 0 2 6 0 0 0 0 4 3 0 0 0 0 0 3 5 0 0 0 0 0 2 5 0 4 5 0 7 1 0 0 0 0 4 2 4 2 0 0 0 0 0 2 9 9 3 0 0 0 0 0 0 0 0 1 8 9 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 1 0 0 0 0 0 0 0 0 4 5 5 3 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 8 9 2 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 2 9 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 2 0 0 0 0 0 0 0 0 0 0 2 9 3 0 0 0 0 0 0 0 0 0 5 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 1 5 0 0 0 3 6 0 0 0 0 0 7 1 0 4 3 0 0 0 0 0 0 0 0 7 1 0 0 0 0 4 3 0 0 7 1 0 0 0 0 2 5 0 0 0 0 1 5 0 0 0 0 0 0 7 1 0 0 0 1 8 1 0 0 1 5 0 0 0 0 5 1 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 1 5 0 8 6 0 0 0 0 0 0 0 3 3 0 0 0 0 1 5 0 0 4 3 0 1 7 1 0 0 8 2 0 0 0 5 1 0 0 5 1 0 0 0 0 1 5 0 0 3 6 0 0 0 0 0 7 1 0 7 1 0 0 0 0 5 3 0 0 0 0 4 5 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 1 5 0 0 0 0 5 1 0 0 0 3 3 0 0 1 5 0 0 0 3 2 0 5 6 0 2 3 0 0 0 0 4 6 4 5 0 0 0 0 0 4 6 0 0 0 5 2 0 0 0 0 0 0 4 6 0 0 0 0 0 3 3 0 0 3 5 0 0 0 0 8 9 9 9 9 3 0 0 1 5 0 0 0 0 0 0 0 0 0 0 5 1 0 0 0 0 3 3 0 0 8 9 9 8 1 0 0 0 0 0 2 9 9 9 5 0 0 0 1 5 0 0 0 0 0 0 0 0 0 1 8 9 9 9 9 9 1 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 8 9 2 8 3 0 0 0 0 0 0 7 1 0 0 0 0 0 0 7 1 2 8 5 2 1 5 0 0 2 5 0 4 5 0 4 3 0 1 5 0 0 0 0 0 0 3 3 0 0 2 5 0 0 0 2 6 0 1 5 0 0 0 0 0 0 3 3 0 0 8 9 9 9 6 0 0 0 0 0 0 8 9 9 9 1 0 0 0 0 0 0 4 2 0 0 0 0 0 2 5 0 0 0 0 3 3 0 0 0 2 5 0 0 4 2 0 0 0 5 1 3 3 5 1 4 2 0 0 0 0 0 8 3 0 0 0 0 0 0 0 0 4 8 1 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 2 5 0 0 0 0 0 0 0 0 1 8 2 0 0 0 0 0 0 0 0 0 4 6 0 0 0 1 8 1 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 7 5 3 0 0 5 3 0 0 0 0 0 0 2 5 0 0 0 0 0 4 1 0 0 1 3 0 0 0 0 0 5 9 9 6 2 5 0 0 2 6 0 0 0 0 4 3 0 0 0 0 0 3 5 0 0 0 0 0 2 5 0 7 1 0 7 1 0 0 0 0 5 1 4 2 0 0 0 0 0 0 0 0 5 9 5 0 0 0 3 9 9 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 8 1 0 0 0 0 0 0 1 7 1 1 8 1 0 0 0 0 0 0 2 8 1 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 9 9 9 9 9 9 8 1 0 4 9 9 9 9 9 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 1 0 0 0 0 0 2 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 8 9 9 9 9 6 0 0 0 3 3 0 0 0 0 0 5 1 0 4 2 0 0 0 0 0 0 0 1 5 0 0 0 0 0 3 3 0 1 8 9 9 9 9 9 9 6 0 0 0 0 1 5 0 0 0 0 0 1 5 0 0 0 0 0 7 1 0 0 1 5 0 0 0 0 5 1 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 1 8 9 5 0 0 0 0 0 0 0 0 3 3 0 0 0 0 1 5 0 0 4 3 0 1 7 1 0 0 7 1 0 0 0 5 1 0 0 7 1 0 0 0 0 0 7 1 0 3 3 0 0 0 0 0 5 1 1 5 0 0 0 0 0 3 3 0 0 0 0 4 2 0 0 0 0 0 0 0 2 9 9 9 9 3 0 0 0 0 1 5 0 0 0 0 0 0 0 1 5 0 0 0 0 5 1 0 0 0 0 7 1 0 4 2 0 0 0 2 5 1 5 5 1 4 2 0 0 0 0 0 5 6 0 0 0 0 0 0 0 8 2 0 2 5 0 0 0 0 0 0 4 5 0 0 0 0 0 0 8 9 9 9 9 8 1 0 0 0 7 1 0 0 0 5 6 0 1 7 1 0 0 0 0 0 0 0 0 0 5 1 0 0 0 0 4 2 0 0 7 1 0 7 1 0 0 0 0 0 2 5 0 3 5 0 0 0 1 7 1 0 0 8 9 9 9 5 0 1 5 0 0 0 0 5 1 0 0 0 0 0 3 3 0 0 0 0 0 4 2 0 0 0 1 5 0 0 0 0 7 1 0 0 5 2 0 0 0 0 0 7 1 0 0 0 5 1 0 7 1 0 5 8 1 1 5 0 0 2 5 0 0 8 2 4 3 0 0 7 1 0 0 0 0 0 4 3 0 0 2 9 9 9 9 8 1 0 0 7 1 0 0 0 0 0 3 3 0 0 7 1 0 2 9 1 0 0 0 0 0 0 0 0 1 8 3 0 0 0 0 0 4 2 0 0 0 0 0 2 5 0 0 0 0 3 3 0 0 0 0 5 1 0 7 1 0 0 0 4 2 5 1 3 3 4 2 0 0 0 0 5 2 5 3 0 0 0 0 0 0 0 2 5 0 0 0 0 0 0 0 3 5 0 0 3 3 0 0 0 0 0 2 5 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 0 0 0 0 0 4 2 0 0 3 9 9 9 9 9 8 1 0 0 0 0 0 0 0 0 3 3 0 0 0 8 5 0 0 0 0 5 0 0 0 0 0 0 4 2 0 0 0 0 1 5 0 0 0 0 7 1 0 0 0 0 0 0 0 0 3 3 0 0 1 5 0 0 0 0 4 2 0 0 0 0 0 3 3 0 0 0 0 0 2 5 0 5 2 0 7 1 0 0 4 9 9 9 9 9 8 1 0 0 0 0 0 0 0 0 5 1 0 0 0 0 0 3 9 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 8 5 6 0 8 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 8 1 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 9 9 9 9 9 9 9 9 2 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 8 6 0 0 0 0 0 0 0 0 0 0 5 8 1 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 2 0 0 0 0 0 0 0 0 0 0 2 9 3 0 0 0 0 0 0 3 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 8 2 0 0 0 1 5 0 0 0 3 6 0 0 0 0 0 7 1 0 4 3 0 0 0 0 0 0 0 0 7 1 0 0 0 0 4 3 0 0 7 1 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 7 1 0 0 0 1 8 1 0 0 1 5 0 0 0 0 5 1 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 1 7 1 4 3 0 0 0 0 0 0 0 3 3 0 0 0 0 1 5 0 0 4 3 0 1 7 1 0 0 7 1 0 0 0 5 1 0 0 5 1 0 0 0 0 1 5 0 0 3 5 0 0 0 0 0 7 1 0 7 1 0 0 0 0 4 3 0 0 0 0 4 2 0 0 0 0 0 0 0 0 0 0 0 0 8 2 0 0 0 1 5 0 0 0 0 0 0 0 1 5 0 0 0 0 5 1 0 0 0 0 4 3 1 5 0 0 0 0 0 5 3 3 3 3 5 1 0 0 0 0 8 3 3 8 1 0 0 0 0 0 3 6 0 5 1 0 0 0 0 0 4 5 0 0 0 0 0 0 2 5 0 0 0 0 4 3 0 0 0 7 1 0 0 0 0 5 1 0 5 2 0 0 0 0 0 0 0 0 0 5 1 0 0 0 0 5 1 0 0 7 1 0 0 0 0 0 0 0 0 2 5 0 0 0 0 0 0 0 5 1 0 0 0 0 1 5 0 0 1 5 0 0 0 0 5 1 0 0 0 0 0 3 3 0 0 0 0 0 4 2 0 0 0 2 6 0 0 0 0 7 1 0 0 2 6 0 0 0 0 0 7 1 0 0 0 5 1 0 7 1 0 0 0 0 1 5 0 0 2 5 0 0 2 6 4 3 0 0 4 3 0 0 0 0 1 8 1 0 0 2 5 0 0 0 0 0 0 0 4 2 0 0 0 0 0 7 1 0 0 7 1 0 0 2 6 0 0 0 4 2 0 0 0 0 2 5 0 0 0 0 0 4 2 0 0 0 0 0 2 5 0 0 0 0 4 3 0 0 0 0 3 3 3 5 0 0 0 0 3 5 5 0 2 5 5 1 0 0 0 5 3 0 0 8 2 0 0 0 0 0 0 2 5 0 0 0 0 0 0 2 6 0 0 0 3 3 0 0 0 0 0 2 5 0 0 0 0 0 0 2 9 1 0 0 0 0 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 3 3 0 0 0 4 3 0 0 0 0 5 0 0 0 0 0 0 7 1 0 0 0 0 1 4 0 0 0 0 5 1 0 0 0 0 0 0 0 0 7 1 0 0 0 7 1 0 0 0 5 1 0 0 0 0 0 0 0 0 0 0 0 0 2 5 0 1 8 9 9 3 0 0 0 0 7 1 5 1 0 0 0 0 2 6 0 0 0 1 7 1 0 0 0 0 1 5 0 0 4 2 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 5 6 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 3 9 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 8 1 0 0 0 0 0 0 0 8 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 7 1 0 0 0 8 6 0 0 0 3 6 6 0 0 0 5 3 0 0 1 8 1 0 0 0 3 9 1 0 3 5 0 0 0 3 9 3 0 0 2 6 0 0 0 0 3 6 0 0 0 0 1 5 0 0 0 0 0 0 3 5 0 0 0 5 9 1 0 0 1 5 0 0 0 0 5 1 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 1 5 0 0 5 3 0 0 0 0 0 0 3 3 0 0 0 0 1 5 0 0 4 3 0 1 7 1 0 0 7 1 0 0 0 5 1 0 0 2 6 0 0 0 0 5 2 0 0 3 9 5 0 0 0 3 3 0 0 3 3 0 0 0 3 9 3 0 0 0 0 4 2 0 0 0 0 0 0 3 5 0 0 0 0 5 2 0 0 0 1 7 1 0 0 8 5 0 0 1 8 1 0 0 8 9 1 0 0 0 0 1 7 5 2 0 0 0 0 0 5 6 1 1 6 5 0 0 0 1 8 2 0 0 2 9 1 0 0 0 0 0 8 4 5 0 0 0 0 0 5 5 0 0 0 3 3 0 0 5 1 0 0 0 0 1 5 0 0 0 7 1 0 0 0 2 8 1 0 0 7 1 0 0 0 3 8 1 0 0 5 1 0 0 0 3 3 0 0 0 7 1 0 0 0 1 4 0 0 0 2 5 0 0 0 0 0 0 0 2 6 0 0 0 0 1 5 0 0 1 5 0 0 0 0 5 1 0 0 0 0 0 3 3 0 0 0 0 0 3 8 1 0 0 8 3 0 0 0 0 7 1 0 0 0 5 2 0 0 0 0 7 1 0 0 0 5 1 0 7 1 0 0 0 0 1 5 0 0 2 5 0 0 0 5 9 3 0 0 0 8 5 0 0 1 8 3 0 0 0 2 5 0 0 0 0 0 0 0 1 8 2 0 0 1 8 3 0 0 0 7 1 0 0 0 4 3 0 0 4 9 1 0 0 0 4 2 0 0 0 0 0 4 2 0 0 0 0 0 0 8 2 0 0 1 7 1 0 0 0 0 0 7 6 1 0 0 0 0 2 6 3 0 0 8 8 1 0 0 4 5 0 0 0 0 8 2 0 0 0 0 0 2 5 0 0 0 0 0 1 7 1 0 0 0 3 3 0 0 0 0 0 2 5 0 0 0 0 0 3 8 1 0 0 1 5 0 0 0 4 6 0 0 0 3 8 1 0 0 0 0 0 0 0 7 1 0 0 0 3 8 1 0 0 2 9 1 0 0 0 1 8 1 0 0 5 2 0 0 0 0 0 2 6 0 0 0 0 0 0 7 1 0 0 2 6 0 0 0 0 0 0 0 0 8 2 0 0 0 0 3 3 0 0 2 5 0 0 0 0 0 0 8 8 1 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 7 1 5 1 0 0 0 0 2 9 9 9 9 9 1 0 0 0 0 0 1 5 0 0 4 2 0 0 0 0 0 0 0 0 0 0 0 0 0 5 3 0 1 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 3 0 0 0 0 0 0 2 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 5 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 1 0 0 0 0 0 0 0 2 9 9 1 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 2 9 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 3 0 0 0 0 3 9 2 0 0 0 0 0 0 0 0 1 8 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 2 9 9 9 9 3 8 8 1 4 9 3 3 9 9 9 3 0 0 0 0 1 8 9 9 9 6 0 0 0 0 4 9 9 9 6 3 9 5 0 0 3 9 9 9 9 6 0 0 0 2 9 9 9 9 9 9 2 0 0 0 4 9 9 9 4 7 1 0 3 9 9 8 1 0 8 9 9 3 0 4 9 9 9 9 9 9 5 0 0 0 0 0 0 1 7 1 0 0 0 3 9 6 0 1 8 9 9 3 0 4 9 9 9 9 9 9 5 2 9 9 9 1 4 9 5 1 8 9 3 9 9 9 1 0 8 9 9 3 0 0 3 9 9 9 9 3 0 0 0 3 3 4 9 9 9 6 0 0 0 0 5 9 9 9 6 3 3 0 0 4 9 9 9 9 9 6 0 0 0 3 9 9 9 9 9 3 0 0 0 0 0 3 9 9 9 3 0 0 0 0 3 9 9 9 2 5 9 3 0 0 0 0 4 8 1 0 0 0 0 0 3 6 0 0 5 3 0 0 2 9 9 9 1 1 8 9 9 2 0 0 0 0 3 9 1 0 0 0 0 2 9 9 9 9 9 9 3 0 8 9 9 6 0 0 4 9 9 9 4 9 9 9 9 9 9 9 1 0 0 0 1 8 9 9 9 6 0 0 0 8 9 9 9 9 9 6 0 0 2 9 9 9 9 9 9 9 5 0 0 4 9 9 9 9 2 0 0 0 0 0 3 9 9 9 9 9 1 0 2 9 9 9 1 1 8 9 9 2 0 2 9 9 9 9 9 9 2 0 0 0 2 9 9 9 3 0 0 0 2 9 9 9 5 0 0 3 9 5 2 9 9 9 9 9 9 9 9 2 8 9 9 5 0 0 5 9 9 9 4 9 9 9 2 0 1 8 3 0 0 0 0 5 9 9 9 2 0 0 0 3 9 9 9 9 3 0 0 0 0 0 1 8 9 9 9 2 0 0 1 8 9 9 5 0 0 1 8 5 0 4 3 8 9 9 9 5 0 0 0 0 8 9 9 9 9 5 0 0 0 0 0 8 9 9 9 1 0 0 0 0 0 0 4 6 0 0 0 0 0 2 9 2 0 0 4 6 0 0 4 9 9 6 0 0 8 9 9 3 0 0 4 9 9 9 9 6 0 0 0 2 9 9 9 9 9 9 3 0 0 1 8 9 9 9 9 9 3 0 0 5 9 9 9 9 9 6 0 0 0 0 3 9 9 9 8 1 0 0 0 0 0 0 5 9 9 8 1 0 0 0 2 9 9 9 8 1 0 0 0 0 0 1 8 9 9 3 0 0 0 0 0 0 4 3 0 0 0 0 0 0 1 8 9 9 8 1 0 0 0 0 5 9 9 9 2 0 0 0 0 0 0 5 9 9 8 1 0 0 0 0 0 0 8 8 1 0 0 0 0 0 4 6 0 0 3 5 0 0 0 0 0 5 0 7 1 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 3 9 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 9 9 6 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 6 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 1 8 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 5 0 0 0 0 0 0 0 0 2 9 9 1 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 2 9 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 0 0 5 6 0 0 0 0 0 0 0 0 0 0 1 8 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 9 9 8 1 4 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 9 8 1 0 0 0 0 1 5 0 7 1 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 4 2 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 1 0 3 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 6 0 0 0 0 4 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 9 5 0 0 0 0 3 9 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 6 0 0 0 0 0 0 8 8 1 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 9 9 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 9 9 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 9 9 5 0 0 0 0 0 0 0 0 0 0 4 9 9 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 9 9 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 9 9 9 9 9 9 9 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 1 0 0 0 0 0 0 0 0 0 4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 9 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 4 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 9 9 9 9 1 0 0 0 0 0 0 0 0 0 0 0 2 9 9 1 0 0 0 0 0 0 0 0 0 0 0 3 9 1 0 0 0 0 0 0 0 0 0 4 6 0 0 0 0 0 1 8 9 1 0 0 0 0 0 0 0 0 5 9 9 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 2 0 0 0 0 0 0 0 3 9 9 9 3 0 0 0 0 0 0 3 9 9 9 5 0 0 0 0 0 0 0 0 1 8 6 0 0 0 0 0 3 9 9 9 9 9 5 0 0 0 0 0 0 0 2 9 9 9 5 0 0 3 9 9 9 9 9 9 8 1 0 0 0 0 3 9 9 9 2 0 0 0 0 0 0 2 9 9 9 2 0 0 0 0 0 0 1 8 9 9 3 0 0 0 0 0 0 0 2 9 2 0 0 0 0 0 0 0 8 3 0 0 8 3 0 0 0 0 0 0 4 2 1 7 1 0 0 0 0 0 3 9 9 9 9 5 0 0 0 0 1 8 9 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 0 0 4 5 0 0 0 0 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 9 9 3 0 0 0 0 3 9 9 9 1 0 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 9 5 0 0 0 0 0 0 0 0 8 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 2 9 5 0 0 0 0 0 0 4 9 1 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 3 0 4 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 9 9 9 6 0 0 0 0 3 9 9 9 9 9 9 8 1 0 0 0 0 0 3 9 9 9 5 2 6 0 1 8 9 9 9 9 9 8 1 0 0 2 9 9 9 9 9 9 9 9 2 0 0 3 9 9 9 9 9 9 9 9 5 0 0 0 3 9 9 9 6 2 6 0 0 8 9 9 5 0 4 9 9 8 1 0 1 8 9 9 9 9 9 9 1 0 0 0 0 0 8 9 9 9 9 9 8 2 9 9 9 9 1 0 8 9 9 3 1 8 9 9 9 9 5 0 0 0 0 5 9 9 1 0 0 0 2 9 9 5 5 9 9 3 0 0 4 9 9 9 6 0 0 0 3 9 9 9 5 0 0 0 0 3 9 9 9 9 9 9 3 0 0 0 0 0 3 9 9 9 5 0 0 0 3 9 9 9 9 9 9 5 0 0 0 0 0 0 4 9 9 9 3 5 1 0 1 8 9 9 9 9 9 9 9 6 0 4 9 9 9 5 0 4 9 9 9 5 5 9 9 9 3 0 3 9 9 9 6 5 9 9 9 3 0 0 8 9 9 9 6 9 9 6 0 0 3 9 9 8 1 2 9 9 9 1 0 0 8 9 9 3 0 0 5 9 9 9 9 9 8 1 0 0 1 8 9 6 8 2 0 0 0 0 0 0 4 8 1 0 0 8 5 0 0 0 0 8 6 0 0 0 8 5 0 0 0 0 0 0 0 5 6 6 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 0 0 5 8 1 0 0 0 0 0 2 5 0 0 0 0 2 8 1 0 0 0 3 6 0 0 0 8 2 0 0 0 0 2 9 1 0 0 8 2 0 0 0 0 1 8 1 0 0 5 1 0 0 0 0 0 0 2 9 2 0 0 0 0 0 0 4 5 0 0 0 2 6 0 0 0 0 0 0 5 2 2 6 0 0 0 0 0 3 6 0 0 0 4 5 0 0 0 0 7 1 0 3 5 0 0 0 0 0 0 0 0 2 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 0 2 9 1 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 5 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 5 0 0 0 0 0 0 0 0 5 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 5 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 6 0 1 8 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 9 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 2 8 1 0 0 0 0 0 5 2 0 0 0 3 8 1 0 0 0 5 6 0 0 0 4 9 6 0 0 0 4 3 0 0 0 2 9 1 0 0 0 7 1 0 0 0 0 5 2 0 0 0 1 7 1 0 0 0 0 3 5 0 0 5 6 0 0 0 3 9 6 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 5 2 0 0 0 2 9 1 0 0 0 0 4 2 0 0 0 0 0 0 0 7 6 3 0 0 0 4 9 6 0 0 2 7 8 1 0 0 0 5 2 0 0 0 5 8 1 0 0 5 6 0 0 0 0 0 7 1 0 0 0 5 5 0 0 0 4 6 0 0 0 4 8 1 0 0 0 7 1 0 0 0 4 6 0 0 0 0 4 5 0 0 0 5 9 1 0 1 7 1 0 1 5 0 0 2 6 0 0 2 6 0 0 0 0 0 5 2 0 0 4 5 0 0 0 0 0 4 5 0 0 7 1 0 0 0 0 0 0 5 2 0 3 6 0 0 0 0 3 6 0 0 0 0 7 1 0 0 0 0 7 1 0 0 0 5 2 0 0 0 3 6 0 0 0 0 0 0 0 7 1 0 0 0 0 0 1 8 1 0 0 0 0 7 1 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 2 8 2 6 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 0 4 5 0 0 0 0 0 0 0 0 0 0 0 0 0 3 6 0 0 0 0 7 1 0 0 0 2 5 0 0 0 0 4 3 0 0 0 0 7 1 0 0 0 4 2 0 0 0 1 5 0 0 0 0 0 0 2 9 2 0 0 0 0 0 0 8 2 0 0 0 1 7 1 0 0 0 0 0 5 1 2 6 0 0 0 0 0 5 1 0 0 0 0 0 0 0 0 1 5 0 0 0 5 0 0 0 0 0 0 0 0 7 2 7 1 0 0 0 0 0 0 0 5 9 9 8 1 0 0 0 1 8 9 9 9 9 9 8 1 0 0 0 0 0 0 0 4 5 0 0 0 0 0 0 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 3 0 0 0 0 0 0 0 0 4 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 1 0 3 9 2 0 0 0 0 0 0 0 0 0 0 8 8 1 0 4 8 1 0 0 0 0 0 0 0 0 0 5 6 0 0 0 8 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 4 9 9 9 9 1 0 0 0 0 2 6 0 5 9 9 8 1 0 0 0 0 0 2 9 9 9 8 3 5 0 0 0 0 8 9 9 8 1 5 2 0 0 0 2 9 9 9 9 3 0 0 0 0 0 8 9 9 9 9 9 9 1 0 0 0 1 8 9 9 5 1 8 9 2 0 0 7 1 5 9 9 5 0 0 0 0 0 5 9 9 9 1 0 0 0 0 0 4 9 9 9 9 9 5 0 0 0 0 0 0 7 1 0 0 8 9 9 6 0 0 0 0 0 7 1 0 0 0 1 8 9 3 8 9 2 3 9 9 2 0 1 8 9 2 3 9 9 6 0 0 0 0 0 0 4 9 9 9 5 0 0 0 4 9 6 0 8 9 9 8 1 0 0 0 0 0 8 9 9 8 1 5 9 5 0 1 8 9 6 0 2 9 9 2 0 0 0 0 5 9 9 9 2 7 1 0 0 8 9 9 9 9 9 9 8 1 0 2 9 9 1 0 0 8 9 8 1 0 4 9 9 9 5 0 2 9 9 9 9 9 9 9 3 0 0 0 3 9 9 6 0 5 9 9 5 0 3 9 9 6 0 1 8 9 9 1 0 0 4 9 9 5 0 0 5 9 9 9 9 9 9 1 0 0 0 0 3 6 0 5 3 0 0 0 0 0 5 2 0 0 0 0 5 2 0 0 3 5 0 0 0 0 0 3 6 0 0 0 4 3 0 0 0 0 2 6 0 0 0 7 1 0 0 0 0 5 2 0 0 0 1 7 1 0 0 0 0 3 5 0 3 5 0 0 0 0 0 2 6 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 5 2 0 0 3 8 1 0 0 0 0 0 4 2 0 0 0 0 0 0 0 7 3 8 1 0 0 8 4 6 0 0 2 6 4 5 0 0 0 5 2 0 0 3 6 0 0 0 0 0 5 5 0 0 0 0 7 1 0 0 0 0 7 1 0 3 6 0 0 0 0 0 3 5 0 0 0 7 1 0 0 0 0 7 1 0 0 0 7 1 0 0 0 0 7 1 0 1 7 1 0 1 5 0 0 2 6 0 0 2 6 0 0 0 0 0 5 2 0 0 1 7 1 0 0 0 0 8 2 0 0 5 2 0 0 0 0 0 0 5 1 0 0 4 5 0 0 1 8 1 0 0 0 0 2 6 0 0 0 4 3 0 0 0 0 5 2 0 0 1 8 1 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 5 2 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 5 3 2 6 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 3 0 0 0 0 7 1 0 0 0 2 6 0 0 0 0 5 2 0 0 0 0 4 2 0 0 0 7 1 0 0 0 0 7 1 0 0 0 0 0 1 8 2 0 0 0 0 0 1 7 1 0 1 8 9 8 1 0 0 0 0 0 7 1 3 5 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 3 5 0 0 0 0 0 0 0 4 3 0 3 6 0 0 0 0 0 0 4 6 0 0 0 0 0 0 0 0 0 0 2 9 2 0 0 0 0 0 0 0 0 0 1 8 2 0 0 0 0 0 0 1 8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 0 4 9 6 0 0 0 0 0 0 0 0 3 9 2 0 0 0 0 0 0 0 0 3 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 5 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 5 9 6 0 0 0 0 0 0 8 5 0 0 5 6 0 0 0 0 0 0 0 0 0 5 9 3 0 0 0 0 2 9 8 1 0 0 0 0 0 0 0 5 2 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 1 8 2 0 0 0 2 7 8 3 0 0 3 9 1 0 0 0 4 8 1 0 0 2 9 5 0 0 1 8 2 0 0 2 9 6 2 0 0 2 9 2 0 0 0 5 5 0 0 0 0 0 0 7 1 0 0 0 0 0 0 1 8 1 0 0 4 7 7 1 0 0 0 7 6 3 0 0 5 5 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 7 1 0 0 3 5 0 0 0 0 0 0 0 7 1 0 0 0 0 0 5 9 1 2 9 8 1 2 6 0 0 0 5 4 6 0 0 5 5 0 0 0 0 5 5 0 0 0 4 6 0 0 0 2 7 8 2 0 0 2 9 1 0 0 1 8 2 0 0 2 8 6 2 0 0 0 0 2 6 3 8 1 1 7 1 0 0 5 5 0 0 0 8 8 1 0 0 0 0 7 1 0 0 0 0 0 0 0 0 7 1 0 0 0 1 7 1 0 0 1 7 1 0 0 0 0 5 3 0 0 8 2 0 0 0 0 0 1 7 1 0 0 5 3 0 0 0 3 6 0 0 0 1 7 1 0 0 0 0 4 5 0 0 0 5 2 0 0 0 3 5 0 0 0 0 0 5 3 0 3 6 0 0 0 0 0 5 2 0 0 0 0 5 2 0 0 5 1 0 0 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 7 1 0 0 7 1 0 3 3 0 0 0 0 0 0 1 7 1 0 4 2 0 0 0 0 5 2 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 5 2 0 3 6 0 0 0 0 0 0 0 4 2 0 0 0 0 0 0 0 7 2 8 2 0 3 6 2 6 0 0 2 6 1 8 1 0 0 5 2 0 0 8 2 0 0 0 0 0 1 7 1 0 0 0 7 1 0 0 0 0 7 1 0 5 1 0 0 0 0 0 0 7 1 0 0 7 1 0 0 0 0 7 1 0 0 0 7 1 0 0 0 0 0 0 0 0 7 1 0 1 5 0 0 2 5 0 0 2 6 0 0 0 0 0 5 2 0 0 0 5 2 0 0 0 2 8 1 0 0 4 2 0 2 9 5 0 0 7 1 0 0 0 5 3 0 8 2 0 0 0 0 0 0 4 3 0 2 6 0 0 0 0 0 5 2 0 0 5 3 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 5 5 0 0 0 0 0 2 6 0 2 6 0 0 0 0 0 3 6 8 9 9 6 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 4 3 0 0 0 4 3 0 0 0 0 5 3 0 0 0 0 5 3 0 0 1 7 1 0 0 0 0 5 2 0 0 0 0 0 1 8 1 0 0 0 0 0 1 7 1 1 8 3 1 7 1 0 0 0 8 9 9 9 9 9 9 6 0 0 0 4 2 0 0 0 0 0 0 0 0 0 1 8 9 8 1 0 0 0 0 0 0 3 6 0 0 0 5 3 0 0 0 0 0 4 2 0 0 0 0 0 0 0 0 0 0 2 9 6 0 0 0 0 0 0 0 0 0 3 8 1 0 0 0 0 0 0 0 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 9 9 9 9 9 9 9 9 2 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 3 6 0 0 0 0 0 0 0 0 0 0 8 9 9 1 0 0 0 0 0 0 0 3 9 1 0 0 0 0 0 0 0 0 3 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 8 9 9 1 0 0 0 0 2 9 1 0 2 9 1 0 0 0 0 0 0 0 4 9 3 0 0 0 0 0 0 0 0 3 9 6 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 3 5 0 0 0 2 9 3 0 0 0 0 3 5 0 0 1 7 1 0 0 0 0 3 5 0 0 4 2 0 0 0 0 2 9 2 0 0 5 2 0 0 0 0 0 5 1 0 0 0 0 0 7 1 0 0 0 0 0 0 5 2 0 0 0 0 4 8 1 0 0 0 8 5 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 7 1 0 4 6 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 5 3 0 0 8 2 0 1 7 1 0 0 5 8 1 0 0 1 7 1 0 0 3 5 0 0 0 0 0 4 3 0 0 2 9 2 0 0 0 0 2 5 0 0 5 2 0 0 0 0 2 9 2 0 0 0 0 2 9 6 0 0 0 0 0 0 0 5 3 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 7 1 0 0 0 1 7 1 0 0 0 5 3 0 0 0 1 7 1 0 0 4 3 0 0 0 0 0 2 5 0 0 0 0 5 3 0 3 6 0 0 0 0 0 4 3 0 0 0 0 8 2 0 0 0 0 0 0 0 2 6 0 0 0 0 0 1 8 1 0 1 8 1 0 0 0 0 5 2 0 0 0 4 6 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 4 3 0 0 7 1 0 3 3 0 0 0 0 0 0 1 7 1 0 4 3 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 5 2 4 6 0 0 0 0 0 0 0 0 4 2 0 0 0 0 0 0 0 7 1 3 5 0 5 3 2 6 0 0 2 6 0 4 5 0 0 5 2 0 1 7 1 0 0 0 0 0 0 5 1 0 0 0 7 1 0 0 0 1 7 1 1 7 1 0 0 0 0 0 0 5 2 0 0 7 1 0 0 0 4 5 0 0 0 0 4 8 1 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 2 6 0 0 0 0 0 5 2 0 0 0 3 6 0 0 0 4 3 0 0 0 4 3 0 4 5 7 1 1 7 1 0 0 0 0 7 6 3 0 0 0 0 0 0 0 0 7 2 8 1 0 0 0 0 0 0 0 0 4 5 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 2 9 9 5 0 0 0 0 0 0 5 2 0 2 6 0 0 0 0 0 3 9 2 0 0 3 5 0 0 0 0 4 3 3 9 9 9 1 0 0 0 0 0 0 0 0 3 6 0 0 0 0 0 0 5 9 9 9 5 0 0 0 0 0 2 9 1 0 1 8 6 3 0 0 1 7 1 0 0 0 0 5 2 0 0 0 0 0 1 8 1 0 0 0 0 0 1 7 1 3 5 0 1 7 1 0 0 0 0 1 7 1 3 3 0 0 0 0 0 0 8 9 9 1 0 0 0 0 0 0 0 0 0 4 9 9 9 3 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 1 8 2 8 2 0 0 0 0 0 0 0 0 4 6 0 0 0 0 0 0 0 0 4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 4 9 6 0 0 0 0 0 0 0 0 2 9 1 0 0 0 0 0 0 0 0 2 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 4 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 5 0 0 0 0 0 0 0 0 0 0 0 0 4 9 5 0 0 0 0 0 0 0 0 0 2 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 3 5 0 0 0 2 8 1 0 0 0 0 1 7 1 0 3 5 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 5 2 0 1 7 1 0 0 0 0 0 3 3 0 0 0 0 0 7 1 0 0 0 0 0 1 7 1 0 0 0 0 2 8 1 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 7 1 5 5 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 5 2 0 0 7 1 0 1 7 1 0 0 5 2 0 0 0 1 7 1 0 0 5 2 0 0 0 0 0 2 6 0 0 2 8 1 0 0 0 0 1 7 1 1 7 1 0 0 0 0 0 5 2 0 0 0 0 2 8 1 0 0 0 0 0 0 0 0 5 9 9 9 5 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 7 1 0 0 0 1 7 1 0 0 0 3 6 0 0 0 4 5 0 0 0 3 5 0 2 9 2 0 4 3 0 0 0 0 0 5 6 6 0 0 0 0 0 0 2 8 1 0 0 3 6 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 3 6 0 0 0 5 3 0 0 0 0 5 9 9 9 9 9 1 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 4 3 0 0 8 9 9 9 3 0 0 0 0 0 0 1 8 9 9 9 3 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 8 9 9 9 9 9 8 1 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 5 9 5 3 9 1 0 0 0 0 0 0 4 2 0 0 0 0 0 0 0 7 1 1 8 2 8 1 2 6 0 0 2 6 0 1 8 1 0 5 2 0 1 7 1 0 0 0 0 0 0 5 2 0 0 0 7 1 0 0 0 5 3 0 1 7 1 0 0 0 0 0 0 5 2 0 0 8 9 9 9 9 5 0 0 0 0 0 0 2 9 9 9 5 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 2 6 0 0 0 0 0 5 2 0 0 0 1 8 1 0 0 7 1 0 0 0 3 5 0 5 1 5 2 1 5 0 0 0 0 0 3 8 1 0 0 0 0 0 0 0 0 2 9 3 0 0 0 0 0 0 0 0 2 8 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 0 8 3 0 0 0 0 3 6 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 5 6 9 1 0 1 8 1 0 0 0 0 0 0 0 5 3 0 0 0 0 0 2 3 0 0 0 4 1 0 0 0 0 0 2 9 9 9 2 3 5 0 0 1 7 1 0 0 0 0 5 2 0 0 0 0 0 1 8 1 0 0 0 0 0 1 7 1 3 5 0 1 7 1 0 0 0 0 1 7 1 4 3 0 0 0 0 0 0 0 0 1 8 9 2 0 0 0 2 9 9 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 9 3 0 0 0 0 0 0 0 0 8 2 0 3 8 1 0 0 0 0 0 0 0 5 5 0 0 0 0 0 0 0 0 3 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 9 9 9 9 9 9 9 1 0 3 9 9 9 9 9 9 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 3 0 0 0 0 0 0 4 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 2 9 9 9 9 9 5 0 0 0 2 6 0 0 0 0 0 0 7 1 0 4 3 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 5 2 0 1 8 9 9 9 9 9 9 9 5 0 0 0 0 0 7 1 0 0 0 0 0 1 7 1 0 0 0 0 1 7 1 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 8 9 6 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 5 2 0 0 7 1 0 1 7 1 0 0 5 2 0 0 0 1 7 1 0 0 7 1 0 0 0 0 0 2 6 0 0 2 6 0 0 0 0 0 0 7 1 1 7 1 0 0 0 0 0 5 2 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 4 5 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 8 2 0 0 7 1 0 0 0 2 6 0 4 9 5 0 5 2 0 0 0 0 0 2 9 2 0 0 0 0 0 0 0 5 3 0 0 5 2 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 8 2 0 0 0 3 6 0 0 0 0 5 2 0 0 0 1 8 5 0 0 7 1 0 0 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 4 2 0 0 7 1 0 3 3 0 0 0 0 0 0 1 7 1 0 4 3 0 0 0 0 7 1 0 0 3 9 9 9 9 3 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 4 3 0 0 0 0 3 5 0 0 0 0 5 3 0 0 2 6 0 0 0 0 0 0 4 2 0 0 0 0 7 1 0 7 1 0 5 6 6 0 2 6 0 0 2 6 0 0 4 5 0 5 2 0 1 7 1 0 0 0 0 0 0 5 1 0 0 0 8 9 9 9 9 3 0 0 1 7 1 0 0 0 0 0 0 5 1 0 0 7 1 0 0 5 3 0 0 0 0 0 0 0 0 0 0 4 9 1 0 0 0 0 0 1 5 0 0 0 0 0 0 2 6 0 0 0 0 0 5 2 0 0 0 0 5 3 0 2 6 0 0 0 0 2 5 1 7 1 3 3 2 6 0 0 0 0 2 8 3 6 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 8 2 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 5 2 0 0 0 5 9 1 0 0 0 1 4 0 0 0 0 0 0 0 8 2 0 0 0 0 0 7 1 0 0 0 2 5 0 0 0 0 0 0 0 0 0 0 4 3 0 0 1 7 1 0 0 0 0 5 2 0 0 0 0 0 0 7 1 0 0 0 0 0 1 7 1 2 9 1 1 7 1 0 0 3 9 9 9 9 9 9 9 3 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 5 9 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 8 2 8 1 3 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 0 0 0 0 0 0 0 0 3 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 9 9 9 9 9 9 9 9 2 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 6 0 0 0 0 0 0 0 0 3 9 1 0 0 0 0 0 0 0 0 0 0 2 9 2 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 5 0 0 0 0 0 0 0 0 0 0 0 0 4 9 5 0 0 0 0 0 0 0 8 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 3 9 1 0 0 0 3 5 0 0 0 2 8 1 0 0 0 0 1 7 1 0 3 5 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 5 2 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 1 7 1 0 0 0 0 2 8 1 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 8 2 4 5 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 5 2 0 0 7 1 0 1 7 1 0 0 5 2 0 0 0 1 7 1 0 0 5 2 0 0 0 0 0 2 6 0 0 2 8 1 0 0 0 0 1 7 1 1 7 1 0 0 0 0 0 5 2 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 2 0 0 0 0 7 1 0 0 0 0 0 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 4 5 0 3 6 0 0 0 0 0 7 1 8 4 8 1 7 1 0 0 0 0 2 8 1 8 2 0 0 0 0 0 0 2 8 1 2 8 1 0 0 0 0 0 0 5 2 0 0 0 0 0 0 2 9 9 9 9 9 9 9 1 0 0 0 5 2 0 0 0 0 0 7 1 0 5 2 0 0 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 7 1 0 0 7 1 0 3 3 0 2 6 0 0 0 1 7 1 0 4 2 0 0 0 0 5 1 0 0 0 0 0 2 6 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 4 3 0 0 0 0 3 5 0 0 0 0 5 2 0 0 0 5 2 0 0 0 0 0 4 2 0 0 0 0 7 1 0 7 1 0 2 9 2 0 2 6 0 0 2 6 0 0 0 8 2 5 2 0 0 8 2 0 0 0 0 0 1 7 1 0 0 0 7 1 0 0 0 0 0 0 0 5 1 0 0 0 0 0 1 7 1 0 0 7 1 0 0 0 5 2 0 0 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 1 5 0 0 0 0 0 0 2 6 0 0 0 0 0 5 2 0 0 0 0 3 6 0 5 3 0 0 0 0 2 6 3 5 0 2 6 3 5 0 0 0 0 7 1 0 4 5 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 4 5 0 0 0 5 2 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 1 8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 4 2 0 0 1 8 9 9 9 9 9 8 1 0 0 0 0 0 0 0 0 0 5 2 0 0 0 4 6 0 0 0 0 2 6 0 0 0 0 0 0 2 8 1 0 0 0 0 1 7 1 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 7 1 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 7 1 0 3 9 9 9 2 0 0 0 0 2 6 0 5 2 0 0 0 0 2 6 0 0 0 0 1 7 1 0 0 0 0 0 4 3 0 1 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 2 0 4 5 5 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 6 0 0 0 0 0 0 0 0 4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 1 8 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 3 0 0 0 0 0 0 0 0 3 9 6 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 5 2 0 0 0 0 3 5 0 0 0 2 9 3 0 0 0 0 3 5 0 0 2 8 1 0 0 0 0 0 0 0 0 4 2 0 0 0 0 2 9 2 0 0 4 2 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 5 2 0 0 0 0 4 8 1 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 7 1 0 5 3 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 5 2 0 0 7 1 0 1 7 1 0 0 5 2 0 0 0 1 7 1 0 0 3 5 0 0 0 0 0 4 3 0 0 2 9 2 0 0 0 0 2 5 0 0 5 2 0 0 0 0 2 9 2 0 0 0 0 2 6 0 0 0 0 0 0 0 2 6 0 0 0 0 0 5 2 0 0 0 0 7 1 0 0 0 0 0 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 1 7 1 5 3 0 0 0 0 0 5 3 8 1 7 3 6 0 0 0 0 1 8 1 0 0 8 2 0 0 0 0 0 0 5 3 5 3 0 0 0 0 0 0 5 3 0 0 0 0 0 0 0 4 5 0 0 0 0 0 5 3 0 0 0 5 2 0 0 0 0 0 7 1 0 2 6 0 0 0 0 0 1 8 1 0 0 4 3 0 0 0 0 2 6 0 0 0 7 1 0 0 0 0 2 6 0 0 0 1 7 1 0 0 0 0 0 0 0 3 5 0 0 0 0 0 2 6 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 4 3 0 0 0 0 5 3 0 0 0 0 5 2 0 0 0 2 6 0 0 0 0 0 4 2 0 0 0 0 7 1 0 7 1 0 0 0 0 0 2 6 0 0 2 6 0 0 0 3 6 5 2 0 0 3 6 0 0 0 0 0 5 5 0 0 0 0 7 1 0 0 0 0 0 0 0 3 5 0 0 0 0 0 4 5 0 0 0 7 1 0 0 0 1 7 1 0 0 3 5 0 0 0 0 0 4 3 0 0 0 0 0 1 5 0 0 0 0 0 0 1 7 1 0 0 0 0 7 1 0 0 0 0 1 8 2 8 1 0 0 0 0 1 7 5 2 0 0 7 4 3 0 0 0 5 3 0 0 0 5 3 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 3 6 0 0 0 0 5 2 0 0 0 0 0 0 7 1 0 0 0 0 0 0 2 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 1 8 1 0 0 0 1 4 0 0 0 0 0 0 4 5 0 0 0 0 0 0 7 1 0 0 0 1 5 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 4 2 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 0 2 5 0 5 2 0 0 0 0 2 9 2 0 0 0 4 5 0 0 0 0 0 0 5 1 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 2 0 0 8 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 8 1 0 0 0 0 0 0 0 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 2 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 9 9 2 0 0 0 0 0 0 0 0 4 9 6 0 0 0 0 0 0 0 3 6 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 5 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 2 0 0 0 0 1 8 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 4 5 0 0 0 5 9 5 0 0 0 2 7 8 2 0 0 3 8 1 0 0 0 4 8 1 0 0 0 4 8 1 0 1 8 2 0 0 2 9 6 2 0 0 0 8 3 0 0 0 0 4 5 0 0 0 0 0 7 1 0 0 0 0 0 0 1 8 1 0 0 3 7 7 1 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 7 1 0 0 5 2 0 0 0 0 0 0 0 7 1 0 0 0 0 0 5 2 0 0 7 1 0 1 7 1 0 0 5 2 0 0 0 1 7 1 0 0 0 5 3 0 0 0 3 6 0 0 0 2 7 8 2 0 0 1 8 1 0 0 1 7 1 0 0 2 8 6 2 0 0 0 0 2 6 0 0 0 0 0 0 0 2 9 6 0 0 0 3 6 0 0 0 0 0 5 5 0 0 0 8 6 0 0 0 5 2 0 0 1 8 8 1 0 0 0 0 0 5 5 7 1 0 0 0 0 0 4 9 5 0 4 6 5 0 0 0 1 8 1 0 0 0 1 8 1 0 0 0 0 0 2 9 8 1 0 0 0 0 0 4 3 0 0 0 0 5 2 0 0 8 2 0 0 0 0 0 3 6 0 0 0 5 2 0 0 0 0 4 6 0 0 0 3 6 0 0 0 2 9 3 0 0 0 4 3 0 0 0 1 8 1 0 0 0 7 1 0 0 0 0 2 6 0 0 0 1 7 1 0 0 0 0 0 0 0 0 5 3 0 0 0 0 3 6 0 0 0 7 1 0 0 0 1 7 1 0 0 0 0 0 0 7 1 0 0 0 0 0 1 8 3 0 0 3 8 1 0 0 0 0 5 2 0 0 0 0 7 1 0 0 0 0 4 2 0 0 0 0 7 1 0 7 1 0 0 0 0 0 2 6 0 0 2 6 0 0 0 0 8 9 2 0 0 0 5 8 1 0 0 5 6 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 4 6 0 0 0 5 8 1 0 0 0 7 1 0 0 0 0 3 5 0 0 3 9 3 0 0 0 1 7 1 0 0 0 0 0 1 5 0 0 0 0 0 0 0 4 6 0 0 0 4 5 0 0 0 0 0 0 5 6 6 0 0 0 0 0 0 8 9 1 0 0 4 9 3 0 0 4 5 0 0 0 0 1 8 1 0 0 0 0 0 0 7 1 0 0 0 0 0 1 8 1 0 0 0 0 5 2 0 0 0 0 0 0 7 1 0 0 0 0 0 2 8 1 0 0 0 0 5 2 0 0 2 9 2 0 0 0 8 5 0 0 0 0 0 0 0 0 2 6 0 0 0 0 2 9 2 0 0 1 8 3 0 0 0 0 0 3 6 0 0 1 8 1 0 0 0 0 0 0 8 2 0 0 0 0 0 0 3 5 0 0 0 5 2 0 0 0 0 0 0 0 0 4 6 0 0 0 0 0 0 7 1 0 0 5 1 0 0 0 0 0 0 3 9 5 0 0 0 0 0 0 4 5 0 0 0 0 0 0 0 0 0 0 3 5 0 7 1 0 0 0 0 2 7 8 9 9 9 5 0 0 0 0 0 0 0 4 3 0 1 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 8 1 0 8 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 2 0 0 0 0 0 0 1 8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 4 5 0 0 0 0 0 0 4 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 6 0 0 0 0 0 0 0 0 0 8 9 9 1 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 8 9 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 6 0 0 8 9 1 0 0 0 0 0 0 0 0 0 0 0 4 9 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 5 9 9 9 3 3 9 8 1 4 9 6 0 8 9 9 8 1 0 0 0 0 0 3 9 9 9 9 5 0 0 0 0 0 8 9 9 8 1 5 9 5 0 0 0 5 9 9 9 9 5 0 0 0 1 8 9 9 9 9 9 8 1 0 0 0 1 8 9 9 6 1 7 1 0 2 9 9 9 3 0 3 9 9 9 2 0 3 9 9 9 9 9 9 9 3 0 0 0 0 0 0 0 3 5 0 0 0 0 1 8 9 1 0 2 9 9 9 2 0 3 9 9 9 9 9 9 9 3 1 8 9 9 5 0 8 9 3 1 8 9 3 8 9 9 5 0 3 9 9 9 2 0 0 0 5 9 9 9 6 0 0 0 0 2 6 0 8 9 9 9 1 0 0 0 0 2 9 9 9 8 1 5 2 0 0 4 9 9 9 9 9 9 5 0 0 0 2 6 3 9 9 9 6 0 0 0 0 0 0 0 8 9 9 9 3 0 0 0 0 1 8 9 9 9 2 8 9 2 0 0 0 0 2 9 5 0 0 0 0 0 0 3 9 2 0 2 9 2 0 0 1 8 9 9 3 0 3 9 9 9 2 0 0 0 0 0 5 5 0 0 0 0 0 0 8 9 9 9 9 9 9 2 2 9 9 9 9 1 0 1 8 9 9 9 4 9 9 9 9 9 9 9 6 0 0 0 0 0 3 9 9 9 9 1 0 0 1 8 9 9 9 9 9 9 1 0 0 2 9 9 9 9 9 9 9 9 6 0 0 3 9 9 9 9 9 1 0 0 0 0 0 0 5 9 9 9 9 6 0 0 2 9 9 9 5 0 4 9 9 9 2 0 1 8 9 9 9 9 9 9 1 0 0 0 0 8 9 9 8 1 0 0 0 2 9 9 9 9 1 0 0 4 9 6 1 8 9 9 9 9 9 9 9 9 3 9 9 9 9 1 0 2 9 9 9 8 5 9 9 9 6 0 0 3 9 2 0 0 0 0 3 9 9 9 5 0 0 0 0 3 9 9 9 9 9 1 0 0 0 0 0 0 3 9 9 9 5 0 0 0 3 9 9 9 9 1 0 0 0 8 8 1 3 5 5 9 9 9 9 1 0 0 0 0 3 9 9 9 9 9 1 0 0 0 0 0 3 9 9 9 5 0 0 0 0 0 0 0 3 9 3 0 0 0 0 0 0 8 6 0 0 0 3 9 2 0 4 9 9 9 1 0 2 9 9 9 2 0 0 2 9 9 9 9 9 3 0 0 0 1 8 9 9 9 9 9 9 2 0 0 0 8 9 9 9 9 9 9 2 0 0 3 9 9 9 9 9 9 9 2 0 0 0 0 8 9 9 9 3 0 0 0 0 0 0 0 2 9 9 9 6 0 0 0 0 0 8 9 9 9 3 0 0 0 0 0 0 0 3 9 9 9 1 0 0 0 0 0 0 2 8 1 0 0 0 0 0 0 0 4 9 9 9 3 0 0 0 0 0 4 9 9 9 5 0 0 0 0 0 0 0 2 9 9 9 3 0 0 0 0 0 0 0 3 9 5 0 0 0 0 0 0 1 8 3 0 0 4 5 0 0 0 0 0 3 5 0 7 1 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 5 9 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 9 3 4 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 5 0 0 0 0 0 0 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 8 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 2 0 0 0 0 0 0 0 0 0 4 9 6 0 0 0 0 0 0 3 6 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 4 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 9 9 6 0 0 0 0 0 0 4 3 0 7 1 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 4 5 0 0 0 0 2 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 2 0 0 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 1 7 1 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 9 9 3 0 8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 5 0 0 0 0 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 9 9 3 0 0 0 0 3 9 9 9 1 0 0 0 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 9 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 2 5 0 0 0 0 0 0 0 0 0 0 7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 9 2 0 0 8 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 9 5 0 0 0 0 0 0 4 9 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 9 9 9 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 9 9 9 9 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 9 9 9 9 1 0 0 0 0 0 0 0 0 0 0 1 8 9 9 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 9 9 9 9 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 8 9 9 9 9 9 9 9 9 9 9 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
character_images=double(character_images)/9;

function [x_2d,y_2d]=voxelposition_to_imageposition(x,y,z,data)
data.Mview=data.viewer_matrix;
switch (data.render_type)
    case {'slicex'}
        sizeIin=[size(data.volume_original,2) size(data.volume_original,3)];
        M=[data.Mview(1,2) data.Mview(1,3) data.Mview(1,4); data.Mview(2,2) data.Mview(2,3) data.Mview(2,4); 0 0 1];
    case {'slicey'}
        sizeIin=[size(data.volume_original,1) size(data.volume_original,3)];
        M=[data.Mview(1,1) data.Mview(1,3) data.Mview(1,4); data.Mview(2,1) data.Mview(2,3) data.Mview(2,4); 0 0 1];     % Rotate 90
    case {'slicez'}
        sizeIin=[size(data.volume_original,1) size(data.volume_original,2)];
        M=[data.Mview(1,1) data.Mview(1,2) data.Mview(1,4); data.Mview(2,1) data.Mview(2,2) data.Mview(2,4); 0 0 1];
end

switch (data.render_type)
    case {'slicex'}
        Tlocalx=y; Tlocaly=z;
    case {'slicey'}
        Tlocalx=x; Tlocaly=z;
    case {'slicez'}
        Tlocalx=x; Tlocaly=y;
end

% Calculate center of the input image
mean_in=sizeIin/2;

x_2d=zeros(1,length(Tlocalx)); y_2d=zeros(1,length(Tlocalx));

Tlocalx=Tlocalx-mean_in(1);
Tlocaly=Tlocaly-mean_in(2);

for i=1:length(x)
    vector=M*[Tlocalx(i);Tlocaly(i);1];
    x_2d(i)=vector(1);
    y_2d(i)=vector(2);
end

% Calculate center of the output image
mean_out=[data.config.ImageSizeRender data.config.ImageSizeRender]/2;

% Make center of the image coordinates 0,0
x_2d=x_2d+mean_out(1); 
y_2d=y_2d+mean_out(2);

function data=mouseposition_to_voxelposition(data)
data.Mview=data.viewer_matrix;
switch (data.render_type)
    case {'slicex'}
        sizeIin=[size(data.volume_original,2) size(data.volume_original,3)];
        M=[data.Mview(1,2) data.Mview(1,3) data.Mview(1,4); data.Mview(2,2) data.Mview(2,3) data.Mview(2,4); 0 0 1];
    case {'slicey'}
        sizeIin=[size(data.volume_original,1) size(data.volume_original,3)];
        M=[data.Mview(1,1) data.Mview(1,3) data.Mview(1,4); data.Mview(2,1) data.Mview(2,3) data.Mview(2,4); 0 0 1];     % Rotate 90
    case {'slicez'}
        sizeIin=[size(data.volume_original,1) size(data.volume_original,2)];
        M=[data.Mview(1,1) data.Mview(1,2) data.Mview(1,4); data.Mview(2,1) data.Mview(2,2) data.Mview(2,4); 0 0 1];
end
M=inv(M);

% Get the mouse position
x_2d=data.mouse_position(2);
y_2d=data.mouse_position(1); 

% To rendered image position
x_2d=x_2d*data.config.ImageSizeRender; y_2d=y_2d*data.config.ImageSizeRender;

% Calculate center of the input image
mean_in=sizeIin/2;

% Calculate center of the output image
mean_out=[data.config.ImageSizeRender data.config.ImageSizeRender]/2;

% Calculate the Transformed coordinates
x_2d=x_2d - mean_out(1); 
y_2d=y_2d - mean_out(2);

location(1)= mean_in(1) + M(1,1) * x_2d + M(1,2) *y_2d + M(1,3) * 1;
location(2)= mean_in(2) + M(2,1) * x_2d + M(2,2) *y_2d + M(2,3) * 1;

switch (data.render_type)
    case {'slicex'}
        data.VoxelLocation=[data.SliceSelected(1) location(1) location(2)];
    case {'slicey'}
        data.VoxelLocation=[location(1) data.SliceSelected(2) location(2)];
    case {'slicez'}
        data.VoxelLocation=[location(1) location(2) data.SliceSelected(3)];
end
data.VoxelLocation=round(data.VoxelLocation);

data.VoxelLocation(data.VoxelLocation<1)=1;
if(data.VoxelLocation(1)>size(data.volume_original,1)), data.VoxelLocation(1)=size(data.volume_original,1); end
if(data.VoxelLocation(2)>size(data.volume_original,2)), data.VoxelLocation(2)=size(data.volume_original,2); end
if(data.VoxelLocation(3)>size(data.volume_original,3)), data.VoxelLocation(3)=size(data.volume_original,3); end

% --- Executes on mouse motion over figure - except title and menu.
function brightness_contrast_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
if(~isempty(data.handle_contrast)&&ishandle(data.handle_contrast))
    handles_contrast=guidata(data.handle_contrast);
    contrast=get(handles_contrast.slider_contrast,'value');
    brightness=get(handles_contrast.slider_brightness,'value');
    autocontrast=get(handles_contrast.checkbox_autocontrast,'value');
    if(contrast~=data.contrast||brightness~=data.brightness||autocontrast~=data.autocontrast)
        data.contrast=contrast; data.brightness=brightness; data.autocontrast=autocontrast;
        setMyData(data);
        show3d(false,false);
    end
end

function brightness_contrast_pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
handles_contrast=guidata(data.handle_contrast);
set(handles_contrast.slider_contrast,'value',0);
set(handles_contrast.slider_brightness,'value',0);
set(handles_contrast.checkbox_autocontrast,'value',0);
data.contrast=0; data.brightness=0; data.autocontrast=false;
setMyData(data);
show3d(false,false);


% --------------------------------------------------------------------
function menu_config_slicescolor_Callback(hObject, eventdata, handles)
% hObject    handle to menu_config_slicescolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData();
if(data.config.ColorSlice)
    data.config.ColorSlice=false;
else
    data.config.ColorSlice=true;
end    
setMyData(data);
set_menu_checks(data);
show3d(false,true);


% --------------------------------------------------------------------
function menu_measure_landmark_Callback(hObject, eventdata, handles)
% hObject    handle to menu_measure_landmark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.mouse_button='select_landmark';
data.measure_landmark=true;
setMyData(data);
set_mouse_shape('select_landmark',data)
