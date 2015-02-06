function varargout = tom_av2_xmipp_ang_sor_cl(varargin)
% TOM_AV2_XMIPP_ANG_SOR_CL M-file for tom_av2_xmipp_ang_sor_cl.fig
%      TOM_AV2_XMIPP_ANG_SOR_CL, by itself, creates a new TOM_AV2_XMIPP_ANG_SOR_CL or raises the existing
%      singleton*.
%
%      H = TOM_AV2_XMIPP_ANG_SOR_CL returns the handle to a new TOM_AV2_XMIPP_ANG_SOR_CL or the handle to
%      the existing singleton*.
%
%      TOM_AV2_XMIPP_ANG_SOR_CL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_AV2_XMIPP_ANG_SOR_CL.M with the given input arguments.
%
%      TOM_AV2_XMIPP_ANG_SOR_CL('Property','Value',...) creates a new TOM_AV2_XMIPP_ANG_SOR_CL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_av2_xmipp_ang_sor_cl_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_av2_xmipp_ang_sor_cl_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_av2_xmipp_ang_sor_cl

% Last Modified by GUIDE v2.5 13-Aug-2009 18:54:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_av2_xmipp_ang_sor_cl_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_xmipp_ang_sor_cl_OutputFcn, ...
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


% --- Executes just before tom_av2_xmipp_ang_sor_cl is made visible.
function tom_av2_xmipp_ang_sor_cl_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_av2_xmipp_ang_sor_cl (see VARARGIN)

% Choose default command line output for tom_av2_xmipp_ang_sor_cl


global storage_av2_xmipp_ang_cl_sort;

storage_av2_xmipp_ang_cl_sort=[];

tmp_idx=varargin{1};
tmp_idx=tmp_idx{1};
thumbs=varargin{2};
thumbs=thumbs{1};

filt=varargin{3};
filt=filt{1};

set(handles.filt,'String',filt);

storage_av2_xmipp_ang_cl_sort.indx=tmp_idx;

tmp_im=tom_spiderread(tmp_idx.filename{1});
sz=size(tmp_im.Value);

tmp_stack=zeros(sz(1),sz(2),length(tmp_idx.filename));

alg_flag=1;

for i=1:length(tmp_idx.filename)
    tmp_im=tom_spiderread(tmp_idx.filename{i});
    if (tmp_idx.flip(i)==1)
        tmp_im.Value=tom_mirror(tmp_im.Value,'x');
    end;
    if (alg_flag==1)
        tmp_stack(:,:,i)=tom_norm(tom_rotate(tom_shift(tmp_im.Value,[tmp_idx.shift(i,1) tmp_idx.shift(i,2)]),tmp_idx.rot(i)),'mean0+1std');
        if (filt>1)
            tmp_stack2(:,:,i)=tom_norm(tom_filter(tmp_stack(:,:,i),filt),'mean0+1std');
        else
            tmp_stack2(:,:,i)=tom_norm(tmp_stack(:,:,i),'mean0+1std');
        end;
        %tmp_stack(:,:,i)=tom_filter(tmp_stack(:,:,i),2);
    else
        tmp_stack(:,:,i)=tom_norm(tom_filter(tmp_im.Value,2),'mean0+1std'); 
    end;
end;
storage_av2_xmipp_ang_cl_sort.tmp_stack=tmp_stack;

storage_av2_xmipp_ang_cl_sort.sz=sz;

axes(handles.m_ax1);
[gg lookup]=tom_gallery(tmp_stack2);
imhandle=imagesc(gg'); axis image; colormap gray; set(handles.m_ax1,'XTick',[]);
set(handles.m_ax1,'YTick',[]);
set(imhandle, 'Tag', 'main_image_sort');
storage_av2_xmipp_ang_cl_sort.lookup=lookup;

hold on;
for i=1:length(storage_av2_xmipp_ang_cl_sort.indx.sel)
    rect_cent=lookup(i,:);
    rect_cent=rect_cent-(sz(1)./2.5);
    %plot(rect_cent(1),rect_cent(2),'gs','MarkerSize',round(sz(1)./8),'MarkerFaceColor','g'); 
    if ( storage_av2_xmipp_ang_cl_sort.indx.sel(i)==1)
        plot(rect_cent(1),rect_cent(2),'gs','MarkerSize',round(sz(1)./8),'MarkerFaceColor','g');
    else
        plot(rect_cent(1),rect_cent(2),'rs','MarkerSize',round(sz(1)./8),'MarkerFaceColor','r');
    end;

end;
hold off;     


axes(handles.axx_proj);
imagesc(thumbs(:,:,1)'); axis image; colormap gray; 
set(handles.m_ax1,'XTick',[]);
set(handles.m_ax1,'YTick',[]);

axes(handles.axx_mean);
imagesc(thumbs(:,:,2)'); axis image; colormap gray; 
set(handles.m_ax1,'XTick',[]);
set(handles.m_ax1,'YTick',[]);

axes(handles.axx_var);
imagesc(thumbs(:,:,3)'); axis image; colormap gray; 
set(handles.m_ax1,'XTick',[]);
set(handles.m_ax1,'YTick',[]);


handles.output = hObject;

set(handles.output,'Pointer','crosshair');
set(findobj('Tag','main_image_sort'),'buttonDownFcn',@pick_particle_sort);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_av2_xmipp_ang_sor_cl wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_av2_xmipp_ang_sor_cl_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
close(handles.figure1);

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global storage_av2_xmipp_ang_cl_sort;

handles.output=storage_av2_xmipp_ang_cl_sort.indx.sel;
guidata(hObject, handles);

uiresume;

% --- Executes on button press in invert_it_man.
function invert_it_man_Callback(hObject, eventdata, handles)
% hObject    handle to invert_it_man (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global storage_av2_xmipp_ang_cl_sort;


try
    im=tom_spiderread(storage_av2_xmipp_ang_cl_sort.indx(1).filename{1});
catch
    im=tom_spiderread(storage_av2_xmipp_ang_cl_sort.indx(1).filename{2});
end;
sz_im=size(im.Value);


storage_av2_xmipp_ang_cl_sort.indx.sel=(storage_av2_xmipp_ang_cl_sort.indx.sel==0);

lookup=storage_av2_xmipp_ang_cl_sort.lookup;


axes(handles.m_ax1);
hold on;
for i=1:length(storage_av2_xmipp_ang_cl_sort.indx.sel)
    rect_cent=lookup(i,:);
    rect_cent=rect_cent-(sz_im(1)./2.5);
    if ( storage_av2_xmipp_ang_cl_sort.indx.sel(i)==1)
        plot(rect_cent(1),rect_cent(2),'gs','MarkerSize',round(sz_im(1)./8),'MarkerFaceColor','g');
    else
        plot(rect_cent(1),rect_cent(2),'rs','MarkerSize',round(sz_im(1)./8),'MarkerFaceColor','r');
    end;
end;
hold off;     


% --- Executes on button press in recalc.
function recalc_Callback(hObject, eventdata, handles)
% hObject    handle to recalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global storage_av2_xmipp_ang_cl_sort;


tmp_stack=storage_av2_xmipp_ang_cl_sort.tmp_stack(:,:,find(storage_av2_xmipp_ang_cl_sort.indx.sel));

[var_st tmp_mean]=tom_calc_variance2d(tmp_stack,'default','default','para');



axes(handles.axx_mean);
imagesc(tmp_mean'); axis image; colormap gray; 
set(handles.m_ax1,'XTick',[]);
set(handles.m_ax1,'YTick',[]);


axes(handles.axx_var);
imagesc(var_st'); axis image; colormap gray; 
set(handles.m_ax1,'XTick',[]);
set(handles.m_ax1,'YTick',[]);


function pick_particle_sort(a,b)

global storage_av2_xmipp_ang_cl_sort;

point1=get(gca,'currentpoint');
button = get(gcf,'selectiontype');



try
    im=tom_spiderread(storage_av2_xmipp_ang_cl_sort.indx(1).filename{1});
catch
    im=tom_spiderread(storage_av2_xmipp_ang_cl_sort.indx(10).filename{1});
end;

sz_im=size(im.Value);

%button values:
%normal: left mouse button
%alt: right mouse button
%extend: middle mouse buttons


if (strcmp(button,'normal') == true)
    
    [pointidx, pointcoords, distance] = tom_nearestpoint(point1(1,1:2),storage_av2_xmipp_ang_cl_sort.lookup);
    %hold on; plot(pointcoords(1),pointcoords(2),'ro'); hold off;
    
    ind=pointidx;
    rect_cent=pointcoords;
    
    rect_cent=rect_cent-(sz_im./2.5);
    
    storage_av2_xmipp_ang_cl_sort.indx.sel(ind)=storage_av2_xmipp_ang_cl_sort.indx.sel(ind)==0;
    
    hold on;
    if ( storage_av2_xmipp_ang_cl_sort.indx.sel(ind)==1)
        plot(rect_cent(1),rect_cent(2),'gs','MarkerSize',round(sz_im(1)./8),'MarkerFaceColor','g');
    else
        plot(rect_cent(1),rect_cent(2),'rs','MarkerSize',round(sz_im(1)./8),'MarkerFaceColor','r');
    end;
    hold off;
    
end;



function filt_Callback(hObject, eventdata, handles)
% hObject    handle to filt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filt as text
%        str2double(get(hObject,'String')) returns contents of filt as a double

global storage_av2_xmipp_ang_cl_sort;

sz=storage_av2_xmipp_ang_cl_sort.sz;

filt=str2double(get(handles.filt,'String'));

tmp_stack2=storage_av2_xmipp_ang_cl_sort.tmp_stack;

if (filt>1)
    parfor i=1:size(storage_av2_xmipp_ang_cl_sort.tmp_stack,3)
        tmp_stack2(:,:,i)=tom_norm(tom_filter(storage_av2_xmipp_ang_cl_sort.tmp_stack(:,:,i),filt),'mean0+1std');
    end;
end;

axes(handles.m_ax1);
[gg lookup]=tom_gallery(tmp_stack2);
imhandle=imagesc(gg'); axis image; colormap gray; set(handles.m_ax1,'XTick',[]);
set(handles.m_ax1,'YTick',[]);
set(imhandle, 'Tag', 'main_image_sort');
storage_av2_xmipp_ang_cl_sort.lookup=lookup;

hold on;
for i=1:length(storage_av2_xmipp_ang_cl_sort.indx.sel)
    rect_cent=lookup(i,:);
    rect_cent=rect_cent-(sz(1)./2.5);
    %plot(rect_cent(1),rect_cent(2),'gs','MarkerSize',round(sz(1)./8),'MarkerFaceColor','g'); 
    if ( storage_av2_xmipp_ang_cl_sort.indx.sel(i)==1)
        plot(rect_cent(1),rect_cent(2),'gs','MarkerSize',round(sz(1)./8),'MarkerFaceColor','g');
    else
        plot(rect_cent(1),rect_cent(2),'rs','MarkerSize',round(sz(1)./8),'MarkerFaceColor','r');
    end;

end;
hold off;     


set(handles.output,'Pointer','crosshair');
set(findobj('Tag','main_image_sort'),'buttonDownFcn',@pick_particle_sort);




% --- Executes during object creation, after setting all properties.
function filt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

filt=str2double(get(handles.filt,'String'));

if (eventdata.VerticalScrollCount==1)
    filt=filt+1;
else
    filt=filt-1;
end;

if (filt<1)
    filt=1;
end;

set(handles.filt,'String',num2str(filt));

filt_Callback(hObject, eventdata, handles);
