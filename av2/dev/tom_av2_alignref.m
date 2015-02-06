function varargout = tom_av2_alignref(varargin)
%GUI for centering and rotating particles in a reference stack 
%
%SYNTAX
%tom_av2_alignref(stackfile)
%
%
%SEE ALSO
%
%
%Copyright (c) 2006
%TOM toolbox for Electron Tomography
%Max-Planck-Institute for Biochemistry
%Dept. Molecular Structural Biology
%82152 Martinsried, Germany
%http://www.biochem.mpg.de/tom
%
%Created: 13/02/06 AK

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_av2_alignref_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_alignref_OutputFcn, ...
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
%%  Opening Function                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tom_av2_alignref_OpeningFcn(hObject, eventdata, handles, varargin)

handles.filename = varargin{1};

try
    handles.stack = tom_emreadc(handles.filename);
catch
    error('Could not open stack!');
end

if size(handles.stack.Value,3) > 1
    set(handles.alignref_slider,'Value',0);
    set(handles.alignref_slider,'Max',size(handles.stack.Value,3));
    set(handles.alignref_slider,'SliderStep',[1./size(handles.stack.Value,3) 1./size(handles.stack.Value,3)]);
else
    set(handles.alignref_slider,'Visible','off');
end
handles.no = 1;
handles = display_particle(handles);

handles.output = hObject;
guidata(hObject, handles);
uiwait(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_av2_alignref_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
close(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Rotate                                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_alignref_rotate_Callback(hObject, eventdata, handles)

handles.stack.Value(:,:,handles.no) = tom_rotate(handles.stack.Value(:,:,handles.no),str2num(get(handles.alignref_rotateval,'String')));
guidata(hObject, handles);
handles = display_particle(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Shift right                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_alignref_right_Callback(hObject, eventdata, handles)

handles.stack.Value(:,:,handles.no) = tom_shift(handles.stack.Value(:,:,handles.no), [1 0]);
guidata(hObject, handles);
handles = display_particle(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Shift left                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_alignref_left_Callback(hObject, eventdata, handles)

handles.stack.Value(:,:,handles.no) = tom_shift(handles.stack.Value(:,:,handles.no), [-1 0]);
guidata(hObject, handles);
handles = display_particle(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Shift up                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_alignref_up_Callback(hObject, eventdata, handles)

handles.stack.Value(:,:,handles.no) = tom_shift(handles.stack.Value(:,:,handles.no), [0 -1]);
guidata(hObject, handles);
handles = display_particle(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Shift down                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_alignref_down_Callback(hObject, eventdata, handles)

handles.stack.Value(:,:,handles.no) = tom_shift(handles.stack.Value(:,:,handles.no), [0 1]);
guidata(hObject, handles);
handles = display_particle(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Center                                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_alignref_center_Callback(hObject, eventdata, handles)

%pick new center
tmpobj = findobj('Tag','alignref_image');
axes(tmpobj);

[x,y] = ginput(1);

%Calculate shift vector
shiftx = -(x-(floor(size(handles.stack.Value,1)./2)+1));
shifty = -(y-(floor(size(handles.stack.Value,2)./2)+1));

handles.stack.Value(:,:,handles.no) = tom_shift(handles.stack.Value(:,:,handles.no), [shiftx, shifty]);
guidata(hObject, handles);
handles = display_particle(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Filter low                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alignref_filter_low_Callback(hObject, eventdata, handles)

guidata(hObject, handles);
handles = display_particle(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Filter high                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alignref_filter_high_Callback(hObject, eventdata, handles)

guidata(hObject, handles);
handles = display_particle(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Filter enable                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alignref_filter_enable_Callback(hObject, eventdata, handles)

guidata(hObject, handles);
handles = display_particle(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Mask radius                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function aligref_maskradius_Callback(hObject, eventdata, handles)

guidata(hObject, handles);
handles = display_particle(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Mask enable                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alignref_mask_enable_Callback(hObject, eventdata, handles)

guidata(hObject, handles);
handles = display_particle(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Slider                                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alignref_slider_Callback(hObject, eventdata, handles)

handles.no = get(hObject,'Value')+1;
if handles.no > size(handles.stack.Value,3)
    handles.no = size(handles.stack.Value,3);
end

handles = display_particle(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Exit                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_alignref_exit_Callback(hObject, eventdata, handles)

tom_emwrite(handles.filename,handles.stack);
disp(['Stack saved to ' handles.filename]);
uiresume(handles.figure1);


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
%%  display particle                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = display_particle(handles)

tmpobj = findobj('Tag','alignref_image');
axes(tmpobj);

im = handles.stack.Value(:,:,handles.no);

%filter
if get(handles.alignref_filter_enable,'Value') == 1
    low = round(str2num(get(handles.alignref_filter_low,'String')));
    high = round(str2num(get(handles.alignref_filter_high,'String')));
    if isempty(low)
        errordlg('Input Filter low');
        return;
    elseif isempty(high)
        errordlg('Input Filter high');
        return;
    else
        im = tom_bandpass(im,low,high);
    end
end

%mask
if get(handles.alignref_mask_enable,'Value') == 1
    maskradius = round(str2num(get(handles.aligref_maskradius,'String')));
    if isempty(maskradius)
        errordlg('Input Mask radius');
        set(tmpobj,'Tag','alignref_image');
    end
    tom_imagesc(tom_spheremask(im,maskradius),'noinfo');
    
%no mask
else
    tom_imagesc(im,'noinfo');
end

set(tmpobj,'Tag','alignref_image');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alignref_rotateval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function alignref_filter_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function alignref_filter_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function aligref_maskradius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function alignref_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function alignref_rotateval_Callback(hObject, eventdata, handles)