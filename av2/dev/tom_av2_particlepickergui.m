function varargout = tom_av2_particlepickergui(varargin)
%GUI for picking particles from 2D EM images 
%
%SYNTAX
%tom_av2_particlepickergui
%

%SEE ALSO
%
%
%Copyright (c) 2005
%TOM toolbox for Electron Tomography
%Max-Planck-Institute for Biochemistry
%Dept. Molecular Structural Biology
%82152 Martinsried, Germany
%http://www.biochem.mpg.de/tom
%
%Created: 20/11/05 AK

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_av2_particlepickergui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_particlepickergui_OutputFcn, ...
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
function tom_av2_particlepickergui_OpeningFcn(hObject, eventdata, handles, varargin)

global storage_av2_particlepicker;

%Specify directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(varargin,1) == 0
    pathname = uigetdir('Select a directory');

    if isequal(pathname,0) 
        error('Cancel button pressed. No data loaded.');
        return; 
    end;
    storage_av2_particlepicker.path = pathname;
else
    storage_av2_particlepicker.path = varargin{1};
end

set(gcf,'Renderer','zbuffer');
axes(findobj('Tag','image'));axis off;
axes(findobj('Tag','powerspectrum'));axis off;
axes(findobj('Tag','particle_clicked'));axis off;
axes(findobj('Tag','particle_aligned'));axis off;
axes(findobj('Tag','particle_average'));axis off;
axes(findobj('Tag','particle_average_filtered'));axis off;
axes(findobj('Tag','classes_text'));axis off;

[storage_av2_particlepicker.dircell, storage_av2_particlepicker.flagcell] = get_dircontents(storage_av2_particlepicker.path,{},{});
set(findobj('Tag','images_bad'),'String', '0');
set(findobj('Tag','images_good'),'String', num2str(size(storage_av2_particlepicker.flagcell,2)));
storage_av2_particlepicker.imagenumber = 0;
storage_av2_particlepicker.classes = {'default'};
storage_av2_particlepicker.classcolors = {[0 1 0]};
storage_av2_particlepicker.classradii = {''};
storage_av2_particlepicker.classfilter = {[]};
storage_av2_particlepicker.reference_images = {};
storage_av2_particlepicker.resampleval = 1;

update_classes('render');
manage_align('create');

classes_stats();

axes(findobj('Tag','image'));
render_image(1);
set(findobj('Tag','bandpass_high'),'String',num2str(size(storage_av2_particlepicker.image,1)./2))


storage_av2_particlepicker.reference_image = struct();

handles.output = hObject;
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_av2_particlepickergui_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Histogram Set                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_histogram_set_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

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
set(gca,'Xlim',[x(1) x(2)]);

if get(findobj('Tag','histogram_emimage'),'Value') == 1
    storage_av2_particlepicker.DataScale=[x(1) x(2)];
    render_image(storage_av2_particlepicker.imagenumber);
else
    storage_av2_particlepicker.powerScale=[x(1) x(2)];
    render_ps(0,storage_av2_particlepicker.powerScale);
end
set(findobj('Tag','histogram_low'),'String',x(1));
set(findobj('Tag','histogram_high'),'String',x(2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Histogram Reset                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_histogram_reset_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

if get(findobj('Tag','histogram_emimage'),'Value') == 1
    [storage_av2_particlepicker.DataScale(1) storage_av2_particlepicker.DataScale(2)] = calculate_histogram(storage_av2_particlepicker.image);
    render_image(storage_av2_particlepicker.imagenumber);
else
    [storage_av2_particlepicker.powerScale(1) storage_av2_particlepicker.powerScale(2)] = calculate_histogram(storage_av2_particlepicker.powerspectrum);
    render_ps(0,storage_av2_particlepicker.powerScale);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Histogram Set Manually                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_histogram_setmanually_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

min=str2num(get(findobj('Tag','histogram_low'),'String'));
max=str2num(get(findobj('Tag','histogram_high'),'String'));

if get(findobj('Tag','histogram_emimage'),'Value') == 1
    if max>min
        storage_av2_particlepicker.DataScale = [min max];
        set(findobj('Tag','histogram'),'Xlim',[min max]);
    else
        set(findobj('Tag','histogram_low'),'String',num2str(storage_av2_particlepicker.DataScale(1)));
        set(findobj('Tag','histogram_high'),'String',num2str(storage_av2_particlepicker.DataScale(2)));
    end
    render_image(storage_av2_particlepicker.imagenumber);
else
    if max>min
        storage_av2_particlepicker.powerScale = [min max];
        set(findobj('Tag','histogram'),'Xlim',[min max]);
    else
        set(findobj('Tag','histogram_low'),'String',num2str(storage_av2_particlepicker.powerScale(1)));
        set(findobj('Tag','histogram_high'),'String',num2str(storage_av2_particlepicker.powerScale(2)));
    end
    render_ps(0,storage_av2_particlepicker.powerScale);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Histogram Pick Low                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_histogram_picklow_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

[x y] = ginput(1);

clickedimage = get(gca,'Tag');
if strcmp(clickedimage,'image') == 1

    if x > 0 & y > 0 &  x < size(storage_av2_particlepicker.image,2) & y < size(storage_av2_particlepicker.image,1)

        val = storage_av2_particlepicker.image(round(x),round(y));

        if val >= storage_av2_particlepicker.DataScale(2)
            errordlg('lower limit must be smaller than upper limit!','Histrogram error');
        else
            storage_av2_particlepicker.DataScale = [val storage_av2_particlepicker.DataScale(2)];
            set(findobj('Tag','histogram_low'),'String',num2str(round(val)));
            render_image(storage_av2_particlepicker.imagenumber);
        end

    end
    
elseif strcmp(clickedimage,'powerspectrum') == 1
    
    if x > 0 & y > 0 &  x < size(storage_av2_particlepicker.powerspectrum,2) & y < size(storage_av2_particlepicker.powerspectrum,1)

        val = storage_av2_particlepicker.powerspectrum(round(x),round(y));

        if val >= storage_av2_particlepicker.powerScale(2)
            errordlg('lower limit must be greater than upper limit!','Histrogram error');
        else
            storage_av2_particlepicker.powerScale = [num2str(round(val)) storage_av2_particlepicker.powerScale(2)];
            set(findobj('Tag','histogram_low'),'String',num2str(round(val)));
            render_ps(0,storage_av2_particlepicker.powerScale);
        end

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Histogram Pick High                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_histogram_pickhigh_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

[x y] = ginput(1);

clickedimage = get(gca,'Tag');
if strcmp(clickedimage,'image') == 1
    if x > 0 & y > 0 &  x < size(storage_av2_particlepicker.image,2) & y < size(storage_av2_particlepicker.image,1)

        val = storage_av2_particlepicker.image(round(x),round(y));

        if val <= storage_av2_particlepicker.DataScale(1)
            errordlg('upper limit must be greater than lower limit!','Histrogram error');
        else
            storage_av2_particlepicker.DataScale = [storage_av2_particlepicker.DataScale(1) val];
            set(findobj('Tag','histogram_high'),'String',num2str(round(val)));
            render_image(storage_av2_particlepicker.imagenumber);
        end

    end
    
elseif strcmp(clickedimage,'powerspectrum') == 1
    if x > 0 & y > 0 &  x < size(storage_av2_particlepicker.powerspectrum,2) & y < size(storage_av2_particlepicker.powerspectrum,1)

        val = storage_av2_particlepicker.powerspectrum(round(x),round(y));

        if val <= storage_av2_particlepicker.powerScale(1)
            errordlg('upper limit must be greater than lower limit!','Histrogram error');
        else
            storage_av2_particlepicker.powerScale = [storage_av2_particlepicker.powerScale(1) val];
            set(findobj('Tag','histogram_high'),'String',num2str(round(val)));
            render_ps(0,storage_av2_particlepicker.powerScale);
        end

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  histogram emimage                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function histogram_emimage_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

[min,max] = calculate_histogram(storage_av2_particlepicker.image);
storage_av2_particlepicker.DataScale = [min, max];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  histogram power spectrum                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function histogram_power_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

[min,max] = calculate_histogram(storage_av2_particlepicker.powerspectrum);
storage_av2_particlepicker.powerScale = [min, max];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Fileslider                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fileslider_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

position = round(get(hObject,'Value'));
drawnow;
render_image(position);
drawnow;
auto_save();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Reload directory                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_rescandir_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

storage_av2_particlepicker.dircell = get_dircontents(storage_av2_particlepicker.path, storage_av2_particlepicker.dircell);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Goto last image with picked particles                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_goto_lastpicked_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

try
    counter = 1;
    for i = 1:size(storage_av2_particlepicker.points,2)
       if size(storage_av2_particlepicker.points(i).points,2) > 0
        counter = i;
       end
    end    

    storage_av2_particlepicker.imagenumber = counter;
    render_image(storage_av2_particlepicker.imagenumber,1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  File input box                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filename_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

filestring = get(hObject,'String');
result = strmatch(filestring,storage_av2_particlepicker.dircell);
if ~isempty(result)
    render_image(result);
    set(findobj('Tag','fileslider'),'Value',result);
else
    errordlg('Filename not found','File error');
    set(hObject,'String',storage_av2_particlepicker.dircell{storage_av2_particlepicker.imagenumber});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Resample input box                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_resample_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

resampleval = round(str2num(get(hObject,'String')));

if isempty(resampleval) | resampleval < 1 | resampleval > 10
    errordlg('Valid resampling range: 1-10');
    set(hObject,'String',num2str(storage_av2_particlepicker.resampleval));
else
    storage_av2_particlepicker.resampleval = resampleval;
    set(hObject,'String',num2str(storage_av2_particlepicker.resampleval));
end

render_image(storage_av2_particlepicker.imagenumber,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Zoom                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_zoom_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

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
storage_av2_particlepicker.actualaxis=[x(1) x(2) y(1) y(3)];
tmpobj = findobj('Tag','image');axes(tmpobj);
axis([x(1) x(2) y(1) y(3)]);
set(tmpobj,'Tag','image');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Zoom Reset                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_zoomreset_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

dim_x = size(storage_av2_particlepicker.image,2);
dim_y = size(storage_av2_particlepicker.image,1);

storage_av2_particlepicker.actualaxis=[1 dim_x 1 dim_y];
tmpobj = findobj('Tag','image');axes(tmpobj);
axis([1 dim_x 1 dim_y]);
set(tmpobj,'Tag','image');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Bandpass Enable                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bandpass_switch_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

render_image(storage_av2_particlepicker.imagenumber);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Bandpass High                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bandpass_high_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

string = round(str2num(get(hObject,'String')));
if ~isempty(string)
    render_image(storage_av2_particlepicker.imagenumber);
else
    errordlg('Enter a number');
    set(hObject,'String',num2str(round(size(storage_av2_particlepicker.image,1)./2)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Bandpass Low                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bandpass_low_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

string = round(str2num(get(hObject,'String')));
if ~isempty(string)
    render_image(storage_av2_particlepicker.imagenumber);
else
    errordlg('Enter a number');
    set(hObject,'String','0');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Picklist load                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_picklist_load_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

[filename, pathname] = uigetfile('*.mat', 'Load Picklist');
if ischar(filename)
    s = load([pathname '/' filename]);

    if isfield(s,'align2d') ~= 1
        errordlg('This is not a picklist file.','File Error');
    end
    
    s = s.align2d;
    %generate new unique labels for each drawmark
    for i = 1:size(s,2)
        s(1,i).label = tempname;
    end
    storage_av2_particlepicker.align = s;

    %rebuild classes
    out = manage_align('getclasses');
    storage_av2_particlepicker.classes = out{1};
    storage_av2_particlepicker.classcolors = out{2};
    storage_av2_particlepicker.classradii = out{3};
    
    %restore radius input box
    set(findobj('Tag','particle_radius'),'String',num2str(storage_av2_particlepicker.classradii{1}));
    update_classes('render');
    
    manage_align('rebuildindex');

    %build the reference images for all classes
    rebuild_reference_images();
    render_avg();
    render_filtered_avg();
    classes_stats();
    render_points();    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Picklist save                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_picklist_save_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

[filename, pathname] = uiputfile('*.mat', 'Save Picklist as');
if ischar(filename)
    align2d = storage_av2_particlepicker.align;
    try
        align2d = rmfield(align2d,'label');
    end
    save([pathname '/' filename],'align2d');
    disp('Picklist Saved');
end

storage_av2_particlepicker.picklistfilename = [pathname '/' filename];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Picklist autosave                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function autosave_switch_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

if isfield(storage_av2_particlepicker,'picklistfilename') == false
    [filename, pathname] = uiputfile('*.mat', 'Save Picklist as');
    if ischar(filename)
        align2d = storage_av2_particlepicker.align;
        align2d = rmfield(align2d,'label');
        save([pathname '/' filename],'align2d');
        disp('Picklist Saved');
    end
    storage_av2_particlepicker.picklistfilename = [pathname '/' filename];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Classes dropdown box                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function menu_classes_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

contents = get(hObject,'String');
classname = contents{get(hObject,'Value')};

set(findobj('Tag','particle_radius'),'String',num2str(manage_align('getclassradius',classname)));

if isfield(storage_av2_particlepicker.reference_image,classname)
    %display new particle class average image
    radius = manage_align('getclassradius',classname);
    centerpos = radius.*2;
    im = storage_av2_particlepicker.reference_image.(classname)(round(centerpos-radius+1):round(centerpos+radius),round(centerpos-radius+1):round(centerpos+radius));
    axesobj = findobj('Tag','particle_average');
    axes(axesobj);
    [h,n]=tom_hist3d(im);
    Scale=[n(1)  n(size(n,2))]; 
    imagesc(im',Scale);
    axis ij; axis off; colormap gray;
    set(axesobj,'Tag','particle_average');

end

check_radiusbox();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Add class                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_class_add_Callback(hObject, eventdata, handles)

prompt = {'Enter name for new class:'};
dlg_title = 'Add class';
num_lines = 1;
def = {''};
answer = inputdlg(prompt,dlg_title,num_lines,def,'on');
if isempty(answer)
    disp('Aborted.');
    return;
end

if isempty(regexpi(answer{1},'^[a-z][a-z0-9_]*$'))
    errordlg('Class name must start with a letter and can only contain letters, numbers and ''_''.');
    return;
end

color = uisetcolor([0 1 0],'Select Color for class');

update_classes('add',answer{1},color);

check_radiusbox();

%activate class
set(findobj('Tag','menu_classes'),'Value',size(get(findobj('Tag','menu_classes'),'String'),1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Delete class                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_class_delete_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

if size(storage_av2_particlepicker.classes,2) == 1
    errordlg('You must have at least one class in the class list!','Delete Error');
    return;
end

contents = get(findobj('Tag','menu_classes'),'String');
classname = contents{get(findobj('Tag','menu_classes'),'Value')};
button = questdlg(sprintf('Are you sure to delete class ''%s'' and all its particles?',classname),'Confirm delete','yes','no','no');
if strcmp(button,'yes') == 1
    update_classes('delete',classname);
end

check_radiusbox();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Save particle stack                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_savestack_Callback(hObject, eventdata, handles)
global storage_av2_particlepicker;

%select classes
[selection,value] = tom_choosebox('ListString',storage_av2_particlepicker.classes,'Name','Select Classes to export','PromptString','Available Classes','SelectString','Classes to Export');

if value == 0
    disp('Aborted');
    return;
end

classes = storage_av2_particlepicker.classes(selection);
classradii = storage_av2_particlepicker.classradii(selection);

%check for different radii
string = '';
radius = [];
radius = unique(cell2mat(classradii));

%select radius if there are multiple radii
if size(radius,2) > 1
    prompt = {['Enter radius: (class radii: ', mat2str(radius), ')']};
    dlg_title = 'Input radius';
    num_lines = 1;
    def = {num2str(radius(1))};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    radius = str2num(answer{1});
end


%get filename for stack
[filename, pathname] = uiputfile('*.em', 'Save Stack as');
if ischar(filename)
   stackname = [pathname '/' filename]; 
else
    disp('Aborted');
    return;
end
[pathstr, name] = fileparts(stackname);
alignname = [pathstr '/' name '.mat'];

%align stack?
button = questdlg('Save stack aligned?','Alignment');
if strcmp(button,'Yes')
    alignedflag = 1;
elseif strcmp(button,'No')
    alignedflag = 0;
else
    disp('Aborted');
    return;
end

%norm stack?
button = questdlg('Save stack normalized?','Norm','phase','3std','none','none');
if strcmp(button,'phase')
    normflag = 1;
elseif strcmp(button,'none')
    normflag = 0;
elseif strcmp(button,'3std')
    normflag = 2;
else
    disp('Aborted');
    return;
end

create_particle_stack(classes, stackname, alignname, alignedflag, normflag, radius);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  input radius                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function particle_radius_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

idx = get(findobj('Tag','menu_classes'),'Value');
storage_av2_particlepicker.classradii{idx} = round(str2num(get(hObject,'String')));
set(hObject,'String',num2str(storage_av2_particlepicker.classradii{idx}));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Pick particle                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_particle_pick_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

radius = round(str2num(get(findobj('Tag','particle_radius'),'String')));
if isempty(radius) == 1
    errordlg('Enter radius for current class first!','Missing radius');
    set(hObject,'Value',0);
else
    if get(hObject,'Value') == 1
        set(handles.output,'Pointer','crosshair');
        set(findobj('Tag','main_image'),'buttonDownFcn',@pick_particle);
    else
        set(handles.output,'Pointer','arrow');
        set(findobj('Tag','main_image'),'buttonDownFcn',@dummy);
    end
end

check_radiusbox();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Delete particle                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_particle_delete_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

axes(findobj('Tag','image'));
[x,y]=ginput(1);

pointcoords = [];
points = manage_align('getpointsincurrentframe');

for i = points
    pointcoords = [pointcoords;storage_av2_particlepicker.align(1,i).position.x storage_av2_particlepicker.align(1,i).position.y];
end

pointidx = tom_nearestpoint([x,y], pointcoords);
pointnum = points(pointidx);
pointlabel = storage_av2_particlepicker.align(1,pointnum).label;
delete(findobj('Tag',strcat('point_',pointlabel)));
drawmark(storage_av2_particlepicker.align(1,pointnum).position.x/storage_av2_particlepicker.resampleval, storage_av2_particlepicker.align(1,pointnum).position.y/storage_av2_particlepicker.resampleval,[1 0 0],strcat('point_',storage_av2_particlepicker.align(1,pointnum).label),num2str(pointidx));

button = questdlg(sprintf('Are you sure to delete particle ''%s''?',num2str(pointidx)),'Confirm delete','yes','no','no');
if strcmp(button,'yes') == 1
    manage_align('deleteparticle',pointnum);
end

classes_stats();
check_radiusbox();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Delete all particles in image                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_deleteall_Callback(hObject, eventdata, handles)

manage_align('deletepointsincurrentframe');

classes_stats();
check_radiusbox();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Refine                                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_particle_refine_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

%Get list of EM-Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numframes = size(storage_av2_particlepicker.points,2);

numpoints = size(storage_av2_particlepicker.align,2);

contents = get(findobj('Tag','menu_classes'),'String');
val = get(findobj('Tag','menu_classes'),'Value');
classname = contents{val};

radius = manage_align('getclassradius',classname);

if manage_align('countclassmembers',classname) == 0
    errordlg('No particles in this class!');
    return;
end

h = waitbar(0,['Refining class ' classname ]);
ref_before = storage_av2_particlepicker.reference_image.(classname);
runs = str2num(get(findobj('Tag','refinenumber'),'String'));

mask = tom_spheremask(ones(size(ref_before)),size(ref_before,1)./2-size(ref_before,1)./5,size(ref_before,1)./10);
mask_cc_trans=tom_spheremask(ones(size(ref_before)),round(str2num(get(findobj('Tag','alignment_translation'),'String'))),1);

for run = 1:runs
    
    j = 1;
    ref_image = storage_av2_particlepicker.reference_image.(classname);
    storage_av2_particlepicker.reference_image.(classname) = zeros(size(ref_image));


    for i = 1:numframes
        filename = storage_av2_particlepicker.dircell{i};
        image = tom_reademheader([storage_av2_particlepicker.path '/' filename]);

        for k = 1:size(storage_av2_particlepicker.points(i).points,2)

            %add particle to stack if classname matches
            if strmatch(storage_av2_particlepicker.align(1,j).class,classname,'exact');

                x = storage_av2_particlepicker.align(1,j).position.x;
                y = storage_av2_particlepicker.align(1,j).position.y;

                %cut particles
                %taper mode
                imagesize = image.Header.Size;
                if x <= radius*2 | y <= radius*2 | x > imagesize(1)-radius*2-1 | y > imagesize(2)-radius*2-1

                    lowx = x-2*radius;
                    lowy = y-2*radius;
                    highx = x+2*radius-1;
                    highy = y+2*radius-1;

                    %set x or y to 1 if particle is in the left or upper edge
                    if x <= radius*2
                        lowx = 1;
                    end
                    if y <= radius*2
                        lowy = 1;
                    end

                    %set x or y to size of image if particle is in the right or lower edge
                    if x > imagesize(1)-radius*2-1
                        highx = imagesize(1);
                    end
                    if y > imagesize(2)-radius*2-1
                        highy = imagesize(2);
                    end

                    %cut out particle, this will give a non quadratic matrix
                    part_box = tom_emreadc([storage_av2_particlepicker.align(1,j).filename],'subregion',[lowx lowy 1],[highx-lowx highy-lowy 0]);
                    part_box = single(part_box.Value);

                    %taper in x direction
                    if size(part_box,1) < radius*4
                        if lowx == 1
                            stripe = part_box(1,:);
                            while size(part_box,1) < radius*4
                                part_box = cat(1,stripe,part_box);
                            end
                        else
                            stripe = part_box(size(part_box,1),:);
                            while size(part_box,1) < radius*4
                                part_box = cat(1,part_box,stripe);
                            end
                        end
                    end
                    %taper in y direction
                    if size(part_box,2) < radius*4
                        if lowy == 1
                            stripe = part_box(:,1);
                            while size(part_box,2) < radius*4
                                part_box = cat(2,stripe,part_box);
                            end
                        else
                            stripe = part_box(:,size(part_box,2));
                            while size(part_box,2) < radius*4
                                part_box = cat(2,part_box,stripe);
                            end
                        end
                    end

                    particle = part_box;
                else
                    particle = tom_emreadc(storage_av2_particlepicker.align(1,j).filename,'subregion',[x-2*radius y-2*radius 1],[4*radius-1 4*radius-1 0]);
                    particle = single(particle.Value);
                end

                [rot trans ccc moved_part] = tom_av2_align(ref_image,tom_norm(particle,'phase'),mask,'',mask_cc_trans,'',5,0);
                
                %add phase normed aligned particle to particle reference image of current class
                storage_av2_particlepicker.reference_image.(classname) = storage_av2_particlepicker.reference_image.(classname) + moved_part;
                
                %update alignment structure with aligned values of this particle
                manage_align('alignparticle',j,trans(1),trans(2),rot(1),ccc);

            end

            waitbar((j+((run-1)*numpoints))./(numpoints*runs),h,['Refining class ' classname ', run ' num2str(run) ' of ' num2str(runs)]);
            j = j + 1;
        end
    end

    change(run) = 0;
    for k=1:size(ref_before,1)
        for l=1:size(ref_before,2)
            change = change + (ref_before(k,l) - storage_av2_particlepicker.reference_image.(classname)(k,l)).^2;
        end
    end
    ref_before = storage_av2_particlepicker.reference_image.(classname);
    
end

if runs > 3
    figure;plot(change);
end
render_avg();
render_filtered_avg();

chstring = '';
for i=1:size(change,2)
    chstring = [chstring sprintf('%0.2g',change(i)) ', '];
end

set(findobj('Tag','refinechange'),'String',['change: ' chstring(1:end-2)]);

close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Load Reference                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_loadref_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;
[filename, pathname] = uigetfile('*.em', 'Pick reference');
if ischar(filename)
    contents = get(findobj('Tag','menu_classes'),'String');
    val = get(findobj('Tag','menu_classes'),'Value');
    classname = contents{val};
    im = tom_emreadc([pathname '/' filename]);
    storage_av2_particlepicker.reference_image.(classname) = im.Value;
    render_avg();
    render_filtered_avg();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Display particle numbers                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_particlenumbers_Callback(hObject, eventdata, handles)

render_points();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Perform alignment                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alignment_switch_Callback(hObject, eventdata, handles)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average filter enable                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function avg_filter_switch_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    low = str2num(get(findobj('Tag','avg_filter_low'),'String'));
    high = str2num(get(findobj('Tag','avg_filter_high'),'String'));

    if isempty(low) | ~isnumeric(low)
        errordlg('Low filter value must be a positive integer.');
        set(hObject,'Value',0);
        return;
    end

    if isempty(high) | ~isnumeric(high)
        errordlg('High filter value must be a positive integer.');
        set(hObject,'Value',0);
        return;
    end

    if low >= high
        errordlg('High filter value must be bigger than low filter value.');
        set(hObject,'Value',0);
        return;
    end

    render_filtered_avg();
else
    tmpobj = findobj('Tag','particle_average_filtered');
    axes(tmpobj);cla;
    set(tmpobj,'Tag','particle_average_filtered');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average filter low                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function avg_filter_low_Callback(hObject, eventdata, handles)

if get(findobj('Tag','avg_filter_switch'),'Value') == 1
    low = str2num(get(findobj('Tag','avg_filter_low'),'String'));
    high = str2num(get(findobj('Tag','avg_filter_high'),'String'));

    if isempty(low) | ~isnumeric(low)
        errordlg('Low filter value must be a positive integer.');
        set(findobj('Tag','avg_filter_switch'),0);
        return;
    end

    if isempty(high) | ~isnumeric(high)
        errordlg('High filter value must be a positive integer.');
        set(findobj('Tag','avg_filter_switch'),0);
        return;
    end

    if low >= high
        errordlg('High filter value must be bigger than low filter value.');
        set(findobj('Tag','avg_filter_switch'),0);
        return;
    end
    
    render_filtered_avg();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average filter high                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function avg_filter_high_Callback(hObject, eventdata, handles)

if get(findobj('Tag','avg_filter_switch'),'Value') == 1
    low = str2num(get(findobj('Tag','avg_filter_low'),'String'));
    high = str2num(get(findobj('Tag','avg_filter_high'),'String'));

    if isempty(low) | ~isnumeric(low)
        errordlg('Low filter value must be a positive integer.');
        set(findobj('Tag','avg_filter_switch'),0);
        return;
    end

    if isempty(high) | ~isnumeric(high)
        errordlg('High filter value must be a positive integer.');
        set(findobj('Tag','avg_filter_switch'),0);
        return;
    end

    if low >= high
        errordlg('High filter value must be bigger than low filter value.');
        set(findobj('Tag','avg_filter_switch'),0);
        return;
    end
    
    render_filtered_avg();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average rotate                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function avg_rotate_apply_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

contents = get(findobj('Tag','menu_classes'),'String');
classname = contents{get(findobj('Tag','menu_classes'),'Value')};

manage_align('refineallparticles',classname,0,0,str2num(get(findobj('Tag','avg_rotate'),'String')));

storage_av2_particlepicker.reference_image.(classname) = tom_rotate(single(storage_av2_particlepicker.reference_image.(classname)),str2num(get(findobj('Tag','avg_rotate'),'String')));

render_avg();
render_filtered_avg();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average center                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function avg_center_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

%pick new center
axes(findobj('Tag','particle_average'));
[x,y] = ginput(1);

contents = get(findobj('Tag','menu_classes'),'String');
classname = contents{get(findobj('Tag','menu_classes'),'Value')};

%Calculate shift vector
shiftx = -(x-(floor(size(storage_av2_particlepicker.reference_image.(classname),1)./4)+1));
shifty = -(y-(floor(size(storage_av2_particlepicker.reference_image.(classname),2)./4)+1));

%Update alignment structure
manage_align('refineallparticles',classname,shiftx,shifty,0);

%update reference image
storage_av2_particlepicker.reference_image.(classname) = single(tom_shift(storage_av2_particlepicker.reference_image.(classname), [shiftx, shifty]));

%display new images
render_avg();
render_filtered_avg();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Powerspectrum average button                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function powerspectrum_av_Callback(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display Powerspectrum                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_powerspectrum_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    render_ps();
    set(findobj('Tag','histogram_power'),'Enable','on');
else
    tmpobj = findobj('Tag','powerspectrum');
    axes(tmpobj);cla;
    set(tmpobj,'Tag','powerspectrum');
    set(findobj('Tag','histogram_power'),'Enable','off');
    set(findobj('Tag','histogram_power'),'Value',0);
    set(findobj('Tag','histogram_emimage'),'Value',1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  average Powerspectrum                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function powerspectrum_averageone_Callback(hObject, eventdata, handles)

if get(findobj('Tag','display_powerspectrum'),'Value') == 1
    render_ps();
    set(findobj('Tag','histogram_power'),'Enable','on');
else
    tmpobj = findobj('Tag','powerspectrum');
    axes(tmpobj);cla;
    set(tmpobj,'Tag','powerspectrum');
    set(findobj('Tag','histogram_power'),'Enable','off');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Powerspectrum Central region                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_pscenter_Callback(hObject, eventdata, handles)

render_ps();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Powerspectrum in external window                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function powerspectrum_window_Callback(hObject, eventdata, handles)

render_ps(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Selection Tool good image                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_selection_good_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

if storage_av2_particlepicker.imagenumber < size(storage_av2_particlepicker.dircell,2)
    render_image(storage_av2_particlepicker.imagenumber+1);
else
    msgbox('End of image stack reached.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Selection Tool bad image                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_selection_bad_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

storage_av2_particlepicker.flagcell{storage_av2_particlepicker.imagenumber} = 0;

%update image counters
set(findobj('Tag','images_bad'),'String', num2str(str2num(get(findobj('Tag','images_bad'),'String')) + 1));
set(findobj('Tag','images_good'),'String', num2str(str2num(get(findobj('Tag','images_good'),'String')) - 1));

if storage_av2_particlepicker.imagenumber < size(storage_av2_particlepicker.dircell,2)
    render_image(storage_av2_particlepicker.imagenumber+1);
else
    msgbox('End of image stack reached.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Selection Tool bad image                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keypress_Callback(hObject, eventdata, handles)

str=get(gcf,'CurrentCharacter');

% d means good
% k means bad

if strcmp(str,'d')
    button_selection_good_Callback(hObject, eventdata, handles)
elseif strcmp(str,'k')
    button_selection_bad_Callback(hObject, eventdata, handles)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Selection Tool save                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_selectionsave_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

%get filename for stack
[filename, pathname] = uiputfile('*.mat', 'Save Filter as');
if ischar(filename)
   stackname = [pathname '/' filename]; 
   particlepicker.filefilter = storage_av2_particlepicker.flagcell;
   particlepicker.filelist = storage_av2_particlepicker.dircell;
   save([pathname '/' filename],'particlepicker');
   disp('Filter list saved.');
else
    disp('Aborted');
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Selection Tool load                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_selectionload_Callback(hObject, eventdata, handles)

global storage_av2_particlepicker;

[filename, pathname] = uigetfile('*.mat', 'Load Filter');
if ischar(filename)
    s = load([pathname '/' filename]);

    if isfield(s,'particlepicker') ~= 1
        errordlg('This is not a filter file.','File Error');
    end
    
    storage_av2_particlepicker.flagcell = s.particlepicker.filefilter;
    storage_av2_particlepicker.dircell = s.particlepicker.filelist;
    render_image(storage_av2_particlepicker.imagenumber);
    set(findobj('Tag','images_bad'),'String', num2str(size(find(cell2mat(storage_av2_particlepicker.flagcell)==0),2)));
    set(findobj('Tag','images_good'),'String', num2str(size(find(cell2mat(storage_av2_particlepicker.flagcell)==1),2)));
    disp('File filter loaded sucessfully.');
end


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
%%  get all the em-files in a directory                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dircell, flagcell_out] = get_dircontents(directory, dircell, flagcell)

%Get list of EM-Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Getting directory list...');

dirlist = dir(directory);
files = size(dirlist,1);
j = size(dircell,2);
if j == 0
    j = 1;
end

%sort out all the em-files
for i = 1:files
    if dirlist(i).isdir == 0
        if isempty(strmatch(dirlist(i).name,dircell))
            if tom_isemfile([directory '/' dirlist(i).name]) == 1
                dircell{j} = dirlist(i).name;
                try
                    flagcell_out{j} = flagcell{j};
                catch
                    flagcell_out{j} = 1;
                end
                j = size(dircell,2) + 1;
            end
        end
    end
    waitbar(i./files,h,[num2str(i), ' of ', num2str(files), ' files scanned.']);
end

close(h);

if size(dircell,2) == 0
    errordlg('No EM files could be found in this directory.');
    return;
end

set(findobj('Tag','fileslider'),'Max',size(dircell,2),'SliderStep',[1./size(dircell,2) 1./size(dircell,2)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display image                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_image(number, forcereload)

if nargin < 2
    forcereload = 0;
end

global storage_av2_particlepicker;

%skip images marked as bad
if storage_av2_particlepicker.flagcell{number} == 0
    if number >= storage_av2_particlepicker.imagenumber
        if number > size(storage_av2_particlepicker.flagcell)
            msgbox('End of image stack reached.');
        else   
            render_image(number+1);
        end
    else
        if number == 1
            msgbox('Beginning of image stack reached.');
        else
            render_image(number-1);
        end
    end
else

    %read in new image if necessary, otherwise reuse current one that is is
    %memory
    if storage_av2_particlepicker.imagenumber ~= number | forcereload == 1
       im = tom_emreadc([storage_av2_particlepicker.path '/' storage_av2_particlepicker.dircell{number}],'resample',[storage_av2_particlepicker.resampleval storage_av2_particlepicker.resampleval 1]);
       im.Value = single(im.Value);
       storage_av2_particlepicker.image = im.Value;
       storage_av2_particlepicker.header = im.Header;
       [storage_av2_particlepicker.DataScale(1) storage_av2_particlepicker.DataScale(2)] = calculate_histogram(storage_av2_particlepicker.image);
    end

    if (get(findobj('Tag','bandpass_switch'),'Value') == 1)
        im = tom_bandpass(storage_av2_particlepicker.image,str2num(get(findobj('Tag','bandpass_low'),'String')), str2num(get(findobj('Tag','bandpass_high'),'String')));
        DataScale = single(storage_av2_particlepicker.DataScale) ./ (size(im,1) .* size(im,2));
        if str2num(get(findobj('Tag','bandpass_low'),'String')) > 0
            DataScale = DataScale - mean(mean(storage_av2_particlepicker.image)) ./ (size(im,1) .* size(im,2));
        end
    else
        im = storage_av2_particlepicker.image;
        DataScale = storage_av2_particlepicker.DataScale;
    end

    storage_av2_particlepicker.imagenumber = number;
    tmpobj = findobj('Tag','image');
    axes(tmpobj);
    imhandle = imagesc(im',DataScale);
    set(imhandle, 'Tag', 'main_image');
    axis ij; axis off; colormap gray;
    set(tmpobj,'Tag','image');
    set(findobj('Tag','filename_top'),'String',[storage_av2_particlepicker.path '/' storage_av2_particlepicker.dircell{number}]);
    set(findobj('Tag','filename'),'String',storage_av2_particlepicker.dircell{number});
    set(findobj('Tag','fileslider'),'Value',number);

    %render existing points for this image
    render_points();

    if get(findobj('Tag','button_particle_pick'),'Value') == 1
        set(findobj('Tag','main_image'),'buttonDownFcn',@pick_particle);
    end

    %render powerspectrum
    if get(findobj('Tag','display_powerspectrum'),'Value') == 1
        render_ps();
    end

    display_imageinfo();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display powerspectrum                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_ps(externflag,range)

global storage_av2_particlepicker;

if nargin < 1
    externflag = 0;
end

im = storage_av2_particlepicker.image;

%split image
split = str2num(get(findobj('Tag','powerspectrum_averageone'),'String'));

if split ~= 1
    ps_out = zeros(size(im)./split);
    for i=1:split
        if (i==1) 
            im_old=im;
        end
        im = split_image(im_old,split,i);
        im = tom_smooth(im,32);
        ps = log(tom_ps(im));
        ps_out = ps_out + ps;    
    end
else
    ps_out = log(tom_ps(im));
end

storage_av2_particlepicker.powerspectrum = ps_out;

%cut out central part of power spectrum
size_ps = size(ps_out);
center = floor(size_ps ./ 2 + 1);
width = floor(center ./ 2 - 1);

%circular integration over ps
integ_ps = tom_cart2polar(ps_out);
integ_ps = sum(integ_ps,2);
tmpobj = findobj('Tag','ctffit');
axes(tmpobj);
cla;

%plot theoretical ctf
Dz = storage_av2_particlepicker.header.Defocus./10000;
pix_size = storage_av2_particlepicker.header.Objectpixelsize./10;
voltage = storage_av2_particlepicker.header.Voltage./1000;
pixs = size(ps_out,1);
Cs = storage_av2_particlepicker.header.Cs;
alpha = 0.02;
Cc = 2.2;
deltaE = 0.8;
ctf_out = tom_ctf(Dz, pix_size, voltage, pixs, Cs, alpha, Cc, deltaE, 0);
hold on;
plot(tom_norm(integ_ps,3),'g','LineWidth',1.5);
hold off;
set(tmpobj,'Tag','ctffit','XLim',[1 width(1)]);

if get(findobj('Tag','checkbox_pscenter'),'Value') == 1
    ps_out = ps_out(center(1)-width(1):center(1)+width(1)-1,center(2)-width(2):center(2)+width(2)-1);
end

%display power spectrum
if externflag == 0
    tmpobj = findobj('Tag','powerspectrum');
    axes(tmpobj);
    if nargin < 2
        try;tom_imagesc(ps_out','noinfo');end;
    else
        try;imagesc(ps_out',storage_av2_particlepicker.powerScale);end;
    end

    axis ij; axis off; colormap gray;
    set(tmpobj,'Tag','powerspectrum');
else
    if isempty(findobj('Tag','externps'))
        figure;axes;
    else
        axes(findobj('Tag','externps'));
    end
    if nargin < 2
        tom_imagesc(ps_out');
    else
        try;imagesc(ps_out',storage_av2_particlepicker.powerScale);end;
    end
    
    set(gca,'Tag','externps');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate histogram                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [min,max] = calculate_histogram(image)

[mean max min std] = tom_dev(image,'noinfo');
min = mean-2.*std;
max = mean+2.*std;

[h,n] = tom_hist3d(image);
h = 200 .* h ./ (100.*size(image,1) .* size(image,2));
axesobj = findobj('Tag','histogram');
axes(axesobj); 
bar(n,h);
axis auto;
set(axesobj,'Tag','histogram');
set(findobj('Tag','histogram_low'),'String',num2str(round(min)));
set(findobj('Tag','histogram_high'),'String',num2str(round(max)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  split image for ps averaging                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function split_image=split_image(im,number_of_splits,split_nr)
im_sz=size(im,1);

inkre=round(im_sz./number_of_splits);
start=((split_nr-1).*inkre)+1;
stop=((split_nr).*inkre);
split_image=im(start:stop,start:stop);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  update classes                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_classes(varargin)

global storage_av2_particlepicker;

if strcmp(varargin{1},'render')
    
    set(findobj('Tag','menu_classes'),'String',storage_av2_particlepicker.classes);
    
elseif strcmp(varargin{1},'add')
    
    storage_av2_particlepicker.classes{size(storage_av2_particlepicker.classes,2) + 1} = varargin{2};
    storage_av2_particlepicker.classcolors{size(storage_av2_particlepicker.classcolors,2) + 1} = varargin{3};
    set(findobj('Tag','menu_classes'),'String',storage_av2_particlepicker.classes);
    
    classes_stats();
    
elseif strcmp(varargin{1},'delete')
    
    %delete class label and color from cells
    entry = strmatch(varargin{2},storage_av2_particlepicker.classes);
    j = 1;
    for i=1:size(storage_av2_particlepicker.classes,2)
        if i ~= entry
            tmpcell{j} = storage_av2_particlepicker.classes{i};
            tmpcell2{j} = storage_av2_particlepicker.classcolors{i};
            try
                tmpcell3{j} = storage_av2_particlepicker.classradii{i};
            end
            j = j + 1;
        end
    end
    
    storage_av2_particlepicker.classes = tmpcell;
    storage_av2_particlepicker.classcolors = tmpcell2;
    try
        storage_av2_particlepicker.classradii = tmpcell3;
    end
    
    if get(findobj('Tag','menu_classes'),'Value') > size(storage_av2_particlepicker.classes,2)
        set(findobj('Tag','menu_classes'),'Value', size(storage_av2_particlepicker.classes,2));
    end
    set(findobj('Tag','menu_classes'),'String',storage_av2_particlepicker.classes);
    manage_align('deletepointsinclass',varargin{2});
    
    %update radius input box
    contents = get(findobj('Tag','menu_classes'),'String');
    classname = contents{get(findobj('Tag','menu_classes'),'Value')};
    set(findobj('Tag','particle_radius'),'String',num2str(manage_align('getclassradius',classname)));
    
    classes_stats();
    
else
    errormsg('Undefined class action');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  pick particles button down function                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pick_particle(a,b)

global storage_av2_particlepicker;

point1=get(gca,'currentpoint');
button = get(gcf,'selectiontype');

%button values:
%normal: left mouse button
%alt: right mouse button
%extend: middle mouse buttons

pt = point1(1,1:2);
x1 = round(pt(1));
y1 = round(pt(2));

%Handle pick particle event (left mouse button)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(button,'normal') == true
    contents = get(findobj('Tag','menu_classes'),'String');
    classname = contents{get(findobj('Tag','menu_classes'),'Value')};
    filename = get(findobj('Tag','filename_top'),'String'); 
    dataset = '';

    pointlabel = manage_align('addparticle',x1*storage_av2_particlepicker.resampleval,y1*storage_av2_particlepicker.resampleval,classname,filename,dataset);
    numpoints = manage_align('countfilemembers',get(findobj('Tag','filename_top'),'String'));
    if get(findobj('Tag','display_particlenumbers'),'Value') == 1
        label = num2str(numpoints);
    else
        label = '';
    end
    drawmark(x1,y1,manage_align('getclasscolor',classname),strcat('point_',pointlabel),label);
    
    %cut out the particle and display
    x1 = x1*storage_av2_particlepicker.resampleval;
    y1 = y1*storage_av2_particlepicker.resampleval;
    radius = round(str2num(get(findobj('Tag','particle_radius'),'String')));
    
    %test if particle is too close to the edge
     imagesize = size(storage_av2_particlepicker.image)*storage_av2_particlepicker.resampleval;
%     if x1 <= radius | y1 <= radius | x1 > imagesize(1)-radius | y1 > imagesize(2)-radius
%         errordlg('Particle to too close to the border.','Error');
%         manage_align('deleteparticle',manage_align('countallparticles'));
%         render_points();
%         return;
%     end

    %taper mode
    centerpos = radius.*2;
    if x1 <= radius*2 | y1 <= radius*2 | x1 > imagesize(1)-radius*2-1 | y1 > imagesize(2)-radius*2-1
        
        lowx = x1-2*radius;
        lowy = y1-2*radius;
        highx = x1+2*radius-1;
        highy = y1+2*radius-1;
        
        %set x or y to 1 if particle is in the left or upper edge
        if x1 <= radius*2
            lowx = 1;
        end
        if y1 <= radius*2
            lowy = 1;
        end
        
        %set x or y to size of image if particle is in the right or lower edge        
        if x1 > imagesize(1)-radius*2-1
            highx = size(storage_av2_particlepicker.image,1)*storage_av2_particlepicker.resampleval;
        end
        if y1 > imagesize(2)-radius*2-1
            highy = size(storage_av2_particlepicker.image,2)*storage_av2_particlepicker.resampleval;
        end
        
        %cut out particle, this will give a non quadratic matrix
        part_box = tom_emreadc([storage_av2_particlepicker.path '/' storage_av2_particlepicker.dircell{storage_av2_particlepicker.imagenumber}],'subregion',[lowx lowy 1],[highx-lowx highy-lowy 0]);
        part_box = single(part_box.Value);
        
        %taper in x direction
        if size(part_box,1) < radius*4
            if lowx == 1
                stripe = part_box(1,:);
                while size(part_box,1) < radius*4
                    part_box = cat(1,stripe,part_box);
                end
            else
                stripe = part_box(size(part_box,1),:);
                while size(part_box,1) < radius*4
                    part_box = cat(1,part_box,stripe);
                end
            end
        end
        %taper in y direction
        if size(part_box,2) < radius*4
            if lowy == 1
                stripe = part_box(:,1);
                while size(part_box,2) < radius*4
                    part_box = cat(2,stripe,part_box);
                end
            else
                stripe = part_box(:,size(part_box,2));
                while size(part_box,2) < radius*4
                    part_box = cat(2,part_box,stripe);
                end
            end
        end
        
        
    else
        %normal cutout mode
        part_box = tom_emreadc([storage_av2_particlepicker.path '/' storage_av2_particlepicker.dircell{storage_av2_particlepicker.imagenumber}],'subregion',[x1-2*radius y1-2*radius 1],[4*radius-1 4*radius-1 0]);
        part_box = single(part_box.Value);
    end
    
    clickedparticle=part_box(round(centerpos-radius+1):round(centerpos+radius),round(centerpos-radius+1):round(centerpos+radius));

    axesobj = findobj('Tag','particle_clicked');
    axes(axesobj);
    imagesc(clickedparticle',storage_av2_particlepicker.DataScale);
    axis ij; axis off; colormap gray;
    set(axesobj,'Tag','particle_clicked');
    
    %Align particle
    if get(findobj('Tag','alignment_switch'),'Value') == 1

        if isfield(storage_av2_particlepicker.reference_image,classname) == 0
            storage_av2_particlepicker.reference_image.(classname) = tom_norm(part_box,'phase');
            %storage_av2_particlepicker.reference_image.(classname) = part_box;
        end
        ref_image = storage_av2_particlepicker.reference_image.(classname);
        
        mask = tom_spheremask(ones(size(ref_image)),size(ref_image,1)./2-size(ref_image,1)/5,size(ref_image,1)/10);
        mask_cc_trans=tom_spheremask(ones(size(ref_image)),round(str2num(get(findobj('Tag','alignment_translation'),'String'))),1);
        [rot trans ccc moved_part] = tom_av2_align(ref_image,tom_norm(part_box,'phase'),mask,'',mask_cc_trans,'',5,0);
        
        %moved_part = tom_norm(moved_part,'phase');
        
        %display aligned particle    
        axesobj = findobj('Tag','particle_aligned');
        axes(axesobj);
        [h,n]=tom_hist3d(moved_part(round(centerpos-radius+1):round(centerpos+radius),round(centerpos-radius+1):round(centerpos+radius)));
        Scale=[n(1)  n(size(n,2))]; 
        imagesc(moved_part(round(centerpos-radius+1):round(centerpos+radius),round(centerpos-radius+1):round(centerpos+radius))',Scale);
        axis ij; axis off; colormap gray;
        set(axesobj,'Tag','particle_aligned');
        
        %display average particle
        render_avg();
        
        %add phase normed aligned particle to particle reference image of current class
        storage_av2_particlepicker.reference_image.(classname) = storage_av2_particlepicker.reference_image.(classname) + moved_part;

        %display new particle class average image
        render_filtered_avg();
        
        %update alignment structure with aligned values of this particle
        manage_align('alignparticle',numpoints,trans(1),trans(2),rot(1),ccc);
        
    end
    
    check_radiusbox();
    classes_stats();
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  check if radius input field should be disabled or enabled          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_radiusbox()
    
contents = get(findobj('Tag','menu_classes'),'String');
classname = contents{get(findobj('Tag','menu_classes'),'Value')};
if manage_align('countclassmembers',classname) > 0
    %deactivate radius input field for classes with defined radius
    set(findobj('Tag','particle_radius'),'Enable','off'); 
else
    %enable for classes with no members
    set(findobj('Tag','particle_radius'),'Enable','on'); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display filtered average image                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_avg()

global storage_av2_particlepicker;

contents = get(findobj('Tag','menu_classes'),'String');
classname = contents{get(findobj('Tag','menu_classes'),'Value')};

radius = round(str2num(get(findobj('Tag','particle_radius'),'String')));
centerpos = radius.*2;
im = storage_av2_particlepicker.reference_image.(classname)(round(centerpos-radius+1):round(centerpos+radius),round(centerpos-radius+1):round(centerpos+radius));

axesobj = findobj('Tag','particle_average');
axes(axesobj);
[h,n]=tom_hist3d(im);
Scale=[n(1)  n(size(n,2))]; 
imagesc(im',Scale);
axis ij; axis off; colormap gray;
set(axesobj,'Tag','particle_average');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display filtered average image                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_filtered_avg()

global storage_av2_particlepicker;

contents = get(findobj('Tag','menu_classes'),'String');
classname = contents{get(findobj('Tag','menu_classes'),'Value')};

if get(findobj('Tag','avg_filter_switch'),'Value') == 1 & isfield(storage_av2_particlepicker.reference_image,classname)

    radius = round(str2num(get(findobj('Tag','particle_radius'),'String')));
    centerpos = radius.*2;
    im = storage_av2_particlepicker.reference_image.(classname)(round(centerpos-radius+1):round(centerpos+radius),round(centerpos-radius+1):round(centerpos+radius));
    
    %display new particle class average filtered image
    axesobj = findobj('Tag','particle_average_filtered');
    axes(axesobj);

    low = round(str2num(get(findobj('Tag','avg_filter_low'),'String')));
    high = round(str2num(get(findobj('Tag','avg_filter_high'),'String')));

    [h,n]=tom_hist3d(im);
    Scale=[n(1)  n(size(n,2))]; 
    
    if low > 0
        Scale = Scale - mean(mean(im));
    end
    Scale = Scale ./ size(im,1)^2;

    imagesc(tom_bandpass(im',low,high),Scale);
    axis ij; axis off; colormap gray;
    set(axesobj,'Tag','particle_average_filtered');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  drawmark paints a mark                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = drawmark(x,y,Color,Tag,Label)

global storage_av2_particlepicker;

if nargin == 3 
    Tag = 'a';
    Label = '';
end
if nargin == 4 
    Label = '';
end

if ~strcmp(get(gca,'Tag'),'image')
    axes(findobj('Tag','image'));
end

hold on;
Center= x + y*sqrt(-1);
%Radius = 5;
Radius = 30./storage_av2_particlepicker.resampleval;
%Gridpt = 100;
%[u,v]=circle(Center,Radius,Gridpt);
%line(u,v,'LineWidth',1,'Color',[1 0 0]);
uu = [x x x x-Radius x+Radius];
vv = [y-Radius y+Radius y y y];
h = line(uu,vv,'LineWidth',2,'color',Color,'Tag',Tag);
if ~isempty(Label)
    t = text(x+Radius,y-Radius,Label,'Color',Color,'Tag',Tag,'FontSize',14,'FontWeight','bold');
end
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  align structure management                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = manage_align(varargin)

global storage_av2_particlepicker;

if strcmp(varargin{1},'create')
    storage_av2_particlepicker.align = struct();
    storage_av2_particlepicker.points = struct();
    out = '';
    
%addparticle: x,y,class,filename,dataset
elseif strcmp(varargin{1},'addparticle')
    pointnumber = size(storage_av2_particlepicker.align,2)+1;
    if pointnumber ~= 1
        if ~isfield(storage_av2_particlepicker.align(1,pointnumber-1),'filename')
            pointnumber = pointnumber-1;
        end
    end
    storage_av2_particlepicker.align(1,pointnumber).dataset = varargin{6};
    storage_av2_particlepicker.align(1,pointnumber).filename = varargin{5};
    storage_av2_particlepicker.align(1,pointnumber).position.x = varargin{2};
    storage_av2_particlepicker.align(1,pointnumber).position.y = varargin{3};
    storage_av2_particlepicker.align(1,pointnumber).class = varargin{4};
    storage_av2_particlepicker.align(1,pointnumber).label = tempname;
    storage_av2_particlepicker.align(1,pointnumber).radius = round(str2num(get(findobj('Tag','particle_radius'),'String')));
    storage_av2_particlepicker.align(1,pointnumber).color = storage_av2_particlepicker.classcolors{get(findobj('Tag','menu_classes'),'Value')};
    storage_av2_particlepicker.align(1,pointnumber).shift.x = 0;
    storage_av2_particlepicker.align(1,pointnumber).shift.y = 0;
    storage_av2_particlepicker.align(1,pointnumber).angle = 0;
    storage_av2_particlepicker.align(1,pointnumber).isaligned = 0;
    storage_av2_particlepicker.align(1,pointnumber).ccc = 0;
    storage_av2_particlepicker.align(1,pointnumber).quality = 0;
    storage_av2_particlepicker.align(1,pointnumber).normed = 'none';
    
    [pathstr, name, ext] = fileparts(varargin{5});
    pn = strmatch(strcat(name,ext),storage_av2_particlepicker.dircell);
    
    try
        oldpoints = storage_av2_particlepicker.points(pn).points;
    catch
        oldpoints = [];
    end
    
    storage_av2_particlepicker.points(pn).points = [oldpoints,pointnumber];

    out = storage_av2_particlepicker.align(1,pointnumber).label;

% alignparticle (number,x,y,angle,ccc)
elseif strcmp(varargin{1},'alignparticle')
    pointnumber = varargin{2};
    shiftx = varargin{3};
    shifty = varargin{4};
    angle = varargin{5};
    ccc = varargin{6};
    storage_av2_particlepicker.align(1,pointnumber).shift.x = shiftx;
    storage_av2_particlepicker.align(1,pointnumber).shift.y = shifty;
    storage_av2_particlepicker.align(1,pointnumber).angle = angle;
    storage_av2_particlepicker.align(1,pointnumber).isaligned = 1;
    storage_av2_particlepicker.align(1,pointnumber).ccc = ccc;

    % refineparticle (number,x,y,angle,ccc)
elseif strcmp(varargin{1},'refineparticle')
    pointnumber = varargin{2};
    shiftx = varargin{3};
    shifty = varargin{4};
    angle = varargin{5};
    ccc = varargin{6};
    if storage_av2_particlepicker.align(1,pointnumber).isaligned == 1
        storage_av2_particlepicker.align(1,pointnumber).shift.x = storage_av2_particlepicker.align(1,pointnumber).shift.x + shiftx;
        storage_av2_particlepicker.align(1,pointnumber).shift.y = storage_av2_particlepicker.align(1,pointnumber).shift.x + shifty;
        storage_av2_particlepicker.align(1,pointnumber).angle = storage_av2_particlepicker.align(1,pointnumber).angle + angle;
        storage_av2_particlepicker.align(1,pointnumber).ccc = ccc;
    end
    
% refineallparticles (classname,x,y,angle)
elseif strcmp(varargin{1},'refineallparticles')
    classname = varargin{2};
    shiftx = varargin{3};
    shifty = varargin{4};
    angle = varargin{5};
    for pointnumber=1:size(storage_av2_particlepicker.align,2)
        %apply refinement to selected class and only prealigned particles
        if storage_av2_particlepicker.align(1,pointnumber).isaligned == 1 & strcmp(storage_av2_particlepicker.align(1,pointnumber).class,classname)
            storage_av2_particlepicker.align(1,pointnumber).shift.x = storage_av2_particlepicker.align(1,pointnumber).shift.x + shiftx;
            storage_av2_particlepicker.align(1,pointnumber).shift.y = storage_av2_particlepicker.align(1,pointnumber).shift.x + shifty;
            storage_av2_particlepicker.align(1,pointnumber).angle = storage_av2_particlepicker.align(1,pointnumber).angle + angle;
        end
    end
    
elseif strcmp(varargin{1},'deleteparticle')
    storage_av2_particlepicker.align(varargin{2}) = []; 
    manage_align('rebuildindex');
    render_points();
    
elseif strcmp(varargin{1},'deletepointsinclass')
    classname = varargin{2};
    if manage_align('countclassmembers',classname) > 0
        tmpstruct = storage_av2_particlepicker.align;
        for i=size(storage_av2_particlepicker.align,2):-1:1
            if strcmp(storage_av2_particlepicker.align(1,i).class,classname)
                tmpstruct(i) = [];
            end
        end
        storage_av2_particlepicker.align = tmpstruct;
        manage_align('rebuildindex');
        render_points();
    end
    
elseif strcmp(varargin{1},'countclassmembers')
    classname = varargin{2};
    counter = 0;
    for i=1:size(storage_av2_particlepicker.align,2)
        try
            if strcmp(storage_av2_particlepicker.align(1,i).class,classname)
                counter = counter + 1;
            end
        end
    end
    out = counter;

elseif strcmp(varargin{1},'countallparticles')
    out = size(storage_av2_particlepicker.align,2);
    
elseif strcmp(varargin{1},'countfilemembers')
    filename = varargin{2};
    counter = 0;
    for i=1:size(storage_av2_particlepicker.align,2)
        if strcmp(storage_av2_particlepicker.align(1,i).filename,filename)
            counter = counter + 1;
        end
    end
    out = counter;

elseif strcmp(varargin{1},'getpointsincurrentframe')
    pn = strmatch(get(findobj('Tag','filename'),'String'),storage_av2_particlepicker.dircell,'exact');
    out = storage_av2_particlepicker.points(pn).points;

elseif strcmp(varargin{1},'deletepointsincurrentframe')    
    for i = sort(manage_align('getpointsincurrentframe'),'descend');
        storage_av2_particlepicker.align(i) = [];
    end
    manage_align('rebuildindex');
    render_points();

elseif strcmp(varargin{1},'getclasscolor')
    classname = varargin{2};
    %if no particle is found in the alignment list of this class, try to
    %get color from the storage cell
    idx = strmatch(classname,storage_av2_particlepicker.classes,'exact');
    out = storage_av2_particlepicker.classcolors(idx);
    out = out{1};
    
elseif strcmp(varargin{1},'getclassradius')
    classname = varargin{2};
    idx = strmatch(classname,storage_av2_particlepicker.classes,'exact');
    try
        out = storage_av2_particlepicker.classradii(idx);
        out = out{1};
    catch
        out = '';
    end

elseif strcmp(varargin{1},'getclasses')
    newclasses = {};
    newcolors = {};
    newradii = {};
    for i=1:size(storage_av2_particlepicker.align,2)
        if strmatch(storage_av2_particlepicker.align(1,i).class,newclasses,'exact')
        else
            newclasses{size(newclasses,2)+1} = storage_av2_particlepicker.align(1,i).class;
            newcolors{size(newcolors,2)+1} = storage_av2_particlepicker.align(1,i).color;
            newradii{size(newradii,2)+1} = storage_av2_particlepicker.align(1,i).radius;
        end
    end
    out = {newclasses newcolors newradii};
    
elseif strcmp(varargin{1},'rebuildindex')
    newidx = struct();
    for i=1:size(storage_av2_particlepicker.align,2)
        [pathstr, name, ext] = fileparts(storage_av2_particlepicker.align(1,i).filename);
        fileidx = strmatch(strcat(name,ext),storage_av2_particlepicker.dircell);
        try
            oldpoints = newidx(fileidx).points;
        catch
            oldpoints = [];
        end
        newidx(fileidx).points = [oldpoints,i];
    end
    storage_av2_particlepicker.points = newidx;
    
else
    errordlg('Unknown align management action');
    out = '';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  render all points in an image                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_points()

global storage_av2_particlepicker;

counter = 0;

delete(findobj('-regexp','Tag','point_*'));

try
    for i=manage_align('getpointsincurrentframe')
        counter = counter+1;
        if get(findobj('Tag','display_particlenumbers'),'Value') == 1
            label = num2str(counter);
        else
            label = '';
        end
        drawmark(storage_av2_particlepicker.align(1,i).position.x/storage_av2_particlepicker.resampleval, storage_av2_particlepicker.align(1,i).position.y/storage_av2_particlepicker.resampleval,storage_av2_particlepicker.align(1,i).color,strcat('point_',storage_av2_particlepicker.align(1,i).label),label);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Picklist autosave function                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function auto_save()

global storage_av2_particlepicker;

if isfield(storage_av2_particlepicker,'picklistfilename') & ischar(storage_av2_particlepicker.picklistfilename) & get(findobj('Tag','autosave_switch'),'Value') == 1
    align2d = storage_av2_particlepicker.align;
    align2d = rmfield(align2d,'label');
    save(storage_av2_particlepicker.picklistfilename,'align2d');
    disp('Picklist autosaved');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Rebuild reference images                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rebuild_reference_images()

global storage_av2_particlepicker;

%Get list of EM-Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Reconstructing reference images, please be patient...');

%delete previous reference images
storage_av2_particlepicker.reference_image = {};
out = manage_align('getclasses');
for i=1:size(out{1},2)
   storage_av2_particlepicker.reference_image.(out{1}{i}) = zeros(out{3}{i}*4);
end

numframes = size(storage_av2_particlepicker.points,2);
j = 1;
numpoints = size(storage_av2_particlepicker.align,2);

for i = 1:numframes
    filename = storage_av2_particlepicker.dircell{i};

    if size(storage_av2_particlepicker.points(i).points,2) > 1
        image = tom_reademheader([storage_av2_particlepicker.path '/' filename]);
    end
    
    for k = 1:size(storage_av2_particlepicker.points(i).points,2)
        %only add aligned particles to reference image
        if storage_av2_particlepicker.align(1,j).isaligned == 1
            radius = storage_av2_particlepicker.align(1,j).radius;
            x = storage_av2_particlepicker.align(1,j).position.x;
            y = storage_av2_particlepicker.align(1,j).position.y;
            classname = storage_av2_particlepicker.align(1,j).class;
            %cut out bigger particle, then rotate it, then cut it to the
            %final size, then norm

            %taper mode
            imagesize = image.Header.Size;
            if x <= radius*2 | y <= radius*2 | x > imagesize(1)-radius*2-1 | y > imagesize(2)-radius*2-1

                lowx = x-2*radius;
                lowy = y-2*radius;
                highx = x+2*radius-1;
                highy = y+2*radius-1;

                %set x or y to 1 if particle is in the left or upper edge
                if x <= radius*2
                    lowx = 1;
                end
                if y <= radius*2
                    lowy = 1;
                end

                %set x or y to size of image if particle is in the right or lower edge
                if x > imagesize(1)-radius*2-1
                    highx = imagesize(1);
                end
                if y > imagesize(2)-radius*2-1
                    highy = imagesize(2);
                end

                %cut out particle, this will give a non quadratic matrix
                part_box = tom_emreadc(storage_av2_particlepicker.align(1,j).filename,'subregion',[lowx lowy 1],[highx-lowx highy-lowy 0]);
                part_box = single(part_box.Value);

                %taper in x direction
                if size(part_box,1) < radius*4
                    if lowx == 1
                        stripe = part_box(1,:);
                        while size(part_box,1) < radius*4
                            part_box = cat(1,stripe,part_box);
                        end
                    else
                        stripe = part_box(size(part_box,1),:);
                        while size(part_box,1) < radius*4
                            part_box = cat(1,part_box,stripe);
                        end
                    end
                end
                %taper in y direction
                if size(part_box,2) < radius*4
                    if lowy == 1
                        stripe = part_box(:,1);
                        while size(part_box,2) < radius*4
                            part_box = cat(2,stripe,part_box);
                        end
                    else
                        stripe = part_box(:,size(part_box,2));
                        while size(part_box,2) < radius*4
                            part_box = cat(2,part_box,stripe);
                        end
                    end
                end
                
                particle = part_box;
            else
                particle = tom_emreadc(storage_av2_particlepicker.align(1,j).filename,'subregion',[x-2*radius y-2*radius 1],[4*radius-1 4*radius-1 0]);
                particle = single(particle.Value);
            end

            particle = tom_rotate(particle,storage_av2_particlepicker.align(1,j).angle);
            particle = tom_shift(particle, [storage_av2_particlepicker.align(1,j).shift.x storage_av2_particlepicker.align(1,j).shift.y]);
            
            particle = tom_norm(particle,'phase');
            
            %add particle to the correct reference image
            if isfield(storage_av2_particlepicker.reference_image,classname) == 0
                storage_av2_particlepicker.reference_image.(classname) = particle;
            else
                storage_av2_particlepicker.reference_image.(classname) = storage_av2_particlepicker.reference_image.(classname) + particle;
            end
        end
        waitbar(j./numpoints,h);
        j = j + 1;
    end
end

close(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  create particle stack                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_particle_stack(classes, stackname, alignname,alignedflag,normflag,radius)

global storage_av2_particlepicker;

%Get list of EM-Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Creating particle stack, please be patient...');

numframes = size(storage_av2_particlepicker.points,2);
j = 1;
numpoints = size(storage_av2_particlepicker.align,2);

stacksize = 0;
for l=1:size(classes,2)
    stacksize = stacksize + manage_align('countclassmembers',classes{l});
end

stack = zeros(radius*2,radius*2,stacksize);
selectedparticles = [];

for i = 1:numframes
    filename = storage_av2_particlepicker.dircell{i};
    image = tom_reademheader([storage_av2_particlepicker.path '/' filename]);
    
    for k = 1:size(storage_av2_particlepicker.points(i).points,2)

        %add particle to stack if classname matches
        if strmatch(storage_av2_particlepicker.align(1,j).class,classes,'exact');
            %classname = storage_av2_particlepicker.align(1,j).class;
            x = storage_av2_particlepicker.align(1,j).position.x;
            y = storage_av2_particlepicker.align(1,j).position.y;

            %aligned particles
            if alignedflag == 1

                %particle is already aligned
                if storage_av2_particlepicker.align(1,j).isaligned == 0
                    %FIXME: align remaining unaligned particles
                end

            end

            %cut particles
            %taper mode
            imagesize = image.Header.Size;
            if x <= radius*2 | y <= radius*2 | x > imagesize(1)-radius*2-1 | y > imagesize(2)-radius*2-1

                lowx = x-2*radius;
                lowy = y-2*radius;
                highx = x+2*radius-1;
                highy = y+2*radius-1;

                %set x or y to 1 if particle is in the left or upper edge
                if x <= radius*2
                    lowx = 1;
                end
                if y <= radius*2
                    lowy = 1;
                end

                %set x or y to size of image if particle is in the right or lower edge
                if x > imagesize(1)-radius*2-1
                    highx = imagesize(1);
                end
                if y > imagesize(2)-radius*2-1
                    highy = imagesize(2);
                end

                %cut out particle, this will give a non quadratic matrix
                part_box = tom_emreadc([storage_av2_particlepicker.align(1,j).filename],'subregion',[lowx lowy 1],[highx-lowx highy-lowy 0]);
                part_box = single(part_box.Value);

                %taper in x direction
                if size(part_box,1) < radius*4
                    if lowx == 1
                        stripe = part_box(1,:);
                        while size(part_box,1) < radius*4
                            part_box = cat(1,stripe,part_box);
                        end
                    else
                        stripe = part_box(size(part_box,1),:);
                        while size(part_box,1) < radius*4
                            part_box = cat(1,part_box,stripe);
                        end
                    end
                end
                %taper in y direction
                if size(part_box,2) < radius*4
                    if lowy == 1
                        stripe = part_box(:,1);
                        while size(part_box,2) < radius*4
                            part_box = cat(2,stripe,part_box);
                        end
                    else
                        stripe = part_box(:,size(part_box,2));
                        while size(part_box,2) < radius*4
                            part_box = cat(2,part_box,stripe);
                        end
                    end
                end
                
                particle = part_box;
            else
                particle = tom_emreadc(storage_av2_particlepicker.align(1,j).filename,'subregion',[x-2*radius y-2*radius 1],[4*radius-1 4*radius-1 0]);
                particle = single(particle.Value);
            end
            
            %rotate if particle is aligned
            if storage_av2_particlepicker.align(1,j).angle ~= 0 & alignedflag == 1
                particle = tom_rotate(particle,storage_av2_particlepicker.align(1,j).angle);
            end
            %shift particle
            if alignedflag == 1
                particle = tom_shift(particle,[storage_av2_particlepicker.align(1,j).shift.x storage_av2_particlepicker.align(1,j).shift.y]);
            end
            
            %reduce to final size
            particle = particle(radius:3*radius-1,radius:3*radius-1);
            %norm particle
            if normflag == 1
                particle = tom_norm(particle,'phase');
                storage_av2_particlepicker.align(1,j).normed='phase';
            elseif normflag == 2
                particle = tom_norm(particle,'3std');
                storage_av2_particlepicker.align(1,j).normed='3std';
            end
            %apply to stack
            stack(:,:,j) = particle;
            selectedparticles = [selectedparticles, j];
        end
        
        waitbar(j./numpoints,h);
        j = j + 1;
    end
end

tom_emwrite(stackname,stack);
align2d = storage_av2_particlepicker.align(1,selectedparticles);
align2d = rmfield(align2d,'label');
save(alignname,'align2d');
disp('Picklist Saved');

close(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display class statistics                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function classes_stats()

global storage_av2_particlepicker;

string = '';

classnames = storage_av2_particlepicker.classes;
classcolors = storage_av2_particlepicker.classcolors;

for i = 1:size(classnames,2)
    numparticles = manage_align('countclassmembers',classnames{i});
    string = strvcat(string, ['\color[rgb]{',num2str(classcolors{i}(1)),' ',num2str(classcolors{i}(2)), ' ' ,num2str(classcolors{i}(3)),'}',classnames{i},': ', num2str(numparticles), ' particles']);
end
tmpobj = findobj('Tag','classes_text');axes(tmpobj);axis off;cla;
text(1,1,string,'FontSize',12,'FontWeight','bold','VerticalAlignment','Top');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display image info                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_imageinfo()

global storage_av2_particlepicker;

string = '';
string = strvcat(string,['Image Size: ' num2str(storage_av2_particlepicker.header.Size(1)) ' x ' num2str(storage_av2_particlepicker.header.Size(2))]);
string = strvcat(string,['Pixel Size: ' num2str(storage_av2_particlepicker.header.Objectpixelsize./10) ' nm']);
string = strvcat(string,['Voltage: ' num2str(storage_av2_particlepicker.header.Voltage./1000) ' kV']);
string = strvcat(string,['Defocus: ' num2str(storage_av2_particlepicker.header.Defocus./10000) ' um']);
string = strvcat(string,['Cs: ' num2str(storage_av2_particlepicker.header.Cs) ' mm']);
set(findobj('Tag','image_infotext'),'String',string);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function powerspectrum_avnumber_Callback(hObject, eventdata, handles)
function alignment_iterations_Callback(hObject, eventdata, handles)
function alignment_translation_Callback(hObject, eventdata, handles)
function alignment_translation_y_Callback(hObject, eventdata, handles)
function alignment_angle_Callback(hObject, eventdata, handles)
function refinenumber_Callback(hObject, eventdata, handles)
function histogram_low_Callback(hObject, eventdata, handles)
function histogram_high_Callback(hObject, eventdata, handles)
function avg_rotate_Callback(hObject, eventdata, handles)
function histogram_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function histogram_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function fileslider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function filename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function bandpass_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function bandpass_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function menu_classes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function particle_radius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function refinenumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function alignment_iterations_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function alignment_translation_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function alignment_translation_y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function alignment_angle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function avg_rotate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function avg_shift_x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function avg_shift_y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function avg_filter_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function avg_filter_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function powerspectrum_avnumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function powerspectrum_averageone_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_resample_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dummy(a,b)
