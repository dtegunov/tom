function varargout = tom_av2_autopickergui(varargin)
%TOM_AV2_AUTOPICKERGUI M-file for tom_av2_autopickergui.fig
%
%   tom_amira_createisosurface(filelabel,threshold, label, color, transformmatrix, iconposition, host)
%
%   TOM_AV2_AUTOPICKERGUI, by itself, creates a new TOM_AV2_AUTOPICKERGUI or raises the existing
%   singleton*.
%
%   H = TOM_AV2_AUTOPICKERGUI returns the handle to a new TOM_AV2_AUTOPICKERGUI or the handle to
%    the existing singleton*.
%
%   TOM_AV2_AUTOPICKERGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%    function named CALLBACK in TOM_AV2_AUTOPICKERGUI.M with the given input arguments.
%
%   TOM_AV2_AUTOPICKERGUI('Property','Value',...) creates a new TOM_AV2_AUTOPICKERGUI or raises the
%    existing singleton*.  Starting from the left, property value pairs are
%    applied to the GUI before tom_av2_autopickergui_OpeningFunction gets called.  An
%    unrecognized property name or invalid value makes property application
%    stop.  All inputs are passed to tom_av2_autopickergui_OpeningFcn via varargin.
%
%    *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%     instance to run (singleton)".
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout           ...
%
%EXAMPLE
%   ... = tom_av2_autopickergui(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   GUIDE, GUIDATA, GUIHANDLES
%
% Edit the above text to modify the response to help tom_av2_autopickergui
%
% Last Modified by GUIDE v2.5 20-Dec-2006 09:30:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_av2_autopickergui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_autopickergui_OutputFcn, ...
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
%%  opening function                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tom_av2_autopickergui_OpeningFcn(hObject, eventdata, handles, varargin)


handles.imagecell = {};
handles.parallelstruct = struct();
axes(handles.axes_autopicker);axis off;
handles.imfilter.filter = struct();
% Choose default command line output for tom_av2_autopickergui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  output function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_av2_autopickergui_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  input directory                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_autopicker_dir_Callback(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse input directory                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_autopicker_browsedir_Callback(hObject, eventdata, handles)

pathname = uigetdir('Input directory');
if ischar(pathname)
    set(handles.input_autopicker_dir,'String',pathname);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  input good/bad file                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_autopicker_goodbadfile_Callback(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse good/bad file                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_autopicker_browsegoodbadfile_Callback(hObject, eventdata, handles)

[filename,pathname] = uigetfile('*.mat','Open File selection file');
if ischar(filename)
    set(handles.input_autopicker_goodbadfile,'String',[pathname, '/', filename]);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  template stack                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_autopicker_browsetemplate_Callback(hObject, eventdata, handles)

[filename,pathname] = uigetfile('*.em','Open template stack');
if ischar(filename)
    set(handles.input_autopicker_template,'String',[pathname, '/', filename]);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse fft file                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_autopicker_browsefftfilename_Callback(hObject, eventdata, handles)

[filename,pathname] = uigetfile('*.em','Open fft file');
if ischar(filename)
    set(handles.input_autopicker_fftfilename,'String',[pathname, '/', filename]);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse output align file                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_autopicker_browsealign_Callback(hObject, eventdata, handles)

[filename,pathname] = uiputfile('*.mat','Save alignment to file...');
if ischar(filename)
    set(handles.input_autopicker_align,'String',[pathname, '/', filename]);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  parallel options                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_autopicker_paralleloptions_Callback(hObject, eventdata, handles)


if ~isfield(handles.parallelstruct,'jobmanager')
    handles.parallelstruct = tom_parallelsettings();
else
    handles.parallelstruct = tom_parallelsettings(handles.parallelstruct);
end

infostring = ['host: ' handles.parallelstruct.jobmanager];
infostring = strvcat(infostring,['number of tasks: ' num2str(handles.parallelstruct.number_of_tasks)]);
infostring = strvcat(infostring,['workers: ',num2str(handles.parallelstruct.workers.min) ' - ' num2str(handles.parallelstruct.workers.max)]);
infostring = strvcat(infostring,['timeout: ',num2str(handles.parallelstruct.timeout)]);

set(handles.text_autopicker_parallelinfo,'String',infostring);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load images                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_autopicker_loadimages_Callback(hObject, eventdata, handles)

handles = create_imagecell(handles);

im = tom_reademheader(handles.imagecell{1});

infotext = ['number of images: ' num2str(size(handles.imagecell,2))];
infotext = strvcat(infotext,['image dimensions: ' num2str(im.Header.Size(1)) ' x ' num2str(im.Header.Size(2))]);
set(handles.text_autopicker_infoimages,'String',infotext);

set(handles.slider_autopicker_file,'Min',1,'Max',size(handles.imagecell,2),'Value',1,'SliderStep',[1./size(handles.imagecell,2) 1./size(handles.imagecell,2)]);

handles = render_image(handles,1);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  file input                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_autopicker_imagename_Callback(hObject, eventdata, handles)

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  file slider                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider_autopicker_file_Callback(hObject, eventdata, handles)

position = round(get(hObject,'Value'));
handles = render_image(handles, position);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  test pick                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_autopicker_testpick_Callback(hObject, eventdata, handles)

%get angles
%%%%%%%%%%%%%%%%%%
angle_start = str2double(get(handles.input_autopicker_anglestart,'String'));
angle_stop = str2double(get(handles.input_autopicker_anglestop,'String'));
angle_incr = str2double(get(handles.input_autopicker_angleincrement,'String'));

if isempty(angle_start); errordlg('Missing start angle.');return;end;
if isempty(angle_stop); errordlg('Missing stop angle.');return;end;
if isempty(angle_incr); errordlg('Missing increment angle.');return;end;
if angle_start > angle_stop; errordlg('stop angle > start angle!');return;end;
if angle_incr == 0; errordlg('angle increment is not valid');return;end;

%get options
%%%%%%%%%%%%%%%%%%
binning = round(str2double(get(handles.input_autopicker_binning,'String')));
numberofparticles = round(str2double(get(handles.input_autopicker_noparticles,'String')));

if isempty(binning); errordlg('Missing binning value.');return;end;
if isempty(numberofparticles); errordlg('Missing number of particles per image.');return;end;
if binning < 0; errordlg('binning value must be positive or zero');return;end;
if numberofparticles < 1; errordlg('number of particles per image must be a positive number');return;end;

%load template
%%%%%%%%%%%%%%%%%%
template = tom_emreadc(get(handles.input_autopicker_template,'String'));
handles.templatesize = template.Header.Size;

%get fft filename
%%%%%%%%%%%%%%%%%%

fftfilename = get(handles.input_autopicker_fftfilename,'String');
reusefft = get(handles.checkbox_reusefftfile,'Value');


%start picking
%%%%%%%%%%%%%%%%%%
if isfield(handles,'currentmask')
    [handles.align2d,handles.ccfmap] = tom_av2_autopicker(template.Value, {handles.currentimage}, [], [], angle_start, angle_stop, angle_incr, numberofparticles, binning, fftfilename, reusefft, {handles.currentmask},1,handles.imfilter.filter,handles.parallelstruct);
else
    [handles.align2d,handles.ccfmap] = tom_av2_autopicker(template.Value, {handles.currentimage}, [], [], angle_start, angle_stop, angle_incr, numberofparticles, binning, fftfilename, reusefft, {},1,handles.imfilter.filter,handles.parallelstruct);
end

handles = render_points(handles);
set(handles.checkbox_autopicker_ccfmap,'Value',0);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  show ccfmap                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_autopicker_ccfmap_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    handles = render_ccfmap(handles);
else
    handles = render_image(handles, handles.currimno);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plot points                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_autopicker_peaks_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    handles = render_points(handles);
else
    if get(handles.checkbox_autopicker_ccfmap,'Value') == 1
        handles = render_ccfmap(handles);
    else
        handles = render_image(handles, handles.currimno);
    end

end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  start picking                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_autopicker_start_Callback(hObject, eventdata, handles)

%get angles
%%%%%%%%%%%%%%%%%%
angle_start = str2double(get(handles.input_autopicker_anglestart,'String'));
angle_stop = str2double(get(handles.input_autopicker_anglestop,'String'));
angle_incr = str2double(get(handles.input_autopicker_angleincrement,'String'));

if isempty(angle_start); errordlg('Missing start angle.');return;end;
if isempty(angle_stop); errordlg('Missing stop angle.');return;end;
if isempty(angle_incr); errordlg('Missing increment angle.');return;end;
if angle_start > angle_stop; errordlg('stop angle > start angle!');return;end;
if angle_incr == 0; errordlg('angle increment is not valid');return;end;

%get options
%%%%%%%%%%%%%%%%%%
binning = round(str2double(get(handles.input_autopicker_binning,'String')));
numberofparticles = round(str2double(get(handles.input_autopicker_noparticles,'String')));

if isempty(binning); errordlg('Missing binning value.');return;end;
if isempty(numberofparticles); errordlg('Missing number of particles per image.');return;end;
if binning < 0; errordlg('binning value must be positive or zero');return;end;
if numberofparticles < 1; errordlg('number of particles per image must be a positive number');return;end;

%load template
%%%%%%%%%%%%%%%%%%
template = tom_emreadc(get(handles.input_autopicker_template,'String'));
handles.templatesize = template.Header.Size;

%get fft filename
%%%%%%%%%%%%%%%%%%

fftfilename = get(handles.input_autopicker_fftfilename,'String');
reusefft = get(handles.checkbox_reusefftfile,'Value');

%get output alignment file
%%%%%%%%%%%%%%%%%%
outalign = get(handles.input_autopicker_align,'String');
if isempty(outalign); errordlg('Missing output aligment file.');return;end;

%start picking
%%%%%%%%%%%%%%%%%%
if get(handles.checkbox_autopicker_parallelmode,'Value') == 1
    align2d = tom_av2_autopicker(template.Value, handles.imagecell, outalign, [], angle_start, angle_stop, angle_incr, numberofparticles, binning, fftfilename, reusefft, handles.maskcell, 0, handles.imfilter.filter,handles.parallelstruct);
else
    align2d = tom_av2_autopicker(template.Value, handles.imagecell, outalign, [], angle_start, angle_stop, angle_incr, numberofparticles, binning, fftfilename, reusefft, handles.maskcell, 0, handles.imfilter.filter);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  set image filter                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_autopicker_filterimages_Callback(hObject, eventdata, handles)

in_struct.filter.types = {'bandpass','kernel'};
if ~isempty(fieldnames(handles.imfilter.filter))
    handles.imfilter = tom_filtergui('filter',in_struct,handles.imfilter);
else
    handles.imfilter = tom_filtergui('filter',in_struct);
end

guidata(hObject, handles);


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
%%  display ccfmap                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = render_ccfmap(handles)

axes(handles.axes_autopicker);
imagesc(handles.ccfmap');axis ij;colormap gray;axis off;
impixelinfo();
figure;surfc(double(tom_binc(handles.ccfmap,1))');shading interp;colormap hot;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display points                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = render_points(handles)

axes(handles.axes_autopicker);
color = tom_colorpalette(handles.templatesize(3));
binning = round(str2double(get(handles.input_autopicker_binning,'String')));

if ~isempty(handles.align2d)
    for i=1:size(handles.align2d,2)
        handles = drawmark(handles,handles.align2d(1,i).position.x/2^binning, handles.align2d(1,i).position.y/2^binning,color(str2double(handles.align2d(1,i).class(13:end)),:));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display image                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = render_image(handles, number)

binning = str2double(get(handles.input_autopicker_binning,'String'));
im = tom_emreadc(handles.imagecell{number},'binning',binning);    
%filter image
if ~isempty(fieldnames(handles.imfilter.filter))
    im.Value = tom_apply_filter(im.Value,handles.imfilter.filter);
end
    
[a,b,c,d]=tom_dev(im.Value,'noinfo');
DataScale = [a-2*d a+2*d];
%generate polygon mask
    if size(handles.maskcell,2) >= number && ~isempty(handles.maskcell{number})
        mask_st.Apply=1;
        mask_st.Value.size_x=size(im.Value,1);
        mask_st.Value.size_y=size(im.Value,2);
        mask_st.Polygons = handles.maskcell{number};
        mask_st.Value.binning = binning;
        mask_st.Type='roipoly';
        mask_st.Value.angle_phi = 0;
        mask_st.Value.angle_psi = 0;
        mask_st.Value.angle_theta = 0;
        mask_st.Inverse=1;
        mask=tom_create_mask(mask_st);
        im.Value = im.Value.*mask;
    end

axes(handles.axes_autopicker);
imagesc(im.Value',DataScale);axis ij;colormap gray;axis off;

set(handles.edit_autopicker_imagename,'String',im.Header.Filename);
handles.currentimage = handles.imagecell{number};
try;
    handles.currentmask = handles.maskcell{number};
end
handles.currimno = number;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  create imagecell                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = create_imagecell(handles)

directory = get(handles.input_autopicker_dir,'String');
unix(['chmod -R 777 ' directory]);
nogoodbadfile = 0;
goodbadfile = get(handles.input_autopicker_goodbadfile,'String');

handles.imagecell = {};
handles.maskcell = {};

try;s = load(goodbadfile); if isfield(s,'particlepicker') == 0; nogoodbadfile = 1;end; catch nogoodbadfile = 1;end;

if nogoodbadfile == 0
    j=1;
    files = size(s.particlepicker.filelist,2);
    h = waitbar(0,'Getting directory list...');
    for i=1:files
        if s.particlepicker.filefilter{i} == 1
            if tom_isemfile([directory '/' s.particlepicker.filelist{i}]) == 1
                handles.imagecell{j} = [directory '/' s.particlepicker.filelist{i}];
                if size(s.particlepicker.maskcell,2) >= i
                    handles.maskcell{j} = s.particlepicker.maskcell{i};
                end
                j = j+1;
            end
        end
        waitbar(i./files,h,[num2str(i), ' of ', num2str(files), ' files scanned.']);
    end
    close(h);
else
    [handles.imagecell] = get_dircontents(directory, {},{});
end
    


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
                dircell{j} = [directory '/' dirlist(i).name];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  drawmark paints a mark                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = drawmark(handles,x,y,Color,Tag,Label)

binning = str2double(get(handles.input_autopicker_binning,'String'));

if nargin == 4 
    Tag = 'a';
    Label = '';
end
if nargin == 5
    Label = '';
end

hold on;
Center= x + y*sqrt(-1);
Radius = 5;
%Radius = 30./binning;
%Gridpt = 100;
%[u,v]=circle(Center,Radius,Gridpt);
%line(u,v,'LineWidth',1,'Color',[1 0 0]);
uu = [x x x x-Radius x+Radius];
vv = [y-Radius y+Radius y y y];
line(uu,vv,'LineWidth',2,'color',Color,'Tag',Tag);
if ~isempty(Label)
    text(x+Radius,y-Radius,Label,'Color',Color,'Tag',Tag,'FontSize',14,'FontWeight','bold');
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_autopicker_template_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_autopicker_goodbadfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_autopicker_dir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function checkbox_autopicker_parallelmode_Callback(hObject, eventdata, handles)
function input_autopicker_anglestart_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_autopicker_anglestop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_autopicker_angleincrement_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_autopicker_anglestart_Callback(hObject, eventdata, handles)
function input_autopicker_anglestop_Callback(hObject, eventdata, handles)
function input_autopicker_angleincrement_Callback(hObject, eventdata, handles)
function input_autopicker_noparticles_Callback(hObject, eventdata, handles)
function input_autopicker_noparticles_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_autopicker_binning_Callback(hObject, eventdata, handles)
function input_autopicker_binning_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit11_Callback(hObject, eventdata, handles)
function edit11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_autopicker_align_Callback(hObject, eventdata, handles)
function input_autopicker_align_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_autopicker_template_Callback(hObject, eventdata, handles)
function input_autopicker_fftfilename_Callback(hObject, eventdata, handles)
function input_autopicker_fftfilename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function checkbox_reusefftfile_Callback(hObject, eventdata, handles)
function edit_autopicker_maskstackfile_Callback(hObject, eventdata, handles)
function edit_autopicker_maskstackfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_autopicker_imagename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider_autopicker_file_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



