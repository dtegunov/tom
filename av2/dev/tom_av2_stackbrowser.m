function varargout = tom_av2_stackbrowser(varargin)
%GUI for browsing 2D particle stacks 
%
%SYNTAX
%tom_av2_stackbrowser(stackfile, alignfile)
%
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
%Created: 12/01/06 AK

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_av2_stackbrowser_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_stackbrowser_OutputFcn, ...
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
function tom_av2_stackbrowser_OpeningFcn(hObject, eventdata, handles, varargin)

global storage_av2_stackbrowser;

storage_av2_stackbrowser.stack = [];
storage_av2_stackbrowser.align = [];
storage_av2_stackbrowser.numcols = 10;
storage_av2_stackbrowser.numimages = 100;
storage_av2_stackbrowser.pos = 1;
storage_av2_stackbrowser.imsize = [0 0 0];
storage_av2_stackbrowser.page = 0;
storage_av2_stackbrowser.display.marks = 1;
storage_av2_stackbrowser.display.ccc = 1;
storage_av2_stackbrowser.display.quality = 1;
storage_av2_stackbrowser.quality = 5;
storage_av2_stackbrowser.fontsize = 10;
storage_av2_stackbrowser.bandpass.enable = 0;
storage_av2_stackbrowser.bandpass.low = 0;
storage_av2_stackbrowser.bandpass.high = 0;
storage_av2_stackbrowser.makerefmode = 0;

set(findobj('Tag','input_quality'),'String',num2str(storage_av2_stackbrowser.quality));

axes(findobj('Tag','imchoosebox'));axis off;

%load stack if given as input parameter
if nargin > 3
    storage_av2_stackbrowser.stackfile = varargin{1};
    set(findobj('Tag','input_particlestack'),'String',varargin{1});
    
    if size(varargin,2) == 1 | (size(varargin,2) == 1 & strcmp(varargin{2},'makerefmode') == 1)
        storage_av2_stackbrowser.align = tom_av2_create_alignfromstack(storage_av2_stackbrowser.stackfile);
        storage_av2_stackbrowser.alignmentfile = '';
    else
        if strcmp(varargin{2},'makerefmode') == 1
            storage_av2_stackbrowser.alignmentfile = '';
            storage_av2_stackbrowser.makerefmode = 1;
            set(findobj('Tag','panel_filter'),'Visible','off');
            set(findobj('Tag','panel_save'),'Visible','off');
            set(findobj('Tag','display_ccc'),'Value',0);
            set(findobj('Tag','display_quality'),'Value',0);
            set(findobj('Tag','button_stackbrowse_makeref'),'Visible','on');
            storage_av2_stackbrowser.display.ccc = 0;
            storage_av2_stackbrowser.display.quality = 0;
            storage_av2_stackbrowser.align = tom_av2_create_alignfromstack(storage_av2_stackbrowser.stackfile);
        else
            storage_av2_stackbrowser.alignmentfile = varargin{2};
            s = load(storage_av2_stackbrowser.alignmentfile);
            storage_av2_stackbrowser.align = s.align2d;
            set(findobj('Tag','input_alignmentfile'),'String',varargin{2});
        end
    end

    for i=1:size(storage_av2_stackbrowser.align,2)
         storage_av2_stackbrowser.align(i).partno = i;
    end


    if storage_av2_stackbrowser.makerefmode == 0
        %select all particles by default
        for i=1:size(storage_av2_stackbrowser.align,2)
            storage_av2_stackbrowser.align(i).selected = 1;
        end
    else
        for i=1:size(storage_av2_stackbrowser.align,2)
            storage_av2_stackbrowser.align(i).selected = 0;
        end
    end
    
    loadstack(1);
    numpages = floor(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages);
    if numpages > 1
        set(findobj('Tag','slider_particles'),'Max',numpages,'SliderStep',[1./numpages 1./numpages]);
    else
        set(findobj('Tag','slider_particles'),'Enable','off');
    end
    set(findobj('Tag','text_pageinfo'),'String',['page 1 of ' num2str(ceil(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages))]);
    set(findobj('Tag','button_loadstack'),'Enable','off');

    render_image();
end

handles.output = hObject;

guidata(hObject, handles);

if storage_av2_stackbrowser.makerefmode == 1
    uiwait(handles.figure1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_av2_stackbrowser_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

global storage_av2_stackbrowser;

if storage_av2_stackbrowser.makerefmode == 1
    close(handles.figure1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directy input particle stack                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_particlestack_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse particle stack                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_particlestackbrowse_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

[filename, pathname] = uigetfile({'*.em'}, 'Pick an image stack');
if ischar(filename)
    if tom_isemfile([pathname '/' filename])
        set(findobj('Tag','input_particlestack'),'String',[pathname '/' filename]);
        storage_av2_stackbrowser.stackfile = [pathname '/' filename];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directy input alignment file                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_alignmentfile_Callback(hObject, eventdata, handles)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse alignment file                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_alignmentfilebrowse_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

[filename, pathname] = uigetfile({'*.mat'}, 'Pick an alignment file');
if ischar(filename)
    if exist([pathname '/' filename],'file')
        set(findobj('Tag','input_alignmentfile'),'String',[pathname '/' filename]);
        storage_av2_stackbrowser.alignmentfile = [pathname '/' filename];
    end
end

s = load(storage_av2_stackbrowser.alignmentfile);
storage_av2_stackbrowser.align = s.align2d;

%select all particles by default
for i=1:size(storage_av2_stackbrowser.align,2)
    storage_av2_stackbrowser.align(i).selected = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  delete selected particles                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_particlesdelete_Callback(hObject, eventdata, handles)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  invert selection                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_invertselection_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

%draw marks
delete(findobj('-regexp','Tag','rect_*'));
if storage_av2_stackbrowser.display.marks == 1
    numrows = storage_av2_stackbrowser.numimages ./ storage_av2_stackbrowser.numcols;
    for row=1:numrows
        pos_y = (row-1).*storage_av2_stackbrowser.imsize(2)+storage_av2_stackbrowser.imsize(2).*0.1;
        for column=1:storage_av2_stackbrowser.numcols
            abspartnumber = storage_av2_stackbrowser.page.*storage_av2_stackbrowser.numimages+(row-1).*storage_av2_stackbrowser.numcols+column;
            if abspartnumber <= size(storage_av2_stackbrowser.align,2)
                %mark particle green if selected and red if not selected
                if storage_av2_stackbrowser.align(abspartnumber).selected == 1
                    storage_av2_stackbrowser.align(abspartnumber).selected = 0;
                    color = [1 0 0];
                else
                    storage_av2_stackbrowser.align(abspartnumber).selected = 1;
                    color = [0 1 0];
                end
                pos_x = (column-1).*storage_av2_stackbrowser.imsize(1)+storage_av2_stackbrowser.imsize(1).*0.1;
                r = rectangle('Position',[pos_x,pos_y,storage_av2_stackbrowser.imsize(1).*0.1,storage_av2_stackbrowser.imsize(2).*0.1],'FaceColor',color,'Tag',['rect_' abspartnumber]);
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  slider callback                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider_particles_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

set(findobj('Tag','text_pageinfo'),'String',['page ' num2str(round(get(hObject,'Value'))+1) ' of ' num2str(ceil(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages))]);

position = round(get(hObject,'Value')).*storage_av2_stackbrowser.numimages+1;
storage_av2_stackbrowser.page = round(get(hObject,'Value'));
loadstack(position);
drawnow;
render_image();
drawnow;
set(findobj('Tag','input_nocols'),'Enable','off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load stack                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_loadstack_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

if isempty(storage_av2_stackbrowser.stackfile)
    errordlg('Load stack file first.');
    return;   
end

if isempty(storage_av2_stackbrowser.align)
    storage_av2_stackbrowser.align = tom_av2_create_alignfromstack(storage_av2_stackbrowser.stackfile);
    %select all particles by default
    for i=1:size(storage_av2_stackbrowser.align,2)
        storage_av2_stackbrowser.align(i).selected = 1;
    end
end

for i=1:size(storage_av2_stackbrowser.align,2)
         storage_av2_stackbrowser.align(i).partno = i;
end

loadstack(1);
numpages = floor(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages);
if numpages > 1
    set(findobj('Tag','slider_particles'),'Max',numpages,'SliderStep',[1./numpages 1./numpages]);
else
    set(findobj('Tag','slider_particles'),'Enable','off');
end
set(findobj('Tag','text_pageinfo'),'String',['page 1 of ' num2str(ceil(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages))]);
set(hObject,'Enable','off');

render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  number of columns                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_nocols_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

val = str2num(get(hObject,'String'));

if ~isempty(val) && val > 1
    storage_av2_stackbrowser.numcols = round(val);
    storage_av2_stackbrowser.numimages = storage_av2_stackbrowser.numcols.^2;
end
loadstack(storage_av2_stackbrowser.pos);
render_image();
numpages = floor(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages);
set(findobj('Tag','slider_particles'),'Max',numpages,'SliderStep',[1./numpages 1./numpages]);
set(findobj('Tag','text_pageinfo'),'String',['page 1 of ' num2str(ceil(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages))]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display good/bad marks                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_marks_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

storage_av2_stackbrowser.display.marks = get(hObject,'Value');
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display CCC values                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_ccc_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

storage_av2_stackbrowser.display.ccc = get(hObject,'Value');
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display quality                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_quality_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

storage_av2_stackbrowser.display.quality = get(hObject,'Value');
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  save to stack                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_savetostack_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;


for i=1:size(storage_av2_stackbrowser.align,2)
    if storage_av2_stackbrowser.align(i).selected == 1
        storage_av2_stackbrowser.align(i).quality = storage_av2_stackbrowser.quality;
    end
end
align2d = storage_av2_stackbrowser.align;
align2d = rmfield(align2d,'selected');
align2d = rmfield(align2d,'partno');

outfile = get(findobj('Tag','input_alignmentfile'),'String');
if isempty(outfile)
    [filename, pathname] = uigetfile({'*.em'}, 'Pick an output alignment file');
    if ischar(filename)
        outfile = [pathname '/' filename];
        save(outfile,'align2d');
    end
else
    save(outfile,'align2d');
end
disp('Alignment File saved.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  input quality                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_quality_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

quality = round(str2num(get(hObject,'String')));

if quality > 10 || quality < 0
    errordlg('quality must be between 0 and 10.');
    set(hObject,'String',num2str(storage_av2_stackbrowser.quality));
    return;
end

storage_av2_stackbrowser.quality = quality;
set(hObject,'String',num2str(storage_av2_stackbrowser.quality));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  histogram reset                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_parthisto_reset_Callback(hObject, eventdata, handles)

calc_histogram();
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  histogram set                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_parthisto_set_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

axes(findobj('Tag','partbrowser_histogram'));
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

storage_av2_stackbrowser.Scale=[x(1) x(2)];
render_image();
set(findobj('Tag','input_parthisto_low'),'String',x(1));
set(findobj('Tag','input_parthisto_high'),'String',x(2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  histogram set manually                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_parthisto_manually_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

min=str2num(get(findobj('Tag','input_parthisto_low'),'String'));
max=str2num(get(findobj('Tag','input_parthisto_high'),'String'));

if max>min
    storage_av2_stackbrowser.Scale = [min max];
    set(findobj('Tag','partbrowser_histogram'),'Xlim',[min max]);
else
    set(findobj('Tag','input_parthisto_low'),'String',num2str(storage_av2_stackbrowser.Scale(1)));
    set(findobj('Tag','input_parthisto_high'),'String',num2str(storage_av2_stackbrowser.Scale(2)));
end

render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bandpass low                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partbandpass_low_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

val = str2num(get(hObject,'String'));

if ~isempty(val) & val > 0
    storage_av2_stackbrowser.bandpass.low = round(val);
end

if storage_av2_stackbrowser.bandpass.enable == 1
    render_image();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bandpass high                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partbandpass_high_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

val = str2num(get(hObject,'String'));

if ~isempty(val) & val > 0
    storage_av2_stackbrowser.bandpass.high = round(val);
end

if storage_av2_stackbrowser.bandpass.enable == 1
    render_image();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bandpass enable                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_partbandpass_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

storage_av2_stackbrowser.bandpass.enable = get(hObject,'Value');

render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for output filtered stack                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_stack_browse_Callback(hObject, eventdata, handles)

[filename, pathname] = uiputfile({'*.em'}, 'Pick a stack file');
if ischar(filename)
    set(findobj('Tag','output_stack'),'String',[pathname '/' filename]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for output align file                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_outputalignmentfile_browse_Callback(hObject, eventdata, handles)

[filename, pathname] = uiputfile({'*.mat'}, 'Pick a mat file');
if ischar(filename)
    set(findobj('Tag','output_alignmentfile'),'String',[pathname '/' filename]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter stack                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_filterstack_Callback(hObject, eventdata, handles)

instack = get(findobj('Tag','input_particlestack'),'String');
inalign = get(findobj('Tag','input_alignmentfile'),'String');
outstack = get(findobj('Tag','output_stack'),'String');
outalign = get(findobj('Tag','output_alignmentfile'),'String');
value = str2num(get(findobj('Tag','filterstack_value'),'String'));

contents = get(findobj('Tag','filterstack_fieldname'),'String');
fieldname = contents{get(findobj('Tag','filterstack_fieldname'),'Value')};

contents = get(findobj('Tag','filterstack_operator'),'String');
operator = contents{get(findobj('Tag','filterstack_operator'),'Value')};

if isempty(outstack)
    errordlg('Select output stack file name!');
    return;
end

if isempty(outalign)
    errordlg('Select output alignment file name!');
    return;
end

if isempty(operator)
    errordlg('Input a relational operator!');
    return;
end

if isempty(operator)
    errordlg('Input a relational operator!');
    return;
end

if isempty(fieldname)
    errordlg('Input a fieldname!');
    return;
end

if isempty(value)
    errordlg('Input a value!');
    return;
end

try
    tom_av2_filterstack(instack,inalign,fieldname,operator,value,outstack,outalign);
catch
    errordlg('Filtering failed, please check input!');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  norm particles                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dropdown_normparticles_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

position = round(get(findobj('Tag','slider_particles'),'Value')).*storage_av2_stackbrowser.numimages+1;
loadstack(position);
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  make reference stack                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_stackbrowse_makeref_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

if storage_av2_stackbrowser.makerefmode == 1

    instackfile = get(findobj('Tag','input_particlestack'),'String');
    [filename, pathname] = uiputfile({'*.em'}, 'Save reference stack as');
    if ischar(filename)
        outstackfile = [pathname '/' filename];
    else
        return;
    end

    if tom_isemfile(instackfile) ~= 1
        error('Could not load stack file.');
    end

    h = tom_reademheader(instackfile);
    newstack = [];
    c = 1;

    stack_sz=0;
    %check out size of new stack
    for i=1:size(storage_av2_stackbrowser.align,2)
        if storage_av2_stackbrowser.align(i).selected == 1;
            stack_sz=stack_sz+1;
        end
    end

    if stack_sz > 0
        %allocate some memory
        newstack=zeros(h.Header.Size(1),h.Header.Size(2),stack_sz);

        zz=1;
        for i=1:size(storage_av2_stackbrowser.align,2)
            if storage_av2_stackbrowser.align(i).selected == 1
                im = tom_emreadc(instackfile,'subregion',[1 1 i],[h.Header.Size(1)-1 h.Header.Size(2)-1 0]);
                newstack(:,:,zz)=im.Value;
                zz=zz+1;
            end
        end

        tom_emwrite(outstackfile,newstack);
        handles.output = outstackfile;
        guidata(hObject, handles);
        uiresume(handles.figure1);
    else
        error('New stack would contain no particles, please check your input.');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  sort particles                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_sortgo_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

fieldnamecontents = get(findobj('Tag','dropdown_sortparticles_fieldname'),'String');
fieldname = fieldnamecontents{get(findobj('Tag','dropdown_sortparticles_fieldname'),'Value')};

fieldordercontents = get(findobj('Tag','dropdown_sortparticles_order'),'String');
fieldorder = fieldordercontents{get(findobj('Tag','dropdown_sortparticles_order'),'Value')};

storage_av2_stackbrowser.stack = tom_sortstruct(storage_av2_stackbrowser.stack,fieldname, fieldorder);

position = round(get(findobj('Tag','slider_particles'),'Value')).*storage_av2_stackbrowser.numimages+1;
loadstack(position);
render_image();


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
%%  calculate histogram                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc_histogram(im)

global storage_av2_stackbrowser;

if nargin < 1
    [mean max min std] = tom_dev(storage_av2_stackbrowser.stack,'noinfo');
else
    [mean max min std] = tom_dev(im,'noinfo');
end
min = mean-3.*std;
max = mean+3.*std;

[h,n] = tom_hist3d(storage_av2_stackbrowser.stack);
h = 200 .* h ./ (100.*size(storage_av2_stackbrowser.stack,1) .* size(storage_av2_stackbrowser.stack,2) .* size(storage_av2_stackbrowser.stack,3));
axesobj = findobj('Tag','partbrowser_histogram');
axes(axesobj); 
bar(n,h);
axis auto;
set(axesobj,'Tag','partbrowser_histogram');
set(findobj('Tag','input_parthisto_low'),'String',num2str(min));
set(findobj('Tag','input_parthisto_high'),'String',num2str(max));
storage_av2_stackbrowser.Scale = [min max];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load stack                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadstack(start)

global storage_av2_stackbrowser;

if ischar(start)
  if isequal(start, 'next')
    start = storage_av2_stackbrowser.pos + storage_av2_stackbrowser.numimages;
  elseif isequal(start, 'previous')
    start = storage_av2_stackbrowser.pos - storage_av2_stackbrowser.numimages;
  end
end

stop = storage_av2_stackbrowser.numimages;

header = tom_reademheader(storage_av2_stackbrowser.stackfile);
storage_av2_stackbrowser.imsize = header.Header.Size;

if start+stop > size(storage_av2_stackbrowser.align,2)
    missingpages = start+stop-size(storage_av2_stackbrowser.align,2)-1;
    stop = size(storage_av2_stackbrowser.align,2)-start+1;
end

storage_av2_stackbrowser.stack = zeros(header.Header.Size(1),header.Header.Size(2),stop-start);

for i=start:stop
     im = tom_emreadc(storage_av2_stackbrowser.stackfile,'subregion',[1 1 storage_av2_stackbrowser.align(i).partno], [header.Header.Size(1)-1 header.Header.Size(2)-1 0]);
     storage_av2_stackbrowser.stack(:,:,i) = im.Value;
end

%if start+stop > header.Header.Size(3)
%    missingpages = start+stop-header.Header.Size(3)-1;
%    stop = header.Header.Size(3)-start+1;
%end
%im = tom_emreadc(storage_av2_stackbrowser.stackfile,'subregion',[1 1 start],[header.Header.Size(1)-1 header.Header.Size(2)-1 stop-1]);
%storage_av2_stackbrowser.stack = im.Value;

storage_av2_stackbrowser.pos = start;

%norm stack
normcontents = get(findobj('Tag','dropdown_normparticles'),'String');
normval = normcontents{get(findobj('Tag','dropdown_normparticles'),'Value')};

if ~strcmp(normval,'none')
    for i=1:size(storage_av2_stackbrowser.stack,3)
        storage_av2_stackbrowser.stack(:,:,i) = tom_norm(storage_av2_stackbrowser.stack(:,:,i),normval);
    end
end

calc_histogram();

if storage_av2_stackbrowser.bandpass.high == 0
    set(findobj('Tag','partbandpass_high')','String',num2str(header.Header.Size(1)./2));
    storage_av2_stackbrowser.bandpass.high = header.Header.Size(1)./2;
end

if exist('missingpages') & missingpages > 0
    storage_av2_stackbrowser.stack = cat(3,storage_av2_stackbrowser.stack,zeros(header.Header.Size(1),header.Header.Size(2),missingpages));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display image                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_image()

global storage_av2_stackbrowser;

storage_av2_stackbrowser.fontsize = storage_av2_stackbrowser.imsize(1)./storage_av2_stackbrowser.numcols./1.3;

%bandpass filter
if storage_av2_stackbrowser.bandpass.enable == 1
    im = zeros(size(storage_av2_stackbrowser.stack));
    for i=1:size(storage_av2_stackbrowser.stack,3)
        im(:,:,i) = tom_bandpass(storage_av2_stackbrowser.stack(:,:,i),storage_av2_stackbrowser.bandpass.low,storage_av2_stackbrowser.bandpass.high);
    end
    scale = storage_av2_stackbrowser.Scale  ./(storage_av2_stackbrowser.imsize(1) .* storage_av2_stackbrowser.imsize(2));
    if  storage_av2_stackbrowser.bandpass.low > 0
            scale = scale - mean(mean(mean(storage_av2_stackbrowser.stack))) ./ (storage_av2_stackbrowser.imsize(1) .* storage_av2_stackbrowser.imsize(1));
    end
else
    im = storage_av2_stackbrowser.stack;
    scale = storage_av2_stackbrowser.Scale;
end
tmpobj = findobj('Tag','imchoosebox');
axes(tmpobj);
h = tom_dspcub(im, 0, storage_av2_stackbrowser.numcols,scale);
set(h,'buttonDownFcn',@mark_particle);

%draw ccc values
delete(findobj('Tag','cccstring'));
if storage_av2_stackbrowser.display.ccc == 1
    numrows = storage_av2_stackbrowser.numimages ./ storage_av2_stackbrowser.numcols;
    for row=1:numrows
        pos_y = (row-1).*storage_av2_stackbrowser.imsize(2)+storage_av2_stackbrowser.imsize(2).*0.9;
        for column=1:storage_av2_stackbrowser.numcols
            abspartnumber = storage_av2_stackbrowser.page.*storage_av2_stackbrowser.numimages+(row-1).*storage_av2_stackbrowser.numcols+column;
            if abspartnumber <= size(storage_av2_stackbrowser.align,2)
                cccval = storage_av2_stackbrowser.align(abspartnumber).ccc;
                pos_x = (column-1).*storage_av2_stackbrowser.imsize(1)+storage_av2_stackbrowser.imsize(1).*0.65;
                t = text(pos_x,pos_y,sprintf('%0.2f',double(cccval)));
                set(t,'Tag','cccstring','Color',[1 1 0],'FontSize',storage_av2_stackbrowser.fontsize,'FontWeight','bold');
            end
        end
    end
end

%draw quality values
delete(findobj('Tag','qualitystring'));
if storage_av2_stackbrowser.display.quality == 1
    numrows = storage_av2_stackbrowser.numimages ./ storage_av2_stackbrowser.numcols;
    for row=1:numrows
        pos_y = (row-1).*storage_av2_stackbrowser.imsize(2)+storage_av2_stackbrowser.imsize(2).*0.2;
        for column=1:storage_av2_stackbrowser.numcols
            abspartnumber = storage_av2_stackbrowser.page.*storage_av2_stackbrowser.numimages+(row-1).*storage_av2_stackbrowser.numcols+column;
            if abspartnumber <= size(storage_av2_stackbrowser.align,2)
                qualityval = storage_av2_stackbrowser.align(abspartnumber).quality;
                pos_x = (column-1).*storage_av2_stackbrowser.imsize(1)+storage_av2_stackbrowser.imsize(1).*0.8;
                t = text(pos_x,pos_y,sprintf('%i',qualityval));
                set(t,'Tag','qualitystring','Color',[1 1 0],'FontSize',storage_av2_stackbrowser.fontsize,'FontWeight','bold');
            end
        end
    end
end

%draw marks
delete(findobj('-regexp','Tag','rect_*'));
if storage_av2_stackbrowser.display.marks == 1
    numrows = storage_av2_stackbrowser.numimages ./ storage_av2_stackbrowser.numcols;
    for row=1:numrows
        pos_y = (row-1).*storage_av2_stackbrowser.imsize(2)+storage_av2_stackbrowser.imsize(2).*0.1;
        for column=1:storage_av2_stackbrowser.numcols
            abspartnumber = storage_av2_stackbrowser.page.*storage_av2_stackbrowser.numimages+(row-1).*storage_av2_stackbrowser.numcols+column;
            if abspartnumber <= size(storage_av2_stackbrowser.align,2)
                %mark particle green if selected and red if not selected
                if storage_av2_stackbrowser.align(abspartnumber).selected == 1
                    color = [0 1 0];
                else
                    color = [1 0 0];
                end
                pos_x = (column-1).*storage_av2_stackbrowser.imsize(1)+storage_av2_stackbrowser.imsize(1).*0.1;
                r = rectangle('Position',[pos_x,pos_y,storage_av2_stackbrowser.imsize(1).*0.1,storage_av2_stackbrowser.imsize(2).*0.1],'FaceColor',color,'Tag',['rect_' abspartnumber]);
            end
        end
    end
end

set(tmpobj,'Tag','imchoosebox');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mark particle                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mark_particle(a,b)

global storage_av2_stackbrowser;

point1=get(gca,'currentpoint');
button = get(gcf,'selectiontype');

%button values:
%normal: left mouse button
%alt: right mouse button
%extend: middle mouse buttons

pt = point1(1,1:2);
x1 = round(pt(1));
y1 = round(pt(2));
    
    %determine which particle was picked
    row = ceil(y1./storage_av2_stackbrowser.imsize(2));
    column = ceil(x1./storage_av2_stackbrowser.imsize(1));
    abspartnumber = storage_av2_stackbrowser.page.*storage_av2_stackbrowser.numimages+(row-1).*storage_av2_stackbrowser.numcols+column;

if strcmp(button,'normal') == true

    
    if abspartnumber <= size(storage_av2_stackbrowser.align,2)
        %toggle selection on particle
        if storage_av2_stackbrowser.align(abspartnumber).selected == 1
            storage_av2_stackbrowser.align(abspartnumber).selected = 0;
            color = [1 0 0];
        else
            storage_av2_stackbrowser.align(abspartnumber).selected = 1;
            color = [0 1 0];
        end

        %redraw mark
        delete(findobj('Tag',['rect_' abspartnumber]));
        pos_y = (row-1).*storage_av2_stackbrowser.imsize(2)+storage_av2_stackbrowser.imsize(2).*0.1;
        pos_x = (column-1).*storage_av2_stackbrowser.imsize(1)+storage_av2_stackbrowser.imsize(1).*0.1;
        r = rectangle('Position',[pos_x,pos_y,storage_av2_stackbrowser.imsize(1).*0.1,storage_av2_stackbrowser.imsize(2).*0.1],'FaceColor',color,'Tag',['rect_' abspartnumber]);
    end

elseif strcmp(button,'alt') == true
    
    string = ['Filename: ' storage_av2_stackbrowser.align(abspartnumber).filename];
    string = strvcat(string, ['Position: X ' num2str(storage_av2_stackbrowser.align(abspartnumber).position.x) ', Y ' num2str(storage_av2_stackbrowser.align(abspartnumber).position.y)]);
    string = strvcat(string, ['Class: ' storage_av2_stackbrowser.align(abspartnumber).class]);
    string = strvcat(string, ['Shift: X ' num2str(storage_av2_stackbrowser.align(abspartnumber).shift.x) ', Y ' num2str(storage_av2_stackbrowser.align(abspartnumber).shift.y)]);
    string = strvcat(string, ['Angle: ' num2str(storage_av2_stackbrowser.align(abspartnumber).angle)]);
    string = strvcat(string, ['Isaligned? ', num2str(storage_av2_stackbrowser.align(abspartnumber).isaligned)]);
    string = strvcat(string, ['Dataset: ', storage_av2_stackbrowser.align(abspartnumber).dataset]);
    string = strvcat(string, ['Number in stack: ', num2str(abspartnumber)]);
    msgbox(string);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_particlestack_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_alignmentfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider_particles_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function input_quality_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_parthisto_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_parthisto_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_parthisto_low_Callback(hObject, eventdata, handles)
function input_parthisto_high_Callback(hObject, eventdata, handles)
function input_norows_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_nocols_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function partbandpass_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function partbandpass_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function output_filteredstack_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function filterstack_fieldname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function filterstack_operator_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function filterstack_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function filterstack_operator_Callback(hObject, eventdata, handles)
function filterstack_fieldname_Callback(hObject, eventdata, handles)
function filterstack_value_Callback(hObject, eventdata, handles)
function output_stack_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function output_alignmentfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function output_alignmentfile_Callback(hObject, eventdata, handles)
function output_stack_Callback(hObject, eventdata, handles)
function dropdown_normparticles_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end