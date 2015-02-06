function varargout = tom_av3_stackbrowser(varargin)
%TOM_AV3_STACKBROWSER is a GUI for browsing 3D particle stacks
%
%   varargout = tom_av3_stackbrowser(varargin)
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
%   ... = tom_av3_stackbrowser(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 05/23/06
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
                   'gui_OpeningFcn', @tom_av3_stackbrowser_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av3_stackbrowser_OutputFcn, ...
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
function tom_av3_stackbrowser_OpeningFcn(hObject, eventdata, handles, varargin)

global storage_av3_stackbrowser;

axes(findobj('Tag','av3_stackbrowser_particles'));axis off;

storage_av3_stackbrowser.numcols = 20;
storage_av3_stackbrowser.numrows = 20;
storage_av3_stackbrowser.resamplezval = 4;
storage_av3_stackbrowser.resamplexyval = 2;

storage_av3_stackbrowser.bandpass.high = 0;
storage_av3_stackbrowser.bandpass.low = 0;
storage_av3_stackbrowser.bandpass.enable = 0;
storage_av3_stackbrowser.orientation = 'xy';

storage_av3_stackbrowser.display.marks = 1;
storage_av3_stackbrowser.display.ccc = 1;
storage_av3_stackbrowser.display.quality = 1;
storage_av3_stackbrowser.display.classno = 1;
storage_av3_stackbrowser.display.projclassno = 1;

storage_av3_stackbrowser.refstep = 0;
storage_av3_stackbrowser.page = 0;


if size(varargin,1) == 0
    [filename,pathname] = uigetfile('*.mat','Select a 3D alignment file');

    if isequal(pathname,0)
        error('Cancel button pressed. No data loaded.');
        return;
    end;
    storage_av3_stackbrowser.alignfile = [pathname '/' filename];
else
    storage_av3_stackbrowser.alignfile = varargin{1};
end


%open the alignment file
s = load(storage_av3_stackbrowser.alignfile);
storage_av3_stackbrowser.align = s.Align;

if ~isfield(storage_av3_stackbrowser.align,'quality')
    for i=1:size(storage_av3_stackbrowser.align,2)
        for k=1:size(storage_av3_stackbrowser.align,1)
            storage_av3_stackbrowser.align(k,i).quality = 0;
        end
    end
end
    
%select all particles by default
for i=1:size(storage_av3_stackbrowser.align,2)
    for j=1:size(storage_av3_stackbrowser.align,1)
        storage_av3_stackbrowser.align(j,i).selected = 1;
    end
end

for i=1:size(storage_av3_stackbrowser.align,2)
    for j=1:size(storage_av3_stackbrowser.align,1)
        storage_av3_stackbrowser.align(j,i).partno = i;
    end
end

storage_av3_stackbrowser.numpages = floor(size(storage_av3_stackbrowser.align,2) ./ storage_av3_stackbrowser.numrows);
if storage_av3_stackbrowser.numpages > 1
    set(findobj('Tag','slider_av3_stackbrowser_particles'),'Max',storage_av3_stackbrowser.numpages,'SliderStep',[1./storage_av3_stackbrowser.numpages 1./storage_av3_stackbrowser.numpages],'Enable','on');
else
    set(findobj('Tag','slider_av3_stackbrowser_particles'),'Enable','off');
end
set(findobj('Tag','text_av3_stackbrowser_pageinfo'),'String',['page 1 of ' num2str(ceil(size(storage_av3_stackbrowser.align,2)./storage_av3_stackbrowser.numrows))]);
set(findobj('Tag','button_loadstack'),'Enable','off');

%update refinement step dropdown
refstring = {'0'};
for i=1:size(storage_av3_stackbrowser.align,1)
    refstring{end+1} = num2str(i);
end
set(findobj('Tag','menu_av3_stackbrowser_refstep'),'String',refstring);



loadstack(1);
render_image();

% Choose default command line output for tom_av3_stackbrowser
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_av3_stackbrowser wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_av3_stackbrowser_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  particles slider                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider_av3_stackbrowser_particles_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

set(findobj('Tag','text_av3_stackbrowser_pageinfo'),'String',['page ' num2str(round(get(hObject,'Value'))+1) ' of ' num2str(ceil(size(storage_av3_stackbrowser.align,2)./storage_av3_stackbrowser.numrows))]);

position = round(get(hObject,'Value')).*storage_av3_stackbrowser.numrows+1;
storage_av3_stackbrowser.page = round(get(hObject,'Value'));
loadstack(position);
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  histogram reset                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_stackbrowser_histogramreset_Callback(hObject, eventdata, handles)

calc_histogram();
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  histogram set                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_stackbrowser_histogramset_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

axes(findobj('Tag','av3_stackbrowser_histogram'));
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

storage_av3_stackbrowser.Scale=[x(1) x(2)];
render_image();
set(findobj('Tag','edit_av3_stackbrowser_histogramlow'),'String',x(1));
set(findobj('Tag','edit_av3_stackbrowser_histogramhigh'),'String',x(2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  histogram set manually                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_stackbrowser_histogramsetmanually_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

min=str2num(get(findobj('Tag','edit_av3_stackbrowser_histogramlow'),'String'));
max=str2num(get(findobj('Tag','edit_av3_stackbrowser_histogramhigh'),'String'));

if max>min
    storage_av3_stackbrowser.Scale = [min max];
    set(findobj('Tag','av3_stackbrowser_histogram'),'Xlim',[min max]);
else
    set(findobj('Tag','edit_av3_stackbrowser_histogramlow'),'String',num2str(storage_av3_stackbrowser.Scale(1)));
    set(findobj('Tag','edit_av3_stackbrowser_histogramhigh'),'String',num2str(storage_av3_stackbrowser.Scale(2)));
end

render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  norm                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function menu_av3_stackbrowser_norm_Callback(hObject, eventdata, handles)

loadstack();
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bandpass low                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_av3_stackbrowser_bandpasslow_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

storage_av3_stackbrowser.bandpass.low = str2num(get(hObject,'String'));
loadstack();
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bandpass high                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_av3_stackbrowser_bandpasshigh_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

storage_av3_stackbrowser.bandpass.high = str2num(get(hObject,'String'));
loadstack();
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bandpass enable                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_av3_stackbrowser_bandpassenable_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

storage_av3_stackbrowser.bandpass.enable = get(hObject,'Value');
loadstack();
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display marks                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_av3_stackbrowser_display_marks_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

storage_av3_stackbrowser.display.marks = get(hObject,'Value');
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display ccc                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_av3_stackbrowser_display_ccc_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

storage_av3_stackbrowser.display.ccc = get(hObject,'Value');
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display quality                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_av3_stackbrowser_display_quality_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

storage_av3_stackbrowser.display.quality = get(hObject,'Value');
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display class no                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_av3_stackbrowser_display_classno_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

storage_av3_stackbrowser.display.classno = get(hObject,'Value');
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display projection class no                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_av3_stackbrowser_display_projclassno_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

storage_av3_stackbrowser.display.projclassno = get(hObject,'Value');
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  resampling value in z                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_av3_stackbrowser_resamplez_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

storage_av3_stackbrowser.resamplezval = str2num(get(hObject,'String'));
loadstack();
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  resampling value in xy                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_av3_stackbrowser_resamplexy_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

storage_av3_stackbrowser.resamplexyval = str2num(get(hObject,'String'));
loadstack();
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  resampling value in xy                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function menu_av3_stackbrowser_orientation_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

contents = get(hObject,'String');
val = contents{get(hObject,'Value')};

storage_av3_stackbrowser.orientation = val;

loadstack();
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  refinement step                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function menu_av3_stackbrowser_refstep_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

contents = get(hObject,'String');
val = contents{get(hObject,'Value')};

storage_av3_stackbrowser.refstep = str2num(val);

loadstack();
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  sort particles                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_stackbrowser_sortgo_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

fieldnamecontents = get(findobj('Tag','menu_av3_stackbrowser_sortfield'),'String');
fieldname = fieldnamecontents{get(findobj('Tag','menu_av3_stackbrowser_sortfield'),'Value')};

if strcmp(fieldname,'CCC')
    fieldname = 'CCC';
elseif strcmp(fieldname,'Class')
    fieldname = 'Class';
elseif strcmp(fieldname,'ProjClass')
    fieldname = 'ProjectionClass';
elseif strcmp(fieldname,'quality')
    fieldname = 'quality';
else
    fieldname = 'partno';
end

fieldordercontents = get(findobj('Tag','menu_av3_stackbrowser_sortorder'),'String');
fieldorder = fieldordercontents{get(findobj('Tag','menu_av3_stackbrowser_sortorder'),'Value')};

if strcmp(fieldorder,'ASC')
    fieldorder = 'ascend';
else
    fieldorder = 'descend';
end

storage_av3_stackbrowser.align = tom_sortstruct3d(storage_av3_stackbrowser.align,fieldname, fieldorder,storage_av3_stackbrowser.refstep);

loadstack(1);
render_image();
set(findobj('Tag','slider_particles'),'Value',0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  invert selection                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_stackbrowser_invertselection_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

%draw marks
delete(findobj('-regexp','Tag','rect_*'));
if storage_av3_stackbrowser.display.marks == 1

    for row=1:storage_av3_stackbrowser.numrows
        pos_y = (row-1).*storage_av3_stackbrowser.imsize(2)+storage_av3_stackbrowser.imsize(2).*0.1;
        
            abspartnumber = storage_av3_stackbrowser.page.*storage_av3_stackbrowser.numrows+row;
            if abspartnumber <= size(storage_av3_stackbrowser.align,2)
                %mark particle green if selected and red if not selected
                if storage_av3_stackbrowser.align(abspartnumber).selected == 1
                    storage_av3_stackbrowser.align(abspartnumber).selected = 0;
                    color = [1 0 0];
                else
                    storage_av3_stackbrowser.align(abspartnumber).selected = 1;
                    color = [0 1 0];
                end
                pos_x = storage_av3_stackbrowser.imsize(1).*0.1;
                r = rectangle('Position',[pos_x,pos_y,storage_av3_stackbrowser.imsize(1).*0.2,storage_av3_stackbrowser.imsize(2).*0.2],'FaceColor',color,'Tag',['rect_' abspartnumber]);
            end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  save selection                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_stackbrowser_saveselection_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

storage_av3_stackbrowser.quality = str2num(get(findobj('Tag','input_av3_stackbrowser_newquality'),'String'));

j = size(storage_av3_stackbrowser.align,1);
for i=1:size(storage_av3_stackbrowser.align,2)
    if storage_av3_stackbrowser.align(1,i).selected == 1
        storage_av3_stackbrowser.align(j,i).quality = storage_av3_stackbrowser.quality;
    end
end
Align = storage_av3_stackbrowser.align;
Align = rmfield(Align,'selected');
Align = rmfield(Align,'partno');

outfile = storage_av3_stackbrowser.alignfile;
if isempty(outfile)
    [filename, pathname] = uigetfile({'*.em'}, 'Pick an output alignment file');
    if ischar(filename)
        outfile = [pathname '/' filename];
        save(outfile,'Align');
    end
else
    save(outfile,'Align');
end
disp('Alignment File saved.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for output directory                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_stackbrowser_browseoutdir_Callback(hObject, eventdata, handles)

pathname = uigetdir('Pick an output directory');
if ischar(pathname)
    set(findobj('Tag','input_av3_stackbrowser_outdir'),'String',pathname);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for output alignment file                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_stackbrowser_browseoutalign_Callback(hObject, eventdata, handles)

[filename, pathname] = uiputfile({'*.mat'}, 'Pick a mat file');
if ischar(filename)
    set(findobj('Tag','input_av3_stackbrowser_outalign'),'String',[pathname '/' filename]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter stack                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_stackbrowser_filterstack_Callback(hObject, eventdata, handles)

global storage_av3_stackbrowser;

outdir = get(findobj('Tag','input_av3_stackbrowser_outdir'),'String');
outalign = get(findobj('Tag','input_av3_stackbrowser_outalign'),'String');
value = str2num(get(findobj('Tag','input_av3_stackbrowser_filterval'),'String'));

contents = get(findobj('Tag','menu_av2_stackbrowser_filterfield'),'String');
fieldname = contents{get(findobj('Tag','menu_av2_stackbrowser_filterfield'),'Value')};

contents = get(findobj('Tag','menu_av2_stackbrowser_filteroperator'),'String');
operator = contents{get(findobj('Tag','menu_av2_stackbrowser_filteroperator'),'Value')};

if isempty(outdir)
    errordlg('Select output directory name!');
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

if ~exist(outdir,'dir')
    mkdir(outdir);
end

%Align = struct();
al = storage_av3_stackbrowser.align;
al(1,size(al,2)).Tempfilename = '';
al = rmfield(al,'partno');
al = rmfield(al,'selected');
lauf=1;
j = size(storage_av3_stackbrowser.align,1);
for k=1:size(storage_av3_stackbrowser.align,2)
    t = eval(['storage_av3_stackbrowser.align(j,k).(fieldname) ' operator ' value;']);
    if t==1
       for h=1:j
           Align(h,lauf) = al(h,k);
       end
       
       [pathstr, name, ext] = fileparts(al(j,k).Filename);
       
       %unix(['cp ' al(j,k).Filename ' ' outdir '/particle_' num2str(lauf) ext]);
       
       Align(j,lauf).Tempfilename = al(j,k).Filename;
       Align(j,lauf).Filename = [outdir '/particle_' num2str(lauf) ext];
       Align(j,lauf).Filename = Align(j,lauf).Tempfilename ;
       lauf=lauf+1; 
    end
end
save(outalign,'Align')

disp(' ');

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

global storage_av3_stackbrowser;

if(get(findobj('Tag','constant_scale'),'Value')==0)
    if nargin < 1
        [mean max min std] = tom_dev(storage_av3_stackbrowser.stack,'noinfo');
    else
        [mean max min std] = tom_dev(im,'noinfo');
    end
    min = mean-3.*std;
    max = mean+3.*std;

    [h,n] = tom_hist3d(storage_av3_stackbrowser.stack);
    h = 200 .* h ./ (100.*size(storage_av3_stackbrowser.stack,1) .* size(storage_av3_stackbrowser.stack,2) .* size(storage_av3_stackbrowser.stack,3) .* size(storage_av3_stackbrowser.stack,4));
    axesobj = findobj('Tag','av3_stackbrowser_histogram');
    axes(axesobj); 
    bar(n,h);
    axis auto;
    set(axesobj,'Tag','av3_stackbrowser_histogram');
    set(findobj('Tag','edit_av3_stackbrowser_histogramlow'),'String',num2str(min));
    set(findobj('Tag','edit_av3_stackbrowser_histogramhigh'),'String',num2str(max));
    storage_av3_stackbrowser.Scale = [min max];
end
%storage_av3_stackbrowser.Scale = [mean-1.*std mean+1.*std];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load stack                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadstack(start)

global storage_av3_stackbrowser;

if nargin == 0
    start = storage_av3_stackbrowser.pos;
end
    

if ischar(start)
  if isequal(start, 'next')
    start = storage_av3_stackbrowser.pos + storage_av3_stackbrowser.numrows;
  elseif isequal(start, 'previous')
    start = storage_av3_stackbrowser.pos - storage_av3_stackbrowser.numrows;
  end
end

stop = storage_av3_stackbrowser.numrows;

if start+stop > size(storage_av3_stackbrowser.align,2)
    %missingpages = start+stop-size(storage_av3_stackbrowser.align,2)-1;
    stop = size(storage_av3_stackbrowser.align,2)-start+1;
end

header = tom_reademheader(storage_av3_stackbrowser.align(1,1).Filename);
storage_av3_stackbrowser.imsize(1) = floor(header.Header.Size(1)./storage_av3_stackbrowser.resamplexyval);
storage_av3_stackbrowser.imsize(2) = floor(header.Header.Size(2)./storage_av3_stackbrowser.resamplexyval);
storage_av3_stackbrowser.imsize(3) = floor(header.Header.Size(3)./storage_av3_stackbrowser.resamplezval);

storage_av3_stackbrowser.numcols = storage_av3_stackbrowser.imsize(3);

storage_av3_stackbrowser.stack = zeros(storage_av3_stackbrowser.imsize(1),storage_av3_stackbrowser.imsize(2),storage_av3_stackbrowser.imsize(3),storage_av3_stackbrowser.numrows);

lauf=1;
for i=start:start+stop-1

     if storage_av3_stackbrowser.refstep == 0
         if strcmp(storage_av3_stackbrowser.orientation,'xy')
             im = tom_emreadc(storage_av3_stackbrowser.align(1,i).Filename,'resample',[storage_av3_stackbrowser.resamplexyval storage_av3_stackbrowser.resamplexyval storage_av3_stackbrowser.resamplezval]);
         elseif strcmp(storage_av3_stackbrowser.orientation,'yz')
             im = tom_emreadc(storage_av3_stackbrowser.align(1,i).Filename,'resample',[storage_av3_stackbrowser.resamplezval storage_av3_stackbrowser.resamplexyval storage_av3_stackbrowser.resamplexyval]);
             im.Value = shiftdim(im.Value,1);
         else
             im = tom_emreadc(storage_av3_stackbrowser.align(1,i).Filename,'resample',[storage_av3_stackbrowser.resamplexyval storage_av3_stackbrowser.resamplezval storage_av3_stackbrowser.resamplexyval]);
             im.Value = shiftdim(im.Value,2);
         end
     else
         im = tom_emreadc(storage_av3_stackbrowser.align(1,i).Filename);

         
         shifts = [];
         rotations = [];
%          for j=1:storage_av3_stackbrowser.refstep
%             shifts = [shifts; storage_av3_stackbrowser.align(j,i).Shift.X storage_av3_stackbrowser.align(j,i).Shift.Y storage_av3_stackbrowser.align(j,i).Shift.Z];
%             rotations = [rotations; storage_av3_stackbrowser.align(j,i).Angle.Phi storage_av3_stackbrowser.align(j,i).Angle.Psi storage_av3_stackbrowser.align(j,i).Angle.Theta];
%          end
%          [euler_out shift_out rott]=tom_sum_rotation(rotations,shifts);
         %rotate
         j=storage_av3_stackbrowser.refstep; 
         shift_out=[storage_av3_stackbrowser.align(j,i).Shift.X storage_av3_stackbrowser.align(j,i).Shift.Y storage_av3_stackbrowser.align(j,i).Shift.Z].*-1;
         euler_out = [storage_av3_stackbrowser.align(j,i).Angle.Psi storage_av3_stackbrowser.align(j,i).Angle.Phi  storage_av3_stackbrowser.align(j,i).Angle.Theta].*-1;
         
         if ~(shift_out(1) == 0 && shift_out(2) == 0 && shift_out(3) == 0)
             im.Value = tom_shift(im.Value,shift_out');
         end
         
         if ~(euler_out(1) == 0 && euler_out(2) == 0 && euler_out(3) == 0)
             im.Value = tom_rotate(im.Value,euler_out,'linear');
         end
         %shift
        
         
         if strcmp(storage_av3_stackbrowser.orientation,'yz')
            im.Value = shiftdim(im.Value,1);
         elseif strcmp(storage_av3_stackbrowser.orientation,'xz')
            im.Value = shiftdim(im.Value,2);
         end

         %bin
         im.Value = single(im.Value);
         im.Value  = tom_bininc(im.Value,[storage_av3_stackbrowser.resamplexyval storage_av3_stackbrowser.resamplexyval storage_av3_stackbrowser.resamplezval]);

     end

     storage_av3_stackbrowser.stack(:,:,:,lauf) = im.Value;
     lauf = lauf+1;
     
     
end

storage_av3_stackbrowser.pos = start;

%norm stack
normcontents = get(findobj('Tag','menu_av3_stackbrowser_norm'),'String');
normval = normcontents{get(findobj('Tag','menu_av3_stackbrowser_norm'),'Value')};

if ~strcmp(normval,'none')
    for i=1:size(storage_av3_stackbrowser.stack,4)
        storage_av3_stackbrowser.stack(:,:,:,i) = tom_norm(storage_av3_stackbrowser.stack(:,:,:,i),normval);
    end
end

calc_histogram();

if storage_av3_stackbrowser.bandpass.high == 0
    set(findobj('Tag','edit_av3_stackbrowser_bandpasshigh')','String',num2str(header.Header.Size(1)./2));
    storage_av3_stackbrowser.bandpass.high = header.Header.Size(1)./2;
end

%if exist('missingpages') && missingpages > 0
%    storage_av3_stackbrowser.stack = cat(4,storage_av3_stackbrowser.stack,zeros(storage_av3_stackbrowser.imsize(1),storage_av3_stackbrowser.imsize(2),storage_av3_stackbrowser.imsize(3),missingpages));
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display image                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_image()

global storage_av3_stackbrowser;

storage_av3_stackbrowser.fontsize = 12;

%bandpass filter
 
%if(~get(findobj('Tag','constant_scale'),'Value')) % always the same scale
if storage_av3_stackbrowser.bandpass.enable == 1
    im = zeros(size(storage_av3_stackbrowser.stack));
    for i=1:size(storage_av3_stackbrowser.stack,4)
        im(:,:,:,i) = tom_bandpass(storage_av3_stackbrowser.stack(:,:,:,i),storage_av3_stackbrowser.bandpass.low,storage_av3_stackbrowser.bandpass.high);
    end
    scale = storage_av3_stackbrowser.Scale  ./(storage_av3_stackbrowser.imsize(1) .* storage_av3_stackbrowser.imsize(2) .* storage_av3_stackbrowser.imsize(3));
    if  storage_av3_stackbrowser.bandpass.low > 0
            scale = scale - mean(mean(mean(mean(storage_av3_stackbrowser.stack)))) ./ (storage_av3_stackbrowser.imsize(1) .* storage_av3_stackbrowser.imsize(2) .* storage_av3_stackbrowser.imsize(3));
    end
else
    im = storage_av3_stackbrowser.stack;
    scale = storage_av3_stackbrowser.Scale;
end
tmpobj = findobj('Tag','av3_stackbrowser_particles');
axes(tmpobj);
im = reshape(im,storage_av3_stackbrowser.imsize(1),storage_av3_stackbrowser.imsize(2),[]);
h = tom_dspcub(im, 0, storage_av3_stackbrowser.numcols,scale);
set(h,'buttonDownFcn',@mark_particle);

%draw ccc values
delete(findobj('Tag','cccstring2'));
if storage_av3_stackbrowser.display.ccc == 1
    for row=1:storage_av3_stackbrowser.numrows
        pos_y = (row-1).*storage_av3_stackbrowser.imsize(2)+storage_av3_stackbrowser.imsize(2).*0.9;
        abspartnumber = ((storage_av3_stackbrowser.page).*storage_av3_stackbrowser.numrows)+row;
        if abspartnumber <= size(storage_av3_stackbrowser.align,2)
            if storage_av3_stackbrowser.refstep == 0
                refstep = storage_av3_stackbrowser.refstep+1;
            else
                refstep = storage_av3_stackbrowser.refstep;
            end
            cccval = storage_av3_stackbrowser.align(refstep,abspartnumber).CCC;
            pos_x = storage_av3_stackbrowser.imsize(1).*0.4;
            t = text(pos_x,pos_y,sprintf('%0.2f',double(cccval)));
            set(t,'Tag','cccstring2','Color',[1 1 0],'FontSize',storage_av3_stackbrowser.fontsize,'FontWeight','bold');
        end
    end
end

%draw quality values
delete(findobj('Tag','qualitystring'));
if storage_av3_stackbrowser.display.quality == 1
    for row=1:storage_av3_stackbrowser.numrows
        pos_y = (row-1).*storage_av3_stackbrowser.imsize(2)+storage_av3_stackbrowser.imsize(2).*0.3;
        abspartnumber = ((storage_av3_stackbrowser.page).*storage_av3_stackbrowser.numrows)+row;
        if abspartnumber <= size(storage_av3_stackbrowser.align,2)
            if storage_av3_stackbrowser.refstep == 0
                refstep = storage_av3_stackbrowser.refstep+1;
            else
                refstep = storage_av3_stackbrowser.refstep;
            end
            cccval = storage_av3_stackbrowser.align(refstep,abspartnumber).quality;
            pos_x = storage_av3_stackbrowser.imsize(1).*0.8;
            t = text(pos_x,pos_y,sprintf('%0.0f',double(cccval)));
            set(t,'Tag','cccstring2','Color',[1 1 0],'FontSize',storage_av3_stackbrowser.fontsize,'FontWeight','bold');
        end
    end
end

%draw class no
delete(findobj('Tag','qualitystring'));
if storage_av3_stackbrowser.display.classno == 1
    for row=1:storage_av3_stackbrowser.numrows
        pos_y = (row-1).*storage_av3_stackbrowser.imsize(2)+storage_av3_stackbrowser.imsize(2).*0.3;
        abspartnumber = ((storage_av3_stackbrowser.page).*storage_av3_stackbrowser.numrows)+row;
        if abspartnumber <= size(storage_av3_stackbrowser.align,2)
            if storage_av3_stackbrowser.refstep == 0
                refstep = storage_av3_stackbrowser.refstep+1;
            else
                refstep = storage_av3_stackbrowser.refstep;
            end
            cccval = storage_av3_stackbrowser.align(refstep,abspartnumber).Class;
            pos_x = storage_av3_stackbrowser.imsize(1).*1.2;
            t = text(pos_x,pos_y,sprintf('%0.0f',double(cccval)));
            set(t,'Tag','cccstring2','Color',[1 1 0],'FontSize',storage_av3_stackbrowser.fontsize,'FontWeight','bold');
        end
    end
end

%draw projection class no
delete(findobj('Tag','qualitystring'));
if storage_av3_stackbrowser.display.projclassno == 1
    for row=1:storage_av3_stackbrowser.numrows
        pos_y = (row-1).*storage_av3_stackbrowser.imsize(2)+storage_av3_stackbrowser.imsize(2).*0.9;
        abspartnumber = ((storage_av3_stackbrowser.page).*storage_av3_stackbrowser.numrows)+row;
        if abspartnumber <= size(storage_av3_stackbrowser.align,2)
            if storage_av3_stackbrowser.refstep == 0
                refstep = storage_av3_stackbrowser.refstep+1;
            else
                refstep = storage_av3_stackbrowser.refstep;
            end
            cccval = storage_av3_stackbrowser.align(refstep,abspartnumber).ProjectionClass;
            pos_x = storage_av3_stackbrowser.imsize(1).*1.2;
            t = text(pos_x,pos_y,sprintf('%0.0f',double(cccval)));
            set(t,'Tag','cccstring2','Color',[1 1 0],'FontSize',storage_av3_stackbrowser.fontsize,'FontWeight','bold');
        end
    end
end



%draw marks
delete(findobj('-regexp','Tag','rect_*'));
if storage_av3_stackbrowser.display.marks == 1
    for row=1:storage_av3_stackbrowser.numrows
        pos_y = (row-1).*storage_av3_stackbrowser.imsize(2)+storage_av3_stackbrowser.imsize(2).*0.1;
        abspartnumber = ((storage_av3_stackbrowser.page).*storage_av3_stackbrowser.numrows)+row;
        if abspartnumber <= size(storage_av3_stackbrowser.align,2)
            if storage_av3_stackbrowser.refstep == 0
                refstep = storage_av3_stackbrowser.refstep+1;
            else
                refstep = storage_av3_stackbrowser.refstep;
            end
            %mark particle green if selected and red if not selected
            if storage_av3_stackbrowser.align(1,abspartnumber).selected == 1
                color = [0 1 0];
            else
                color = [1 0 0];
            end
            pos_x = storage_av3_stackbrowser.imsize(1).*0.1;
            r = rectangle('Position',[pos_x,pos_y,storage_av3_stackbrowser.imsize(1).*0.2,storage_av3_stackbrowser.imsize(2).*0.2],'FaceColor',color,'Tag',['rect_' abspartnumber]);
        end
    end
end

set(tmpobj,'Tag','av3_stackbrowser_particles');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mark particle                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mark_particle(a,b)

global storage_av3_stackbrowser;

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
    row = ceil(y1./storage_av3_stackbrowser.imsize(2));
    abspartnumber = storage_av3_stackbrowser.page.*storage_av3_stackbrowser.numrows+(row);

if strcmp(button,'normal') == true

    
    if abspartnumber <= size(storage_av3_stackbrowser.align,2)
        %toggle selection on particle
        
        if storage_av3_stackbrowser.align(1,abspartnumber).selected == 1
            storage_av3_stackbrowser.align(1,abspartnumber).selected = 0;
            color = [1 0 0];
        else
            storage_av3_stackbrowser.align(1,abspartnumber).selected = 1;
            color = [0 1 0];
        end

        %redraw mark
        delete(findobj('Tag',['rect_' abspartnumber]));
        pos_y = (row-1).*storage_av3_stackbrowser.imsize(2)+storage_av3_stackbrowser.imsize(2).*0.1;
        pos_x = storage_av3_stackbrowser.imsize(1).*0.1;
        r = rectangle('Position',[pos_x,pos_y,storage_av3_stackbrowser.imsize(1).*0.2,storage_av3_stackbrowser.imsize(2).*0.2],'FaceColor',color,'Tag',['rect_' abspartnumber]);
    end

elseif strcmp(button,'alt') == true

     if storage_av3_stackbrowser.refstep == 0
         step =1;
     else
        step = storage_av3_stackbrowser.refstep;
     end
     string = ['Filename: ' storage_av3_stackbrowser.align(step,abspartnumber).Filename];
%     string = strvcat(string, ['Position: X ' num2str(storage_av2_stackbrowser.align(abspartnumber).position.x) ', Y ' num2str(storage_av2_stackbrowser.align(abspartnumber).position.y)]);
%     string = strvcat(string, ['Class: ' storage_av2_stackbrowser.align(abspartnumber).class]);
%     string = strvcat(string, ['Shift: X ' num2str(storage_av2_stackbrowser.align(abspartnumber).shift.x) ', Y ' num2str(storage_av2_stackbrowser.align(abspartnumber).shift.y)]);
%     string = strvcat(string, ['Angle: ' num2str(storage_av2_stackbrowser.align(abspartnumber).angle)]);
%     string = strvcat(string, ['Isaligned? ', num2str(storage_av2_stackbrowser.align(abspartnumber).isaligned)]);
%     string = strvcat(string, ['Dataset: ', storage_av2_stackbrowser.align(abspartnumber).dataset]);
     string = strvcat(string, ['Number in stack: ', num2str(abspartnumber)]);
     msgbox(string);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  sort 3d structure                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function structure = tom_sortstruct3d(structure, fieldname, order, step)

if step == 0
    step = 1;
end

if nargin < 3
    order = 'ascend';
end

sortmatrix = zeros(1,size(structure,2));

for i=1:size(structure,2)
    sortmatrix(i) = structure(step,i).(fieldname);
end

[sorted,idx] = sort(sortmatrix,2,order);

for i=1:size(structure,2)
    for j=1:size(structure,1)
        outstruct(j,i) = structure(j,idx(i)); 
    end
end

structure = outstruct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function menu_av3_stackbrowser_sortorder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function menu_av3_stackbrowser_sortfield_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function menu_av3_stackbrowser_norm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_av3_stackbrowser_filterval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function menu_av2_stackbrowser_filteroperator_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function menu_av2_stackbrowser_filterfield_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_av3_stackbrowser_bandpasshigh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_av3_stackbrowser_bandpasslow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_av3_stackbrowser_newquality_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_av3_stackbrowser_outalign_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_av3_stackbrowser_outdir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_av3_stackbrowser_histogramhigh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_av3_stackbrowser_histogramlow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider_av3_stackbrowser_particles_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listbox1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function menu_av2_stackbrowser_filterfield_Callback(hObject, eventdata, handles)
function menu_av2_stackbrowser_filteroperator_Callback(hObject, eventdata, handles)
function input_av3_stackbrowser_filterval_Callback(hObject, eventdata, handles)
function menu_av3_stackbrowser_sortfield_Callback(hObject, eventdata, handles)
function menu_av3_stackbrowser_sortorder_Callback(hObject, eventdata, handles)
function input_av3_stackbrowser_newquality_Callback(hObject, eventdata, handles)
function input_av3_stackbrowser_outdir_Callback(hObject, eventdata, handles)
function input_av3_stackbrowser_outalign_Callback(hObject, eventdata, handles)
function input_av3_stackbrowser_resamplez_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_av3_stackbrowser_resamplexy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function menu_av3_stackbrowser_orientation_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_av3_stackbrowser_histogramlow_Callback(hObject, eventdata, handles)
function edit_av3_stackbrowser_histogramhigh_Callback(hObject, eventdata, handles)
function menu_av3_stackbrowser_refstep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in constant_scale.
function constant_scale_Callback(hObject, eventdata, handles)
% hObject    handle to constant_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of constant_scale
